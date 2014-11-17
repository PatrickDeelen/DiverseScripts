package nl.systemsgenetics.preparebealgevcf;

import gnu.trove.map.hash.TIntObjectHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import org.apache.commons.cli.ParseException;

/**
 * Hello world!
 *
 */
public class App {

	private static Configuration configuration;
	private static final String ENCODING = "ISO-8859-1";
	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private static final int VCF_CHR_COLUMN = 0;
	private static final int VCF_POS_COLUMN = 1;
	private static final int VCF_REF_COLUMN = 3;
	private static final int VCF_ALT_COLUMN = 4;
	private static final int VCF_GENOTYPE_COLUMN = 9;
	private static final int REF_CHR_COLUMN = 0;
	private static final int REF_POS_COLUMN = 1;
	private static final int REF_REF_COLUMN = 2;
	private static final int REF_ALT_COLUMN = 3;

	/**
	 * Filter single sample VCF based on reference SNP list for same alleles and
	 * presence Warm if ref and alt are swapped. If observed change this script
	 * so we can correct for it Keep track if imputation windows is not covered
	 * by any study SNPs, exclude reference SNPs in these windows
	 *
	 *
	 * @param args
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException, Exception {

		if (args.length == 0) {
			Configuration.printHelp();
			System.exit(1);
		}

		try {
			configuration = new Configuration(args);
		} catch (ParseException ex) {
			System.err.println("Configuration error: " + ex.getMessage());
			Configuration.printHelp();
			System.exit(1);
		}

		final File studyVcfFile = new File(configuration.getPathToInputVcf());
		final File referenceVariantListFile = new File(configuration.getReferenceVariantList());

		final File studyFilteredVcfFile = new File(configuration.getPathToOuputVcf());
		final File referenceExcludeFile = new File(configuration.getPathToReferenceExcludeList());

		final int chunkSize = configuration.getChunkSize();

		if (!studyVcfFile.canRead()) {
			throw new FileNotFoundException("Study VCF file not found at: " + studyVcfFile.getAbsolutePath());
		}

		if (!referenceVariantListFile.canRead()) {
			throw new FileNotFoundException("Reference variant file not found at: " + referenceVariantListFile.getAbsolutePath());
		}

		BufferedReader studyVcfReader;
		if (studyVcfFile.getName().endsWith(".gz")) {
			studyVcfReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(studyVcfFile)), ENCODING));
		} else {
			studyVcfReader = new BufferedReader(new InputStreamReader(new FileInputStream(studyVcfFile), ENCODING));
		}

		BufferedReader refernceVariantReader;
		if (referenceVariantListFile.getName().endsWith(".gz")) {
			refernceVariantReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(referenceVariantListFile)), ENCODING));
		} else {
			refernceVariantReader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceVariantListFile), ENCODING));
		}

		final File studyFilteredVcfFolder = studyFilteredVcfFile.getParentFile();
		if (studyFilteredVcfFolder != null) {
			if (!studyFilteredVcfFolder.isDirectory() && !studyFilteredVcfFolder.mkdirs()) {
				throw new IOException("Failed to create dir: " + studyFilteredVcfFolder.getAbsolutePath());
			}
		}

		final File referenceExcludeFolder = referenceExcludeFile.getParentFile();
		if (referenceExcludeFolder != null) {
			if (!referenceExcludeFolder.isDirectory() && !referenceExcludeFolder.mkdirs()) {
				throw new IOException("Failed to create dir: " + referenceExcludeFolder.getAbsolutePath());
			}
		}

		BufferedWriter studyVcfFilteredWriter;
		if (studyFilteredVcfFile.getName().endsWith(".gz")) {
			studyVcfFilteredWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(studyFilteredVcfFile))));
		} else {
			studyVcfFilteredWriter = new BufferedWriter(new FileWriter(studyFilteredVcfFile));
		}
		
		BufferedWriter referenceExcludeWriter = new BufferedWriter(new FileWriter(referenceExcludeFile));

		String vcfLine;
		String refLine = refernceVariantReader.readLine();
		String[] refLineElements;

		if (refLine == null) {
			throw new Exception("Ref file is empty");
		}

		refLineElements = TAB_PATTERN.split(refLine);
		int refLinePos = Integer.valueOf(refLineElements[REF_POS_COLUMN]);
		int refLineChunk = posToChunk(refLinePos, chunkSize);

		int lastVcfPos = 0;
		String lastChr = null;
		int excludeCount = 0;
		int previousChunkMatch = 0;

		TIntObjectHashMap<AtomicInteger> studyChunkIncluded = new TIntObjectHashMap<AtomicInteger>();
		TIntObjectHashMap<AtomicInteger> refChunkExcluded = new TIntObjectHashMap<AtomicInteger>();

		TIntObjectHashMap<HashSet<String>> potentialExclude = new TIntObjectHashMap<HashSet<String>>();

		{
			HashSet<String> chunkPotenialExclude = new HashSet<String>();
			potentialExclude.put(refLineChunk, chunkPotenialExclude);
			chunkPotenialExclude.add(refLineElements[REF_CHR_COLUMN] + ":" + refLineElements[REF_POS_COLUMN] + '\n');
		}

		while ((vcfLine = studyVcfReader.readLine()) != null) {
			if (vcfLine.length() == 0) {
				System.err.println("Found empty line in VCF. Skipping this line");
			} else if (vcfLine.charAt(0) == '#') {
				studyVcfFilteredWriter.write(vcfLine);
				studyVcfFilteredWriter.write('\n');
			} else {
				//Assess this variant line
				String[] vcfLineElements = TAB_PATTERN.split(vcfLine);

				String chr = vcfLineElements[VCF_CHR_COLUMN];
				int pos = Integer.parseInt(vcfLineElements[VCF_POS_COLUMN]);
				String ref = vcfLineElements[VCF_REF_COLUMN];
				String alt = vcfLineElements[VCF_ALT_COLUMN];
				int currentChunk = posToChunk(pos, chunkSize);

				if (lastVcfPos > pos) {
					throw new Exception("VCF not sorted. Filtering failed");
				}

				if (lastChr == null) {
					lastChr = chr;
				} else if (!lastChr.equals(chr)) {
					throw new Exception("Found multiple chr in VCF. This is not supported");
				}

				if (ref.length() != 1 || alt.length() != 1) {
					//only SNPs supported
					lastVcfPos = pos;
					//System.out.println("Excluding " + chr + ":" + vcfLineElements[VCF_POS_COLUMN] + " with non SNP alleles: " + ref + "/" + alt);
					++excludeCount;
					continue;
				}

				if(vcfLineElements[VCF_GENOTYPE_COLUMN].startsWith("./.")){
					lastVcfPos = pos;
					++excludeCount;
					continue;
				}





				while (refLinePos < pos) {
					refLine = refernceVariantReader.readLine();



					if (refLine != null) {

						refLineElements = TAB_PATTERN.split(refLine);

						if (!refLineElements[REF_CHR_COLUMN].equals(chr)) {
							throw new Exception("Error found chr: " + refLineElements[REF_CHR_COLUMN] + " in ref file, only expected: " + chr);
						}


						refLinePos = Integer.valueOf(refLineElements[REF_POS_COLUMN]);
						refLineChunk = posToChunk(refLinePos, chunkSize);

						if (refLineChunk > previousChunkMatch && refLineChunk < currentChunk) {
							referenceExcludeWriter.write(refLineElements[REF_CHR_COLUMN] + ":" + refLineElements[REF_POS_COLUMN] + '\n');

							if (!refChunkExcluded.containsKey(refLineChunk)) {
								refChunkExcluded.put(refLineChunk, new AtomicInteger(1));
							} else {
								refChunkExcluded.get(refLineChunk).incrementAndGet();
							}


						} else if (refLineChunk != previousChunkMatch) {
							HashSet<String> chunkPotenialExclude = potentialExclude.get(refLineChunk);
							if (chunkPotenialExclude == null) {
								chunkPotenialExclude = new HashSet<String>();
								potentialExclude.put(refLineChunk, chunkPotenialExclude);
							}
							chunkPotenialExclude.add(refLineElements[REF_CHR_COLUMN] + ":" + refLineElements[REF_POS_COLUMN] + '\n');
						}


					} else {
						refLineElements = null;
						refLinePos = -1;
						refLineChunk = -1;
						break;
					}

				}

				if (refLine != null && refLinePos == pos) {

					if (ref.equals(refLineElements[REF_REF_COLUMN]) && alt.equals(refLineElements[REF_ALT_COLUMN])) {
						studyVcfFilteredWriter.write(vcfLine);
						studyVcfFilteredWriter.write('\n');
						previousChunkMatch = currentChunk;

						if (!studyChunkIncluded.containsKey(currentChunk)) {
							studyChunkIncluded.put(currentChunk, new AtomicInteger(1));
						} else {
							studyChunkIncluded.get(currentChunk).incrementAndGet();
						}

						//Time to save the ref variants that are on the potential exclude list
						potentialExclude.remove(previousChunkMatch);

					} else if (alt.equals(refLineElements[REF_REF_COLUMN]) && ref.equals(refLineElements[REF_ALT_COLUMN])) {
						System.err.println("Warning REF and ALT are swapped. Variant is excluded since fixing this is not supported");
						++excludeCount;
					} else {
						//System.out.println("Excluding " + chr + ":" + vcfLineElements[VCF_POS_COLUMN] + " with alleles: " + ref + "/" + alt + " ref is " + refLineElements[REF_REF_COLUMN] + "/" + refLineElements[REF_ALT_COLUMN]);
						++excludeCount;
					}

				} else {
					//System.out.println("Excluding " + chr + ":" + vcfLineElements[VCF_POS_COLUMN] + " with alleles: " + ref + "/" + alt + " Not found in ref");
				}

				lastVcfPos = pos;


			}
		}


		//If there are still variant in the potential exclude file, exclude them now.
		for (int chunk : potentialExclude.keys()) {
			
			if (!refChunkExcluded.containsKey(chunk)) {
				refChunkExcluded.put(chunk, new AtomicInteger(1));
			}
			
			for (String excludeKey : potentialExclude.get(chunk)) {
				referenceExcludeWriter.write(excludeKey);//new line is already in this string
				refChunkExcluded.get(chunk).incrementAndGet();
			}
		}

		//Assess if there are ref variants remaining that need to be excluded
		while ((refLine = refernceVariantReader.readLine()) != null) {

			refLineElements = TAB_PATTERN.split(refLine);

			refLinePos = Integer.valueOf(refLineElements[REF_POS_COLUMN]);
			refLineChunk = posToChunk(refLinePos, chunkSize);

			if (refLineChunk > previousChunkMatch) {

				if (!refChunkExcluded.containsKey(refLineChunk)) {
					refChunkExcluded.put(refLineChunk, new AtomicInteger(1));
				} else {
					refChunkExcluded.get(refLineChunk).incrementAndGet();
				}

				referenceExcludeWriter.write(refLineElements[REF_CHR_COLUMN] + ":" + refLineElements[REF_POS_COLUMN] + '\n');
			}

		}
		
		System.out.println("Excluded total " + excludeCount + " variants from study");
		
//		System.out.println("Number of included study variants per chunk");
//		for(int chunk : studyChunkIncluded.keys()){
//			System.out.println(" - " + chunk + "\t" + studyChunkIncluded.get(chunk));
//		}
//		
//		System.out.println("Number of excluded reference variants per chunk");
//		for(int chunk : refChunkExcluded.keys()){
//			System.out.println(" - " + chunk + "\t" + refChunkExcluded.get(chunk));
//		}

		studyVcfFilteredWriter.close();
		studyVcfReader.close();
		refernceVariantReader.close();
		referenceExcludeWriter.close();

	}

	private static int posToChunk(int pos, int chunkSize) {

		return pos / chunkSize;

	}
}
