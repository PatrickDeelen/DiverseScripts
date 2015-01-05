package nl.systemsgenetics.genomestudioexporttoopticall;

import au.com.bytecode.opencsv.CSVReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.StringUtils;

/**
 * 
 *
 */
public class Main {

	private static final Pattern SAMPLE_COUNT_PATTERN = Pattern.compile("Total Samples\\s+(\\d+)");
	private static final Pattern SNP_COUNT_PATTERN = Pattern.compile("Total SNPs\\s+(\\d+)");

	private static final Options OPTIONS;

	static {

		OPTIONS = new Options();

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Final report file");
		OptionBuilder.withLongOpt("report");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('r'));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output folder");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('o'));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Probe update file, optional");
		OptionBuilder.withLongOpt("probe");
		OPTIONS.addOption(OptionBuilder.create('p'));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Samples to include, optional");
		OptionBuilder.withLongOpt("samples");
		OPTIONS.addOption(OptionBuilder.create('s'));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Variants to include, optional");
		OptionBuilder.withLongOpt("variants");
		OPTIONS.addOption(OptionBuilder.create('v'));
		
	}
	
	@SuppressWarnings("RedundantStringConstructorCall")
	public static void main(String[] args) throws Exception {

		int rsColumn = 0;
		int sampleColumn = 1;
		int chrColumn = 18;
		int chrPosColumn = 19;
		int allelesColumn = 22;
		int aSignalColumn = 29;
		int bSignalColumn = 30;
		
		final File finalReportFile;
		final File outputFolder;
		final File probeUpdateFile;
		final File sampleIncludeFile;
		final File variantIncludeFile;
		
		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			finalReportFile = new File(commandLine.getOptionValue('r'));
			outputFolder = new File(commandLine.getOptionValue('o'));
			
			if(commandLine.hasOption('p')){
				probeUpdateFile = new File(commandLine.getOptionValue('p'));
			} else {
				probeUpdateFile = null;
			}
			
			if(commandLine.hasOption('s')){
				sampleIncludeFile = new File(commandLine.getOptionValue('s'));
			} else {
				sampleIncludeFile = null;
			}
			
			if(commandLine.hasOption('v')){
				variantIncludeFile = new File(commandLine.getOptionValue('v'));
			} else {
				variantIncludeFile = null;
			}
			

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		System.out.println("Input file: " + finalReportFile.getAbsolutePath());
		System.out.println("Output folder: " + outputFolder.getAbsolutePath());

		if (probeUpdateFile != null) {
			System.out.println("Probe update file: " + probeUpdateFile.getAbsolutePath());
		}

		if (!outputFolder.isDirectory()) {
			if (!outputFolder.mkdirs()) {
				throw new Exception("Cannot create output dir");
			}
		}

		if (!outputFolder.canWrite()) {
			throw new Exception("Cannot write to output folder");
		}

		final HashMap<String, UpdatedProbeInfo> probeUpdates;
		if (probeUpdateFile != null) {
			probeUpdates = new HashMap<String, UpdatedProbeInfo>();

			CSVReader reader = new CSVReader(new FileReader(probeUpdateFile), '\t', '\0', 1);
			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {

				probeUpdates.put(nextLine[0], new UpdatedProbeInfo(nextLine[0], nextLine[1], nextLine[2], nextLine[3], nextLine[4], nextLine[5], nextLine[6]));

			}
		} else {
			probeUpdates = null;
		}
		
		final Set<String> sampleIncludeSet;
		if(sampleIncludeFile != null){
			sampleIncludeSet = readFilter(sampleIncludeFile);
		} else {
			sampleIncludeSet = null;
		}
		
		final Set<String> variantIncludeSet;
		if(variantIncludeFile != null){
			variantIncludeSet = readFilter(variantIncludeFile);
		} else {
			variantIncludeSet = null;
		}
		
		final HashSet<String> excludedSamplesByIncludeList = new HashSet<String>();
		final HashSet<String> excludedVariantsByIncludeList = new HashSet<String>();
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(finalReportFile), "UTF-8"));

		String line;

		line = reader.readLine();

		if (!line.equals("[Header]")) {
			throw new Exception("Expected: [Header]");
		}

		boolean inHeader = true;
		int sampleCount = 0;
		int snpCount = 0;
		int stringBuilderSize = 0;
		HashMap<String, TreeMap<Integer, StringBuilder>> outputData = new HashMap<String, TreeMap<Integer, StringBuilder>>(25);
		LinkedHashSet<String> samples = null;

		int excludedNaSnps = 0;
		int excludedInvalidSnps = 0;

		String currentSample = "";

		finalReportLines:
		while ((line = reader.readLine()) != null) {

			if (inHeader) {

				Matcher sampleCountMatcher = SAMPLE_COUNT_PATTERN.matcher(line);

				if (sampleCountMatcher.matches()) {
					sampleCount = Integer.parseInt(sampleCountMatcher.group(1));
					continue;
				}

				Matcher snpCountMatcher = SNP_COUNT_PATTERN.matcher(line);

				if (snpCountMatcher.matches()) {
					snpCount = Integer.parseInt(snpCountMatcher.group(1));
					continue;
				}

				if (line.equals("[Data]")) {

					if (sampleCount == 0) {
						throw new Exception("Did not find sample count in header");
					}
					if (snpCount == 0) {
						throw new Exception("Did not find SNP count in header");
					}

					samples = new LinkedHashSet<String>(sampleCount);

					stringBuilderSize = 100 + 12 * sampleCount;


					//skip header of data
					reader.readLine();

					inHeader = false;


				}


			} else {

				String[] elements = StringUtils.splitPreserveAllTokens(line, '\t');

				final String varId = elements[rsColumn];
				final String chr;
				final int pos;
				final char a;
				final char b;
				final String sampleId = elements[sampleColumn];
				
				if(sampleIncludeSet != null && !sampleIncludeSet.contains(sampleId)){
					excludedSamplesByIncludeList.add(sampleId);
					continue finalReportLines;
				}
				
				if(variantIncludeSet != null && !variantIncludeSet.contains(varId)){
					excludedVariantsByIncludeList.add(varId);
					continue finalReportLines;
				}

				if (probeUpdates == null) {

					String alleles = elements[allelesColumn];

					if (alleles.equals("[N/A]")) {
						++excludedNaSnps;
						continue finalReportLines;
					}

					if (alleles.length() != 5) {
						++excludedInvalidSnps;
						System.err.println("Illegal SNP alleles: " + alleles);
						continue finalReportLines;
					}

					chr = elements[chrColumn];
					pos = Integer.valueOf(elements[chrPosColumn]);
					a = alleles.charAt(1);
					b = alleles.charAt(3);
					
				} else {
					
					UpdatedProbeInfo probeUpdate = probeUpdates.get(varId);
					
					if(probeUpdate == null){
						continue;
					} else {
						chr = probeUpdate.getChr();
						pos = probeUpdate.getPos();
						a = probeUpdate.getaAllele();
						b = probeUpdate.getbAllele();
					}
					
				}



				if (!sampleId.equals(currentSample)) {
					if (samples.contains(sampleId)) {
						throw new Exception("Sample order problem");
					}
					samples.add(new String(sampleId));
					currentSample = sampleId;
				}

				TreeMap<Integer, StringBuilder> outputDataChr;
				if (outputData.containsKey(chr)) {
					outputDataChr = outputData.get(chr);
				} else {
					outputDataChr = new TreeMap<Integer, StringBuilder>();
					outputData.put(new String(chr), outputDataChr);
				}

				StringBuilder snpOutput;
				if (!outputDataChr.containsKey(pos)) {
					snpOutput = new StringBuilder(stringBuilderSize);
					snpOutput.append(varId);
					snpOutput.append('\t');
					snpOutput.append(String.valueOf(pos));
					snpOutput.append('\t');
					snpOutput.append(a);
					snpOutput.append(b);
					outputDataChr.put(pos, snpOutput);
				} else {
					snpOutput = outputDataChr.get(pos);
				}
				snpOutput.append('\t');
				snpOutput.append(elements[aSignalColumn]);
				snpOutput.append('\t');
				snpOutput.append(elements[bSignalColumn]);

			}

		}

		if(probeUpdates == null){
			System.out.println("Excluded N/A SNPs: " + excludedNaSnps);
			System.out.println("Excluded invalid alleles SNPs: " + excludedInvalidSnps);
		}
		
		if(sampleIncludeSet != null){
			System.out.println("Excluded samples by include list: " + excludedSamplesByIncludeList.size());
		}
		
		if(variantIncludeSet != null){
			System.out.println("Excluded variants by include list: " + excludedVariantsByIncludeList.size());
		}

		StringBuilder header = new StringBuilder("SNP\tCoor\tAlleles");
		for (String sample : samples) {
			header.append('\t');
			header.append(sample);
			header.append('A');
			header.append('\t');
			header.append(sample);
			header.append('B');
		}
		header.append('\n');

		String headerString = header.toString();

		int varCount = 0;
		for (Map.Entry<String, TreeMap<Integer, StringBuilder>> chrEntry : outputData.entrySet()) {

			String chr = chrEntry.getKey();
			TreeMap<Integer, StringBuilder> outputDataChr = chrEntry.getValue();

			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(outputFolder, "chr_" + chr + "_intensities.tsv")));

			writer.append(headerString);

			for (StringBuilder snpOutput : outputDataChr.values()) {
				snpOutput.append('\n');
				writer.append(snpOutput);
				++varCount;
			}

			writer.close();

		}
		
		System.out.println("Included probes: " + varCount);

	}
	
	private static Set<String> readFilter(File filterFile) throws UnsupportedEncodingException, IOException{
		
		HashSet<String> filter = new HashSet<String>();
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filterFile), "UTF-8"));
		
		String line;
		
		while ((line = reader.readLine()) != null) {
			filter.add(line);
		}
		
		return Collections.unmodifiableSet(filter);
		
	}
	
}
