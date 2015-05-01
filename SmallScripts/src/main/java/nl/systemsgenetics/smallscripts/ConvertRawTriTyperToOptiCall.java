package nl.systemsgenetics.smallscripts;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.GenotypeDataException;
import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.io.trityper.WGAFileMatrixRawData;

/**
 *
 * @author Patrick Deelen
 */
public class ConvertRawTriTyperToOptiCall {

	private static final double twoDividedByPI = 2.0d / Math.PI;

	public static void main(String[] args) throws Exception {


		String trityperFolderPath = args[0];
		String snpInfoFilePath = args[1];
		String outputFolderPath = args[2];
		String sampleFilterFilePath = args.length > 3 ? args[3] : null;

		System.out.println("TriTyper: " + trityperFolderPath);
		System.out.println("SNP info: " + snpInfoFilePath);
		System.out.println("Output folder: " + outputFolderPath);
		System.out.println("Sample filter file: " + sampleFilterFilePath);

		File trityperFolder = new File(trityperFolderPath);
		File outputFolder = new File(outputFolderPath);
		outputFolder.mkdirs();

		File snpFile = new File(trityperFolder, "SNPs.txt");
		File sampleFile = new File(trityperFolder, "Individuals.txt");
		File rawdataFile = new File(trityperFolder, "RawDataMatrix.dat");

		final TObjectIntHashMap<String> allSNPHash = new TObjectIntHashMap<String>();
		final ChrPosTreeMap<String> snps = new ChrPosTreeMap<String>();
		final ArrayList<String> samples = new ArrayList<String>();
		final HashMap<String, String> snpAlleleMap = new HashMap<String, String>();

		BufferedReader snpFileReader = new BufferedReader(new FileReader(snpFile));

		String line;
		int i = 0;
		while ((line = snpFileReader.readLine()) != null) {
			allSNPHash.put(line, i);
			++i;
		}
		snpFileReader.close();

		final HashSet<String> includedSamples;
		if (sampleFilterFilePath == null) {
			includedSamples = null;
		} else {
			includedSamples = new HashSet<String>();
			File sampleFilterFile = new File(sampleFilterFilePath);

			BufferedReader sampleFilterFileReader = new BufferedReader(new FileReader(sampleFilterFile));
			while ((line = sampleFilterFileReader.readLine()) != null) {
				includedSamples.add(line);
			}
			sampleFilterFileReader.close();
		}


		BufferedReader sampleFileReader = new BufferedReader(new FileReader(sampleFile));
		while ((line = sampleFileReader.readLine()) != null) {
			samples.add(line);
		}
		sampleFileReader.close();

		final boolean[] sampleIncluded = new boolean[samples.size()];
		int samplesIncluded = 0;
		if (includedSamples == null) {
			//All are included
			for (int j = 0; j < sampleIncluded.length; ++j) {
				++samplesIncluded;
				sampleIncluded[j] = true;
			}
		} else {
			for (int j = 0; j < sampleIncluded.length; ++j) {
				boolean included = includedSamples.contains(samples.get(j));
				if(included){
					++samplesIncluded;
					sampleIncluded[j] = true;
				}
				
			}
		}

		System.out.println("Included samples: " + samplesIncluded);
		
		BufferedReader snpInfoFileReader = new BufferedReader(new FileReader(snpInfoFilePath));
		String[] snpInfo;
		while ((line = snpInfoFileReader.readLine()) != null) {

			snpInfo = StringUtils.split(line, '\t');

			if (snpInfo.length != 4) {
				throw new GenotypeDataException("Error in reading snp info. Line does not contain 4 elements: " + line);
			}

			if (allSNPHash.containsKey(snpInfo[2])) {

				@SuppressWarnings("RedundantStringConstructorCall")
				String snpId = new String(snpInfo[2]);
				snps.put(snpInfo[0], Integer.valueOf(snpInfo[1]), snpId);
				snpAlleleMap.put(snpId, snpInfo[3].intern());

			}
		}

		if (allSNPHash.size() != snps.size()) {
			System.err.println("Warning " + (allSNPHash.size() - snps.size()) + " snps without a mapping, these will be ignored");
		}

		snpInfoFileReader.close();

		WGAFileMatrixRawData rawdata = new WGAFileMatrixRawData(allSNPHash.size(), samples.size(), rawdataFile, true);



		for (String chr : snps.getChrs()) {

			BufferedWriter outputWriter = new BufferedWriter(new FileWriter(new File(outputFolder, "chr_" + chr + "_intensities.tsv")));

			outputWriter.append("SNP\tCoor\tAlleles");

			
			for (int sampleI = 0; sampleI < samples.size(); ++sampleI) {

				if (!sampleIncluded[sampleI]) {
					continue;
				}

				outputWriter.append('\t');
				outputWriter.append(samples.get(sampleI));
				outputWriter.append('A');
				outputWriter.append('\t');
				outputWriter.append(samples.get(sampleI));
				outputWriter.append('B');
				
				
			}

			outputWriter.write('\n');


			for (Integer pos : snps.getChrPositions(chr)) {

				String snp = snps.get(chr, pos);

				if (!allSNPHash.containsKey(snp)) {
					throw new Exception("SNP without data");
				}

				int snpI = allSNPHash.get(snp);

				outputWriter.write(snp);
				outputWriter.write('\t');
				outputWriter.write(String.valueOf(pos));
				outputWriter.write('\t');
				outputWriter.write(snpAlleleMap.get(snp));

				int countR0T0 = 0;
				
				for (int sampleI = 0; sampleI < samples.size(); ++sampleI) {

					if (!sampleIncluded[sampleI]) {
						continue;
					}

					byte rByte = rawdata.getR(snpI, sampleI);
					byte thetaByte = rawdata.getTheta(snpI, sampleI);
					
					if(rByte == 0 && thetaByte == 0){
						++countR0T0;
					}

					double rValue = (double) (-Byte.MIN_VALUE + (int) rByte) / 50d;
					double theta = (double) (-Byte.MIN_VALUE + (int) thetaByte) / 200d;
					double x = rValue / (1.0d + Math.tan(theta / twoDividedByPI));
					double y = rValue - x;

					outputWriter.write('\t');
					outputWriter.write(Double.toString(x));
					outputWriter.write('\t');
					outputWriter.write(Double.toString(y));

				}

				System.out.println(snp + "\t" + countR0T0);
				
				outputWriter.write('\n');

			}

			outputWriter.close();
		}
		
		System.out.println("Done");

	}
}
