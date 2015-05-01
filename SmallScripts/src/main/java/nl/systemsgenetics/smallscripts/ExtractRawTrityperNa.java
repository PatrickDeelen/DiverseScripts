package nl.systemsgenetics.smallscripts;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.GenotypeDataException;
import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.io.trityper.WGAFileMatrixRawData;

/**
 *
 * @author Patrick Deelen
 */
public class ExtractRawTrityperNa {

	private static final double twoDividedByPI = 2.0d / Math.PI;

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		String trityperFolderPath = args[0];
		String snpInfoFilePath = args[1];
		String outputFolderPath = args[2];

		System.out.println("TriTyper: " + trityperFolderPath);
		System.out.println("SNP info: " + snpInfoFilePath);
		System.out.println("Output file: " + outputFolderPath);

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

		BufferedReader sampleFileReader = new BufferedReader(new FileReader(sampleFile));
		while ((line = sampleFileReader.readLine()) != null) {
			samples.add(line);
		}
		sampleFileReader.close();

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



		BufferedWriter snpOutputWriter = new BufferedWriter(new FileWriter(new File(outputFolder, "snpNaCounts")));
		BufferedWriter genotypeOutputWriter = new BufferedWriter(new FileWriter(new File(outputFolder, "genotypeNa")));
		BufferedWriter rOutputWriter = new BufferedWriter(new FileWriter(new File(outputFolder, "rNa")));
		BufferedWriter tOutputWriter = new BufferedWriter(new FileWriter(new File(outputFolder, "tNa")));

		snpOutputWriter.append("SNP\tCoor\tAlleles\tNaCount\n");

		TObjectIntHashMap<String> sampleNaCounts = new TObjectIntHashMap<String>();

		{
			String snp = "rs604864";
			BufferedWriter lookupWriter = new BufferedWriter(new FileWriter(new File(outputFolder, snp + ".txt")));
			int snpI = allSNPHash.get(snp);
			lookupWriter.append("sample\trByte\tthetaByte\tr\ttheta\ta\tb\n");
			for (int sampleI = 0; sampleI < samples.size(); ++sampleI) {

				byte rByte = rawdata.getR(snpI, sampleI);
				byte thetaByte = rawdata.getTheta(snpI, sampleI);

				double rValue = (double) (-Byte.MIN_VALUE + (int) rByte) / 50d;
				double theta = (double) (-Byte.MIN_VALUE + (int) thetaByte) / 200d;
				double x = rValue / (1.0d + Math.tan(theta / twoDividedByPI));
				double y = rValue - x;

				lookupWriter.append(samples.get(sampleI));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(rByte));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(thetaByte));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(rValue));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(theta));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(x));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(y));
				lookupWriter.append('\n');

			}
			lookupWriter.close();
		}
		
		{
			String snp = "rs4299462";
			BufferedWriter lookupWriter = new BufferedWriter(new FileWriter(new File(outputFolder, snp + ".txt")));
			int snpI = allSNPHash.get(snp);
			lookupWriter.append("sample\trByte\tthetaByte\tr\ttheta\ta\tb\n");
			for (int sampleI = 0; sampleI < samples.size(); ++sampleI) {

				byte rByte = rawdata.getR(snpI, sampleI);
				byte thetaByte = rawdata.getTheta(snpI, sampleI);

				double rValue = (double) (-Byte.MIN_VALUE + (int) rByte) / 50d;
				double theta = (double) (-Byte.MIN_VALUE + (int) thetaByte) / 200d;
				double x = rValue / (1.0d + Math.tan(theta / twoDividedByPI));
				double y = rValue - x;

				lookupWriter.append(samples.get(sampleI));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(rByte));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(thetaByte));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(rValue));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(theta));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(x));
				lookupWriter.append('\t');
				lookupWriter.append(String.valueOf(y));
				lookupWriter.append('\n');

			}
			lookupWriter.close();
		}

		for (String chr : snps.getChrs()) {
			for (Integer pos : snps.getChrPositions(chr)) {

				String snp = snps.get(chr, pos);

				if (!allSNPHash.containsKey(snp)) {
					throw new Exception("SNP without data");
				}

				int snpI = allSNPHash.get(snp);

				snpOutputWriter.write(snp);
				snpOutputWriter.write('\t');
				snpOutputWriter.write(String.valueOf(pos));
				snpOutputWriter.write('\t');
				snpOutputWriter.write(snpAlleleMap.get(snp));

				int countR0T0 = 0;

				for (int sampleI = 0; sampleI < samples.size(); ++sampleI) {

					byte rByte = rawdata.getR(snpI, sampleI);
					byte thetaByte = rawdata.getTheta(snpI, sampleI);

					if (rByte == 0 && thetaByte == 0) {
						++countR0T0;
						sampleNaCounts.adjustOrPutValue(samples.get(sampleI), 1, 1);
						genotypeOutputWriter.append(samples.get(sampleI));
						genotypeOutputWriter.append('\t');
						genotypeOutputWriter.append(snp);
						genotypeOutputWriter.append('\n');
					} else if (thetaByte == 0) {
						tOutputWriter.append(samples.get(sampleI));
						tOutputWriter.append('\t');
						tOutputWriter.append(snp);
						tOutputWriter.append('\n');
					} else if (rByte == 0) {
						rOutputWriter.append(samples.get(sampleI));
						rOutputWriter.append('\t');
						rOutputWriter.append(snp);
						rOutputWriter.append('\n');
					}

				}

				snpOutputWriter.write('\t');
				snpOutputWriter.write(String.valueOf(countR0T0));

				snpOutputWriter.write('\n');

			}
		}

		snpOutputWriter.close();
		genotypeOutputWriter.close();
		tOutputWriter.close();
		rOutputWriter.close();

		BufferedWriter sampleOutputWriter = new BufferedWriter(new FileWriter(new File(outputFolder, "sampleNaCounts")));
		for (String sample : samples) {
			sampleOutputWriter.write(sample);
			sampleOutputWriter.write('\t');
			sampleOutputWriter.write(String.valueOf(sampleNaCounts.get(sample)));
			sampleOutputWriter.write('\n');
		}
		sampleOutputWriter.close();

		System.out.println("Done");


	}
}
