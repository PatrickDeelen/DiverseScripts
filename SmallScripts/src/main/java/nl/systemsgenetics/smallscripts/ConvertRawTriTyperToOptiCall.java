package nl.systemsgenetics.smallscripts;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.util.ChrPos;
import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.io.trityper.WGAFileMatrixRawData;

/**
 *
 * @author Patrick Deelen
 */
public class ConvertRawTriTyperToOptiCall {

	private static final double twoDividedByPI = 2.0d / Math.PI;

	public static void main(String[] args) throws Exception {


		String trityperFolderPath = "";
		String snpInfoFilePath = "";
		String outputFolderPath = "";


		File trityperFolder = new File(trityperFolderPath);
		File outputFolder = new File(outputFolderPath);
		outputFolder.mkdirs();

		File snpFile = new File(trityperFolder, "SNPs.txt");
		File sampleFile = new File(trityperFolder, "Samples.txt");
		File rawdataFile = new File(trityperFolder, "***.txt");

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
				snpAlleleMap.put(snpId, snpInfo[4].intern());

			}
		}

		if (allSNPHash.size() != snps.size()) {
			System.err.println("Warning " + (allSNPHash.size() - snps.size()) + " snps without a mapping, these will be ignored");
		}

		snpInfoFileReader.close();

		WGAFileMatrixRawData rawdata = new WGAFileMatrixRawData(snps.size(), samples.size(), rawdataFile, true);



		for (String chr : snps.getChrs()) {

			BufferedWriter outputWriter = new BufferedWriter(new FileWriter(new File(outputFolder, "chr_" + chr + "_intensities.tsv")));

			outputWriter.append("SNP\tCoor\tAlleles");

			for (String sample : samples) {
				outputWriter.append('\t');
				outputWriter.append(sample);
			}

			outputWriter.write('\n');


			for (Integer pos : snps.getChrPositions(chr)) {

				String snp = snps.get(chr, pos);
				int snpI = allSNPHash.get(snp);

				outputWriter.write(snp);
				outputWriter.write('\t');
				outputWriter.write(String.valueOf(pos));
				outputWriter.write('\t');
				outputWriter.write(snpAlleleMap.get(snp));


				for (int sampleI = 0; sampleI < samples.size(); ++sampleI) {
					byte rByte = rawdata.getR(snpI, sampleI);
					byte thetaByte = rawdata.getTheta(snpI, sampleI);

					double rValue = (double) (-Byte.MIN_VALUE + (int) rByte) / 50d;
					double theta = (double) (-Byte.MIN_VALUE + (int) thetaByte) / 200d;
					double x = rValue / (1.0d + Math.tan(theta / twoDividedByPI));
					double y = rValue - x;

					outputWriter.write('\t');
					outputWriter.write(Double.toString(x));
					outputWriter.write('\t');
					outputWriter.write(Double.toString(y));

				}

				outputWriter.write('\n');

			}

			outputWriter.close();
		}

	}
}
