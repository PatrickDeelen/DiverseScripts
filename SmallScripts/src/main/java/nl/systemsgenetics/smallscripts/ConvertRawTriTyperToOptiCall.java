package nl.systemsgenetics.smallscripts;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.GenotypeDataException;
import umcg.genetica.collections.ChrPosTreeMap;

/**
 *
 * @author Patrick Deelen
 */
public class ConvertRawTriTyperToOptiCall {

	public static void main(String[] args) throws Exception {

		String snpFile = "";
		String snpMapFile = "";
		String sampleFile = "";

		final TObjectIntHashMap<String> allSNPHash = new TObjectIntHashMap<String>();
		final ChrPosTreeMap<String> snps = new ChrPosTreeMap<String>();
		final ArrayList<String> samples = new ArrayList<String>();

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

		BufferedReader snpMapFileReader = new BufferedReader(new FileReader(snpMapFile));

		String[] chrPosId;
		while ((line = snpMapFileReader.readLine()) != null) {

			chrPosId = StringUtils.split(line, '\t');

			if (chrPosId.length != 3) {
				throw new GenotypeDataException("Error in Trityper SNPMappings.txt. Line does not contain 3 elements: ");
			}

			if (allSNPHash.containsKey(chrPosId[2])) {

				snps.put(chrPosId[0], Integer.valueOf(chrPosId[1]), new String(chrPosId[2]));

			}
		}

		snpMapFileReader.close();
		
		
		
		
		

	}
}
