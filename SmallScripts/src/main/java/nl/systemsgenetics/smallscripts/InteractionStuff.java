package nl.systemsgenetics.smallscripts;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedHashSet;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author Patrick Deelen
 */
public class InteractionStuff {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, Exception {

		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream("D:\\UMCG\\Genetica\\Projects\\rp3\\oppositeEffects\\opposite_interactions_metaZ.txt"), "UTF-8"));

		String line;

		TObjectDoubleHashMap<String> covariateMaxAbsZ = new TObjectDoubleHashMap<String>(50000, 0.75f, 0);
		reader.readLine();
		while ((line = reader.readLine()) != null) {
			String[] elements = StringUtils.split(line, '\t');

			final String covarivate = elements[2];
			final double zscore = Math.abs(Double.parseDouble(elements[3]));

			if (zscore > covariateMaxAbsZ.get(covarivate)) {
				covariateMaxAbsZ.put(covarivate, zscore);
			}

		}

		reader.close();
		
		System.out.println("Total covariats: " + covariateMaxAbsZ.size());


		LinkedHashSet<String> covariateIncluded = new LinkedHashSet<String>();
		for (String covaraite : covariateMaxAbsZ.keySet()) {
			if (covariateMaxAbsZ.get(covaraite) >= 11) {
				covariateIncluded.add(covaraite);
			}
		}

		System.out.println("Total covariats passed: " + covariateIncluded.size());

		reader = new BufferedReader(new InputStreamReader(new FileInputStream("D:\\UMCG\\Genetica\\Projects\\rp3\\oppositeEffects\\opposite_interactions_metaZ.txt"), "UTF-8"));

		reader.readLine();
		HashMap<String, TObjectDoubleHashMap<String>> geneCovariatesMetaZ = new HashMap<String, TObjectDoubleHashMap<String>>();
		HashMap<String, String> geneVariant = new HashMap<String, String>();
		while ((line = reader.readLine()) != null) {
			String[] elements = StringUtils.split(line, '\t');

			final String variant = elements[0];
			final String gene = elements[1];
			final String covarivate = elements[2];
			final double zscore = Double.parseDouble(elements[3]);

			TObjectDoubleHashMap<String> metaZ = geneCovariatesMetaZ.get(gene);
			if (metaZ == null) {
				metaZ = new TObjectDoubleHashMap<String>(50000, 0.75f, Double.NaN);
				geneCovariatesMetaZ.put(gene, metaZ);
			}

			metaZ.put(covarivate, zscore);

			String variantCheck = geneVariant.get(gene);
			if (variantCheck == null) {
				geneVariant.put(gene, variant);
			} else if (!variantCheck.equals(variant)) {
				throw new Exception();
			}


		}
		
		System.out.println("Genes: " + geneVariant.size());

		BufferedWriter writer = new BufferedWriter(new FileWriter("D:\\UMCG\\Genetica\\Projects\\rp3\\oppositeEffects\\opposite_interactions_metaZ_topMatrix.txt"));

		writer.write("Variant\tGene");
		for (String covariate : covariateIncluded) {
			writer.write('\t');
			writer.write(covariate);
		}
		writer.write('\n');

		for (String gene : geneVariant.keySet()) {
			String variant = geneVariant.get(gene);
			TObjectDoubleHashMap<String> metaZ = geneCovariatesMetaZ.get(gene);
			writer.write(variant);
			writer.write('\t');
			writer.write(gene);
			for (String covariate : covariateIncluded) {
				writer.write('\t');
				writer.write(String.valueOf(metaZ.get(covariate)));
			}
			writer.write('\n');

			
		}

		writer.close();
	}
}
