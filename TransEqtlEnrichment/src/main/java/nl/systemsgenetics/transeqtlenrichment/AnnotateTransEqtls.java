/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.transeqtlenrichment;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author patri
 */
public class AnnotateTransEqtls {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		CSVReader reader = new CSVReaderBuilder(new InputStreamReader(new FileInputStream("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05.txt"))).withCSVParser(parser).withSkipLines(1).build();

		HashMap<String, ArrayList<String>> cellTypeToIdMapping = new HashMap<>();

		ArrayList<TransEqtl> transEqtls = new ArrayList<>();
		HashMap<String, TransEqtl> transEqtlsMap = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			String snp = nextLine[1];
			String transGene = nextLine[4];

			TransEqtl transEqtl = new TransEqtl(snp, transGene);

			transEqtls.add(transEqtl);
			transEqtlsMap.put(snp + "_" + transGene, transEqtl);

		}

		reader.close();

		ArrayList<TransEqtlExplanationSource> transEqtlExplanationSources = new ArrayList<>();

		transEqtlExplanationSources.add(new TransEqtlExplanationSource("TF_target", new File("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\explanations\\blood_networks_matrix.txt")));
		transEqtlExplanationSources.add(new TransEqtlExplanationSource("TF_predicted_target", new File("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\explanations\\blood_networks_matrix_expandenedTargets.txt")));
		transEqtlExplanationSources.add(new TransEqtlExplanationSource("Predicted_TF_target", new File("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\explanations\\blood_networks_matrix_expandenedTf.txt")));
		transEqtlExplanationSources.add(new TransEqtlExplanationSource("Predicted_TF_predicted_target", new File("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\explanations\\blood_networks_matrix_expandenedTargets_expandenedTf.txt")));
		transEqtlExplanationSources.add(new TransEqtlExplanationSource("Coregulation", new File("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\explanations\\coregulation.txt")));
		transEqtlExplanationSources.add(new TransEqtlExplanationSource("PPI", new File("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\explanations\\ppi.txt")));
		transEqtlExplanationSources.add(new TransEqtlExplanationSource("HiC_contact", new File("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\explanations\\hic.txt")));

		for (TransEqtlExplanationSource transEqtlExplanationSource : transEqtlExplanationSources) {

			reader = new CSVReaderBuilder(new InputStreamReader(new FileInputStream(transEqtlExplanationSource.getPath()))).withCSVParser(parser).withSkipLines(1).build();

			while ((nextLine = reader.readNext()) != null) {

				String snp = nextLine[0];
				String cisGene = nextLine[1];
				String transGene = nextLine[2];

				TransEqtl transEqtl = transEqtlsMap.get(snp + "_" + transGene);

				if (transEqtl != null) {

					transEqtl.addExplanation(new TransEqtlExplanation(cisGene, transEqtlExplanationSource.getName()));

				}

			}

			reader.close();

		}

		CSVWriter writer = new CSVWriter(new FileWriter(new File("C:\\UMCG\\Genetica\\Projects\\eQtlMeta\\2018-04-03-eQTLsFDR-Significant-0.05-WithoutCrossMappingEffects-PrunedFDR-Significant-0.05_explanations.txt")), '\t', '\0', '\0', "\n");
		String[] outputLine = new String[2 + transEqtlExplanationSources.size()];
		int c = 0;
		outputLine[c++] = "SNP";
		outputLine[c++] = "TransGene";
		for (TransEqtlExplanationSource transEqtlExplanationSource : transEqtlExplanationSources) {
			outputLine[c++] = transEqtlExplanationSource.getName();
		}
		writer.writeNext(outputLine);

		for (TransEqtl transEqtl : transEqtls) {
			c = 0;
			outputLine[c++] = transEqtl.getSnp();
			outputLine[c++] = transEqtl.getTransGene();
			

			HashMap<String, ArrayList<TransEqtlExplanation>> transEqtlExplanations = transEqtl.getTransEqtlExplanations();

			for (TransEqtlExplanationSource transEqtlExplanationSource : transEqtlExplanationSources) {
				ArrayList<TransEqtlExplanation> transEqtlExplanationsList = transEqtlExplanations.get(transEqtlExplanationSource.getName());
				StringBuilder explanations = new StringBuilder();
				if (transEqtlExplanationsList != null) {
					boolean first = true;
					for (TransEqtlExplanation transEqtlExplanation : transEqtlExplanationsList) {
						if (first) {
							first = false;
						} else {
							explanations.append(';');
						}
						explanations.append(transEqtlExplanation.getCisGene());
//						explanations.append('-');
//						explanations.append(transEqtlExplanation.getExplanation());
					}
				}
				outputLine[c++] = explanations.toString();
			}

			
			writer.writeNext(outputLine);
		}

		writer.close();

	}

}
