/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.transeqtlenrichment;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.DoubleSorting;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import gnu.trove.map.hash.TDoubleDoubleHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;
import umcg.genetica.math.matrix2.DoubleMatrixDatasetFastSubsetLoader;
import umcg.genetica.math.stats.FisherExactTest;
import umcg.genetica.math.stats.ZScores;

/**
 *
 * @author patri
 */
public class TransEnrich {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, Exception {

		//convertTransPermutationFiles();
		//convertTransFiles();
		//prepareNetworks();
		//correlateNetworkToCoregulation();
		//tTestNetworkToCoregulation();
		//String prefix = "/groups/umcg-wijmenga/scr02/depict2/transEqtl/depict2Enrichment_v38/transEqtl";
		//createCoregulatedMatrix();
		String prefix = args[0];


		//String networkPath = "/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks_matrix";//_expandenedTargets_expandenedTf
		String networkPath = args[1];
		String overlapPath = args[2];

		System.out.println(prefix);
		System.out.println(networkPath);
		System.out.println(overlapPath);


		if (!new File(prefix + "_geneQvalues.dat").exists()) {
			System.out.println("Coverting p to q");
			convertPvalueToQvalue(prefix);
		}

		//expandTargetsUsingCoreg();
		//	expandTfUsingCoreg();
		DoubleMatrixDatasetFastSubsetLoader networkDataLoader = new DoubleMatrixDatasetFastSubsetLoader(networkPath);
//		DoubleMatrixDatasetFastSubsetLoader networkDataLoader = new DoubleMatrixDatasetFastSubsetLoader("/groups/umcg-wijmenga/scr02/depict2/transEqtl/InBio_Map_core_2016_09_12/gene-matrix_19-9-2019_expandenedTargets");

//		DoubleMatrixDatasetFastSubsetLoader networkDataLoader = new DoubleMatrixDatasetFastSubsetLoader("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks_matrix_expandenedTf");
		//	DoubleMatrixDatasetFastSubsetLoader networkDataReversedLoader = new DoubleMatrixDatasetFastSubsetLoader("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks_matrix_reversed");
		DoubleMatrixDatasetFastSubsetLoader genePvaluesLoader = new DoubleMatrixDatasetFastSubsetLoader(prefix + "_geneQvalues");
		//DoubleMatrixDatasetFastSubsetLoader geneCoregulationLoader = new DoubleMatrixDatasetFastSubsetLoader("/groups/umcg-wijmenga/scr02/depict2/transEqtl/depict2Enrichment_v37_2/transEqtl_CoregulationEqtlgen_Enrichment_zscoreExHla.txt");

		double minEvidence = 0;

		ArrayList<String> sharedTfGenes = new ArrayList<>();

		Map<String, Integer> genePvaluesHashGenes = genePvaluesLoader.getOriginalRowMap();
		for (String networkGene : networkDataLoader.getOriginalRowMap().keySet()) {
			if (genePvaluesHashGenes.containsKey(networkGene)) {
				sharedTfGenes.add(networkGene);
			}
		}

		System.out.println("Shared TF genes " + sharedTfGenes.size());

		DoubleMatrixDataset<String, String> networkData = networkDataLoader.loadSubsetOfRowsBinaryDoubleData(sharedTfGenes);

		DoubleMatrixDataset<String, String> genePvalues = genePvaluesLoader.loadSubsetOfRowsBinaryDoubleData(sharedTfGenes);

		HashMap<String, ArrayList<String>> geneSnpUsage = new HashMap<>();

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new InputStreamReader(new FileInputStream(prefix + "_usedVariantsPerGene.txt"))).withCSVParser(parser).build();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			String snp = nextLine[0];
			String gene = nextLine[1];

			ArrayList<String> geneSnps = geneSnpUsage.get(gene);
			if (geneSnps == null) {
				geneSnps = new ArrayList<>();
				geneSnpUsage.put(gene, geneSnps);
			}
			geneSnps.add(snp);

		}

		//DoubleMatrixDataset<String, String> genePvaluesP1 = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData("/groups/umcg-wijmenga/scr02/depict2/transEqtl/depict2Enrichment_v39/transEqtlCisAndWindowPermutation1_genePvalues", sharedTfGenes);
		ArrayList<String> sharedTargetGenes = networkData.getColObjects();
		Iterator<String> sharedColumnsIt = sharedTargetGenes.iterator();
		while (sharedColumnsIt.hasNext()) {
			String col = sharedColumnsIt.next();
			if (!genePvalues.getHashCols().containsKey(col)) {
				sharedColumnsIt.remove();
			}

		}

		networkData = networkData.viewColSelection(sharedTargetGenes);

		ArrayList<String> colGenes = networkData.getColObjects();
		HashSet<String> includedCols = new HashSet<>();
		for (int c = 0; c < networkData.columns(); ++c) {
			int countPassMin = 0;
			for (int r = 0; r < networkData.rows(); ++r) {
				if (networkData.getElementQuick(r, c) > minEvidence) {
					countPassMin++;
				}
			}
			//if here then no elements in row > minEvidence
			if (countPassMin > 0) {
				includedCols.add(colGenes.get(c));
			}
		}

		ArrayList<String> rowGenes = networkData.getRowObjects();
		HashSet<String> includedRows = new HashSet<>();
		for (int r = 0; r < networkData.rows(); ++r) {
			int countPassMin = 0;
			for (int c = 0; c < networkData.columns(); ++c) {
				if (networkData.getElementQuick(r, c) > minEvidence) {
					countPassMin++;
				}
			}
			//if here then no elements in row > minEvidence
			if (countPassMin > 0) {
				includedRows.add(rowGenes.get(r));
			}
		}

		networkData.viewSelection(rowGenes, colGenes);
		genePvalues = genePvalues.viewSelection(rowGenes, colGenes);
		//genePvaluesP1 = genePvaluesP1.viewSelection(rowGenes, colGenes);

		//System.out.println("name" + "\t" + "transEffectAndInNetwork" + "\t" + "transEffectNotInNetwork" + "\t" + "noTransEffectButInNetwork" + "\t" + "noTransEffectNotInNetwork" + "\t" + "Odds ratio" + "\t" + "Fisher p-value");
		//enrichmentRegulatoryCirc(genePvalues, networkData, "Real cis to trans", 0.05, minEvidence);
		enrichmentRegulatoryCirc2(genePvalues, networkData, "Real cis to trans", 0.05, minEvidence, geneSnpUsage, overlapPath);
		//enrichmentRegulatoryCirc(genePvalues, networkDataReversed, "Real trans to cis", 1e-4);

		//enrichmentRegulatoryCirc(genePvaluesP1, networkData, "Permutation cis to trans", 0.05, minEvidence);
		//enrichmentRegulatoryCirc(genePvaluesP1, networkDataReversed, "Permutation trans to cis", 1e-4);
		System.out.println("-------");
	}

	private static void convertPvalueToQvalue(String prefix) throws IOException, Exception {

		DoubleMatrixDataset<String, String> genePvaluesMatrix = DoubleMatrixDataset.loadDoubleBinaryData(prefix + "_genePvalues");

		ArrayList<String> allGenes = genePvaluesMatrix.getRowObjects();
		ArrayList<String> includedGenes = new ArrayList<>();
		for (int i = 0; i < genePvaluesMatrix.rows(); ++i) {
			if (!Double.isNaN(genePvaluesMatrix.getElementQuick(i, 0))) {
				includedGenes.add(allGenes.get(i));
			}
		}

		System.out.println("Number of genes: " + includedGenes.size());

		genePvaluesMatrix = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(prefix + "_genePvalues", includedGenes);

		DoubleMatrix1D genePvalues = genePvaluesMatrix.getMatrix().vectorize();
		DoubleMatrix1D genePvaluesP1 = DoubleMatrixDataset.loadSubsetOfRowsBinaryDoubleData(prefix + "Permutation1_genePvalues", includedGenes).getMatrix().vectorize();

		DoubleMatrix1D genePvaluesSorted = DoubleSorting.quickSort.sort(genePvalues);
		DoubleMatrix1D genePvaluesP1Sorted = DoubleSorting.quickSort.sort(genePvaluesP1);

		TDoubleDoubleHashMap pToqMap = new TDoubleDoubleHashMap();

		int numberPermutedSmallerThenCurrent = 0;

		for (int i = 0; i < genePvaluesSorted.size(); ++i) {

			double currentP = genePvaluesSorted.getQuick(i);

			while (numberPermutedSmallerThenCurrent < genePvaluesP1Sorted.size() && genePvaluesP1Sorted.getQuick(numberPermutedSmallerThenCurrent) <= currentP) {
				numberPermutedSmallerThenCurrent++;
			}

			//note i is current number of results
			double qvalue = (double) (numberPermutedSmallerThenCurrent) / (double) (i + 1);

			pToqMap.put(currentP, qvalue);

		}

		for (int r = 0; r < genePvaluesMatrix.rows(); ++r) {
			for (int c = 0; c < genePvaluesMatrix.columns(); ++c) {
				genePvaluesMatrix.setElementQuick(r, c, pToqMap.get(genePvaluesMatrix.getElementQuick(r, c)));
			}
		}

//		for(double pValue : pToqMap.keys()){
//			System.out.println(pValue + "\t" + pToqMap.get(pValue));
//		}
//		System.out.println("----");
		genePvaluesMatrix.saveBinary(prefix + "_geneQvalues");

	}

	private static void createCoregulatedMatrix() throws IOException {

		DoubleMatrixDataset<String, String> geneCoregulation = DoubleMatrixDataset.loadDoubleBinaryData("/groups/umcg-wijmenga/scr02/depict2/pathwayDatabases/geneCoregulationZscoreEqtlgen");

		System.out.println("Create coregulation matrix");

		double bonfP = 0.05 / (geneCoregulation.rows() * geneCoregulation.columns());
		double bonfZ = Math.abs(ZScores.pToZTwoTailed(bonfP));

		System.out.println("rows: " + geneCoregulation.rows() + " cols: " + geneCoregulation.columns() + " bon P: " + bonfP + " bonfZ: " + bonfZ);

		for (int i = 0; i < geneCoregulation.columns(); ++i) {
			for (int j = 0; j < geneCoregulation.columns(); ++j) {
				if (i == j) {
					continue;
				}
				if (Math.abs(geneCoregulation.getElementQuick(i, j)) >= bonfZ) {
					geneCoregulation.setElementQuick(i, j, 1);
				} else {
					geneCoregulation.setElementQuick(i, j, 0);
				}

			}
		}

		geneCoregulation.saveBinary("/groups/umcg-wijmenga/scr02/depict2/transEqtl/geneCoregulationBinary");

	}

	private static void expandTargetsUsingCoreg() throws IOException {

		DoubleMatrixDataset<String, String> networkData = DoubleMatrixDataset.loadDoubleBinaryData("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks_matrix");
		//DoubleMatrixDataset<String, String> networkData = DoubleMatrixDataset.loadDoubleBinaryData("/groups/umcg-wijmenga/scr02/depict2/transEqtl/InBio_Map_core_2016_09_12/gene-matrix_19-9-2019.dat");
		DoubleMatrixDataset<String, String> geneCoregulation = DoubleMatrixDataset.loadDoubleBinaryData("/groups/umcg-wijmenga/scr02/depict2/pathwayDatabases/geneCoregulationZscoreEqtlgen");

		HashSet<String> sharedTargetGenes = new HashSet<>();

		Map<String, Integer> genePvaluesHashGenes = geneCoregulation.getHashCols();
		for (String networkGene : networkData.getHashCols().keySet()) {
			if (genePvaluesHashGenes.containsKey(networkGene)) {
				sharedTargetGenes.add(networkGene);
			}
		}

		DoubleMatrixDataset<String, String> networkDataSharedCols = networkData.viewColSelection(sharedTargetGenes);
		DoubleMatrixDataset<String, String> geneCoregulationSharedCols = geneCoregulation.viewSelection(sharedTargetGenes, sharedTargetGenes);

		//for columns
		//for correlated colums
		//sum to original column
		for (int i = 0; i < networkDataSharedCols.columns(); ++i) {
			for (int j = 0; j < networkDataSharedCols.columns(); ++j) {
				if (i == j) {
					continue;
				}
				if (Math.abs(geneCoregulationSharedCols.getElementQuick(i, j)) >= 6.432250669562) {
					networkDataSharedCols.getCol(i).assign(networkDataSharedCols.getCol(j), cern.jet.math.tdouble.DoubleFunctions.max);
				}

			}
		}

		System.out.println("Done expanding targets");

		networkData.saveBinary("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks_matrix_expandenedTargets");
		//networkData.saveBinary("/groups/umcg-wijmenga/scr02/depict2/transEqtl/InBio_Map_core_2016_09_12/gene-matrix_19-9-2019_expandenedTargets");

		System.out.println("File saved");

	}

	private static void expandTfUsingCoreg() throws IOException, Exception {

		DoubleMatrixDatasetFastSubsetLoader networkDataLoader = new DoubleMatrixDatasetFastSubsetLoader("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks_matrix");
		DoubleMatrixDatasetFastSubsetLoader geneCoregulationLoader = new DoubleMatrixDatasetFastSubsetLoader("/groups/umcg-wijmenga/scr02/depict2/pathwayDatabases/geneCoregulationZscoreEqtlgen");

		HashSet<String> sharedTfGenes = new HashSet<>();

		Map<String, Integer> genePvaluesHashGenes = geneCoregulationLoader.getOriginalRowMap();
		for (String networkGene : networkDataLoader.getOriginalRowMap().keySet()) {
			if (genePvaluesHashGenes.containsKey(networkGene)) {
				sharedTfGenes.add(networkGene);
			}
		}

		DoubleMatrixDataset<String, String> networkData = networkDataLoader.loadSubsetOfRowsBinaryDoubleData(sharedTfGenes);

		DoubleMatrixDataset<String, String> coregulationData = geneCoregulationLoader.loadSubsetOfRowsBinaryDoubleData(sharedTfGenes);

		DoubleMatrixDataset<String, String> expandedNetworkData = new DoubleMatrixDataset<>(coregulationData.getHashCols(), networkData.getHashCols());

		HashSet<String> tfGenesExpanded = new HashSet<>();

		for (String tfGene : networkData.getRowObjects()) {

			//first put original new data
			//expandedNetworkData.getRow(tfGene).assign(networkData.getRow(tfGene));
			for (String coRegGene : coregulationData.getColObjects()) {
				if (sharedTfGenes.contains(coRegGene)) {
					continue;
				}
				if (Math.abs(coregulationData.getElement(tfGene, coRegGene)) >= 6.432250669562) {

					tfGenesExpanded.add(coRegGene);

					for (String targetGene : networkData.getColObjects()) {

						expandedNetworkData.setElement(coRegGene, targetGene, Math.max(expandedNetworkData.getElement(coRegGene, targetGene), networkData.getElement(tfGene, targetGene)));

					}

				}

			}

		}

		System.out.println("Updated number of TF genes: " + tfGenesExpanded.size());

		expandedNetworkData = expandedNetworkData.viewRowSelection(tfGenesExpanded);

		expandedNetworkData.saveBinary("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks_matrix_expandenedTf");

		System.out.println("Expanding complete");

	}

	private static void enrichmentRegulatoryCirc(DoubleMatrixDataset<String, String> genePvalues, DoubleMatrixDataset<String, String> networkData, String name, double genePvalue, double minEvidence) {

		ArrayList<String> colGenes = genePvalues.getColObjects();
		ArrayList<String> rowGenes = genePvalues.getRowObjects();

		FisherExactTest fet = new FisherExactTest();

		int transEffectAndInNetwork = 0;
		int transEffectNotInNetwork = 0;
		int noTransEffectButInNetwork = 0;
		int noTransEffectNotInNetwork = 0;

		for (int traitI = 0; traitI < genePvalues.columns(); ++traitI) {

			int transEffectAndInNetworkTrait = 0;
			int transEffectNotInNetworkTrait = 0;
			int noTransEffectButInNetworkTrait = 0;
			int noTransEffectNotInNetworkTrait = 0;

			for (int geneI = 0; geneI < genePvalues.rows(); ++geneI) {

				double pvalue = genePvalues.getElementQuick(geneI, traitI);
				double networkScore = networkData.getElementQuick(geneI, traitI);

				if (pvalue <= genePvalue) {
					if (networkScore > minEvidence) {

						//System.out.println(rowGenes.get(geneI) + " - " + colGenes.get(traitI));
						transEffectAndInNetwork++;
						transEffectAndInNetworkTrait++;
					} else {
						transEffectNotInNetwork++;
						transEffectNotInNetworkTrait++;
					}
				} else {
					if (networkScore > minEvidence) {
						noTransEffectButInNetwork++;
						noTransEffectButInNetworkTrait++;
					} else {
						noTransEffectNotInNetwork++;
						noTransEffectNotInNetworkTrait++;
					}
				}
			}

			//System.out.println(traits.get(traitI) + "\t" + transEffectAndInNetworkTrait + "\t" + transEffectNotInNetworkTrait + "\t" + noTransEffectButInNetworkTrait + "\t" + noTransEffectNotInNetworkTrait + "\t" + fet.getFisherPValue(transEffectAndInNetworkTrait, transEffectNotInNetworkTrait, noTransEffectButInNetworkTrait, noTransEffectNotInNetworkTrait));
		}

		double pvalue = fet.getFisherPValue(transEffectAndInNetwork, transEffectNotInNetwork, noTransEffectButInNetwork, noTransEffectNotInNetwork);

		double ratio = (noTransEffectNotInNetwork / (double) noTransEffectButInNetwork) / (transEffectNotInNetwork / (double) transEffectAndInNetwork);

		//System.out.println(name + "\t" + transEffectAndInNetwork + "\t" + transEffectNotInNetwork + "\t" + noTransEffectButInNetwork + "\t" + noTransEffectNotInNetwork + "\t" + ratio + "\t" + pvalue);
		System.out.println(transEffectAndInNetwork + "\t" + noTransEffectButInNetwork + "\n" + transEffectNotInNetwork + "\t" + noTransEffectNotInNetwork + "\n" + ratio + "\t" + pvalue);

	}

	private static void enrichmentRegulatoryCirc2(DoubleMatrixDataset<String, String> genePvalues, DoubleMatrixDataset<String, String> networkData, String name, double genePvalue, double minEvidence, HashMap<String, ArrayList<String>> geneSnpUsage, String overlapPath) throws IOException {

		ArrayList<String> colGenes = genePvalues.getColObjects();
		ArrayList<String> rowGenes = genePvalues.getRowObjects();

		FisherExactTest fet = new FisherExactTest();

		int transEffectAndInNetwork = 0;
		int transEffectNotInNetwork = 0;
		int noTransEffectButInNetwork = 0;
		int noTransEffectNotInNetwork = 0;

		CSVWriter writer = new CSVWriter(new FileWriter(new File(overlapPath)), '\t', '\0', '\0', "\n");
		String[] outputLine = new String[3];
		int c = 0;
		outputLine[c++] = "SNP";
		outputLine[c++] = "CisGene";
		outputLine[c++] = "TransGene";
		writer.writeNext(outputLine);

		for (int traitI = 0; traitI < genePvalues.columns(); ++traitI) {

			int transEffectAndInNetworkTrait = 0;
			int transEffectNotInNetworkTrait = 0;
			int noTransEffectButInNetworkTrait = 0;
			int noTransEffectNotInNetworkTrait = 0;

			for (int geneI = 0; geneI < genePvalues.rows(); ++geneI) {

				double pvalue = genePvalues.getElementQuick(geneI, traitI);
				double networkScore = networkData.getElementQuick(geneI, traitI);

				if (pvalue <= genePvalue) {
					if (networkScore > minEvidence) {

						//System.out.println(rowGenes.get(geneI) + " - " + colGenes.get(traitI));
						String cisGene = rowGenes.get(geneI);

						ArrayList<String> snpsLinkedToGene = geneSnpUsage.get(cisGene);
						if (snpsLinkedToGene != null) {
							for (String snp : snpsLinkedToGene) {
								c = 0;
								outputLine[c++] = snp;
								outputLine[c++] = cisGene;
								outputLine[c++] = colGenes.get(traitI);
								writer.writeNext(outputLine);
							}
						}

						transEffectAndInNetwork++;
						transEffectAndInNetworkTrait++;
					} else {
						transEffectNotInNetwork++;
						transEffectNotInNetworkTrait++;
					}
				} else {
					if (networkScore > minEvidence) {
						noTransEffectButInNetwork++;
						noTransEffectButInNetworkTrait++;
					} else {
						noTransEffectNotInNetwork++;
						noTransEffectNotInNetworkTrait++;
					}
				}
			}

			//System.out.println(traits.get(traitI) + "\t" + transEffectAndInNetworkTrait + "\t" + transEffectNotInNetworkTrait + "\t" + noTransEffectButInNetworkTrait + "\t" + noTransEffectNotInNetworkTrait + "\t" + fet.getFisherPValue(transEffectAndInNetworkTrait, transEffectNotInNetworkTrait, noTransEffectButInNetworkTrait, noTransEffectNotInNetworkTrait));
		}
		
		writer.close();

		double pvalue = fet.getFisherPValue(transEffectAndInNetwork, transEffectNotInNetwork, noTransEffectButInNetwork, noTransEffectNotInNetwork);

		double ratio = (noTransEffectNotInNetwork / (double) noTransEffectButInNetwork) / (transEffectNotInNetwork / (double) transEffectAndInNetwork);

		//System.out.println(name + "\t" + transEffectAndInNetwork + "\t" + transEffectNotInNetwork + "\t" + noTransEffectButInNetwork + "\t" + noTransEffectNotInNetwork + "\t" + ratio + "\t" + pvalue);
		System.out.println(transEffectAndInNetwork + "\t" + noTransEffectButInNetwork + "\n" + transEffectNotInNetwork + "\t" + noTransEffectNotInNetwork + "\n" + ratio + "\t" + pvalue);

	}

	private static void correlateNetworkToCoregulation() throws IOException {

		//
		DoubleMatrixDataset<String, String> networkData = DoubleMatrixDataset.loadDoubleBinaryData("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks_matrix");
		DoubleMatrixDataset<String, String> coregulationData = DoubleMatrixDataset.loadDoubleBinaryData("/groups/umcg-wijmenga/scr02/depict2/transEqtl/depict2Enrichment_v38/transEqtl_CoregulationEqtlgen_Enrichment_zscoreExHla");

		System.out.println("Network rows: " + networkData.rows() + " cols: " + networkData.columns());
		System.out.println("Coregulation rows: " + coregulationData.rows() + " cols: " + coregulationData.columns());

		ArrayList<String> sharedColumns = networkData.getColObjects();
		Iterator<String> sharedColumnsIt = sharedColumns.iterator();
		while (sharedColumnsIt.hasNext()) {
			String col = sharedColumnsIt.next();
			if (!coregulationData.getHashCols().containsKey(col)) {
				sharedColumnsIt.remove();
			}

		}

		networkData = networkData.viewColSelection(sharedColumns);
		coregulationData = coregulationData.viewColSelection(sharedColumns);

		final SimpleRegression regression = new SimpleRegression();

		ArrayList<String> genes = networkData.getColObjects();

		for (String r : networkData.getRowObjects()) {

			regression.clear();

			if (!coregulationData.containsRow(r)) {
				continue;
			}

			DoubleMatrix1D netWorkDataRow = networkData.getRow(r);
			DoubleMatrix1D coregulationDataRow = coregulationData.getRow(r);

			for (int c = 0; c < networkData.columns(); ++c) {

//				if (Double.isNaN(networkData.getElementQuick(c, c))) {
//					System.err.println("Nan in network data");
//				}
//
//				if (Double.isNaN(coregulationData.getElementQuick(c, c))) {
//					System.err.println("Nan in corregulation data");
//				}
				//System.out.println(networkData.getElementQuick(r, c) + "\t" +  coregulationData.getElementQuick(r, c));
				regression.addData(netWorkDataRow.getQuick(c), coregulationDataRow.getQuick(c));

			}

			System.out.println(r + "\t" + regression.getR() + "\t" + regression.getN());

		}
	}

	private static void tTestNetworkToCoregulation() throws IOException {

		//
		DoubleMatrixDataset<String, String> networkData = DoubleMatrixDataset.loadDoubleBinaryData("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks_matrix");
		//DoubleMatrixDataset<String, String> coregulationData = DoubleMatrixDataset.loadDoubleBinaryData("/groups/umcg-wijmenga/scr02/depict2/transEqtl/depict2Enrichment_v38/transEqtl_CoregulationEqtlgen_Enrichment_zscoreExHla");
		DoubleMatrixDataset<String, String> coregulationData = DoubleMatrixDataset.loadDoubleBinaryData("/groups/umcg-wijmenga/scr02/depict2/transEqtl/depict2Enrichment_v38/transEqtlPermutation1_CoregulationEqtlgen_Enrichment_zscoreExHla");

		System.out.println("Network rows: " + networkData.rows() + " cols: " + networkData.columns());
		System.out.println("Coregulation rows: " + coregulationData.rows() + " cols: " + coregulationData.columns());

		ArrayList<String> sharedColumns = networkData.getColObjects();
		Iterator<String> sharedColumnsIt = sharedColumns.iterator();
		while (sharedColumnsIt.hasNext()) {
			String col = sharedColumnsIt.next();
			if (!coregulationData.getHashCols().containsKey(col)) {
				sharedColumnsIt.remove();
			}

		}

		networkData = networkData.viewColSelection(sharedColumns);
		coregulationData = coregulationData.viewColSelection(sharedColumns);

		final SimpleRegression regression = new SimpleRegression();

		ArrayList<String> genes = networkData.getColObjects();

		for (String r : networkData.getRowObjects()) {

			regression.clear();

			if (!coregulationData.containsRow(r)) {
				continue;
			}

			DoubleMatrix1D netWorkDataRow = networkData.getRow(r);
			DoubleMatrix1D coregulationDataRow = coregulationData.getRow(r);

			int genesInNetworkRow = netWorkDataRow.cardinality();
			int genesNotInNetworkRow = (int) netWorkDataRow.size() - genesInNetworkRow;

			double[] coregulationZscoreGenesInNetwork = new double[genesInNetworkRow];
			double[] coregulationZscoreGenesNotInNetwork = new double[genesNotInNetworkRow];

			int a = 0;
			int b = 0;

			for (int c = 0; c < networkData.columns(); ++c) {

				if (netWorkDataRow.getQuick(c) > 0) {
					coregulationZscoreGenesInNetwork[a++] = coregulationDataRow.getQuick(c);
				} else {
					coregulationZscoreGenesNotInNetwork[b++] = coregulationDataRow.getQuick(c);
				}

			}

			TTest ttest = new TTest();
			double pvalue = ttest.homoscedasticTTest(coregulationZscoreGenesInNetwork, coregulationZscoreGenesNotInNetwork);

			System.out.println(r + "\t" + pvalue);

		}
	}

	private static void convertTransPermutationFiles() throws IOException, NumberFormatException, Exception {

		DoubleMatrixDatasetFastSubsetLoader realTrans = new DoubleMatrixDatasetFastSubsetLoader("/groups/umcg-wijmenga/scr02/depict2/transEqtl/ZScoreMatrix.txt.binary");

		Set<String> snps = realTrans.getOriginalRowMap().keySet();
		Set<String> genes = realTrans.getOriginalColMap().keySet();

		HashSet<String> colToExclude = new HashSet<>();
		colToExclude.add("Alleles");
		colToExclude.add("AlleleAssessed");

		DoubleMatrixDataset<String, String> permutedMatrix = DoubleMatrixDataset.loadDoubleTextDoubleDataExlcudeCols("/groups/umcg-wijmenga/scr01/depict2/transData/2018-01-24-CombinedMetaAnalysis/ZScoreMatrix-Permutation1.txt.gz", '\t', colToExclude);

		LinkedHashMap<String, Integer> newHashCol = new LinkedHashMap<>();

		for (Map.Entry<String, Integer> oldEntry : permutedMatrix.getHashCols().entrySet()) {
			String oldGene = oldEntry.getKey();
			String[] oldGeneElements = oldGene.split("_");
			newHashCol.put(oldGeneElements[0], oldEntry.getValue());

		}

		permutedMatrix.setHashCols(newHashCol);

		permutedMatrix = permutedMatrix.viewSelection(snps, genes);

		for (int r = 0; r < permutedMatrix.rows(); ++r) {
			for (int c = 0; c < permutedMatrix.columns(); ++c) {
				if (Double.isNaN(permutedMatrix.getElementQuick(r, c))) {
					permutedMatrix.setElementQuick(r, c, 0);
				}
			}
		}

		permutedMatrix.saveBinary("/groups/umcg-wijmenga/scr01/depict2/transData/2018-01-24-CombinedMetaAnalysis/selectedEffects/ZScoreMatrix-Permutation1SubsetNaIs0");

	}

	private static void convertTransFiles() throws IOException, NumberFormatException, Exception {

		HashSet<String> colToExclude = new HashSet<>();
		colToExclude.add("Alleles");
		colToExclude.add("AlleleAssessed");

		DoubleMatrixDataset<String, String> transNumberOfSamples = DoubleMatrixDataset.loadDoubleTextDoubleDataExlcudeCols("/groups/umcg-wijmenga/scr01/depict2/transData/2018-01-24-CombinedMetaAnalysis/ZScoreMatrixNrSamples.txt.gz", '\t', colToExclude);

		ArrayList<String> includedSnpsInitially = new ArrayList<>();

		for (String snp : transNumberOfSamples.getHashRows().keySet()) {

			DoubleMatrix1D rowCounts = transNumberOfSamples.getRow(snp);

			int countLt10000 = 0;

			for (int i = 0; i < rowCounts.size(); ++i) {
				if (rowCounts.getQuick(i) < 10000) {
					countLt10000++;
				}
			}

			if (countLt10000 <= 10000) {
				includedSnpsInitially.add(snp);
			}

			System.out.println(countLt10000);

		}

		DoubleMatrixDataset<String, String> transNumberOfSamples2 = transNumberOfSamples.viewRowSelection(includedSnpsInitially);

		ArrayList<String> includedGenes = new ArrayList<>();

		for (String gene : transNumberOfSamples2.getHashCols().keySet()) {

			DoubleMatrix1D colCounts = transNumberOfSamples2.getCol(gene);

			int countLt10000 = 0;

			for (int i = 0; i < colCounts.size(); ++i) {
				if (colCounts.getQuick(i) < 5000) {
					countLt10000++;
				}
			}

			if (countLt10000 <= 50) {
				includedGenes.add(gene);
			}

		}

		DoubleMatrixDataset<String, String> transNumberOfSamplesSelectedGenes = transNumberOfSamples2.viewColSelection(includedGenes);

		ArrayList<String> includedSnps = new ArrayList<>();

		for (String snp : transNumberOfSamplesSelectedGenes.getHashRows().keySet()) {

			DoubleMatrix1D rowCounts = transNumberOfSamplesSelectedGenes.getRow(snp);

			int countLt10000 = 0;

			for (int i = 0; i < rowCounts.size(); ++i) {
				if (rowCounts.getQuick(i) < 5000) {
					countLt10000++;
				}
			}

			if (countLt10000 == 0) {
				includedSnps.add(snp);
			}

			System.out.println(countLt10000);

		}

		System.out.println("Including initially " + includedSnpsInitially.size() + " out of " + transNumberOfSamples.rows() + " SNPs");
		System.out.println("Including " + includedGenes.size() + " out of " + transNumberOfSamples.columns() + " genes");
		System.out.println("Including " + includedSnps.size() + " out of " + transNumberOfSamples.rows() + " SNPs");

		DoubleMatrixDataset<String, String> realMatrix = DoubleMatrixDataset.loadDoubleTextData("/groups/umcg-wijmenga/scr01/depict2/transData/2018-01-24-CombinedMetaAnalysis/ZScoreMatrix.txt.gz", '\t').viewSelection(includedSnps, includedGenes);
		realMatrix.saveBinary("/groups/umcg-wijmenga/scr01/depict2/transData/2018-01-24-CombinedMetaAnalysis/selectedEffects/ZScoreMatrixSubset");

	}

	private static void prepareNetworks() throws IOException, NumberFormatException {
		HashMap<String, String> geneMapping = readGenesMapping(new File("/groups/umcg-wijmenga/scr02/depict2/ensgR75.txt"));

		//File networkFolder = new File("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks/");
		File networkFolder = new File("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks/");

		ArrayList<File> networkFiles = new ArrayList<>();

		for (File file : networkFolder.listFiles()) {
			if (file.getName().endsWith(".gz")) {
				networkFiles.add(file);
			}
		}

		DoubleMatrixDataset<String, String> networkInteraction = new DoubleMatrixDataset<>(geneMapping.values(), geneMapping.values());
		DoubleMatrixDataset<String, String> networkInteractionReversed = new DoubleMatrixDataset<>(geneMapping.values(), geneMapping.values());

		HashSet<String> nonZeroRow = new HashSet<>();
		HashSet<String> nonZeroCol = new HashSet<>();

		LinkedHashMap<String, Integer> networkInteractionHash = networkInteraction.getHashRows();

		for (File networkFile : networkFiles) {

			System.out.println("Parsing file: " + networkFile.getName());

			final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
			final CSVReader reader = new CSVReaderBuilder(new InputStreamReader(new GZIPInputStream(new FileInputStream(networkFile)))).withCSVParser(parser).build();

			String[] nextLine;
			while ((nextLine = reader.readNext()) != null) {

				if (nextLine.length != 3) {
					System.out.println("Error on line: " + String.join("\t", nextLine));
				}

				String gene1 = nextLine[0];
				String gene2 = nextLine[1];

				String gene1ensg = geneMapping.get(gene1);
				String gene2ensg = geneMapping.get(gene2);

				if (gene1ensg == null || gene2ensg == null) {
					continue;
				}

				double score = Double.parseDouble(nextLine[2]);

				int gene1e = networkInteractionHash.get(gene1ensg);
				int gene2e = networkInteractionHash.get(gene2ensg);

				nonZeroRow.add(gene1ensg);
				nonZeroCol.add(gene2ensg);

				networkInteraction.setElementQuick(gene1e, gene2e, score + networkInteraction.getElementQuick(gene1e, gene2e));
				networkInteractionReversed.setElementQuick(gene2e, gene1e, score + networkInteraction.getElementQuick(gene2e, gene1e));

			}
		}

		networkInteraction = networkInteraction.viewSelection(nonZeroRow, nonZeroCol);
		networkInteractionReversed = networkInteractionReversed.viewSelection(nonZeroCol, nonZeroRow);

		networkInteraction.saveBinary("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks_matrix");
		networkInteractionReversed.saveBinary("/groups/umcg-wijmenga/scr02/depict2/regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/blood_networks_matrix_reversed");

	}

	private static HashMap<String, String> readGenesMapping(File geneFile) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(geneFile))).withCSVParser(parser).withSkipLines(1).build();

		HashMap<String, String> geneMapping = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			geneMapping.put(nextLine[6], nextLine[0]);

		}

		return geneMapping;

	}

	protected static final List<String> readMatrixAnnotations(File file) throws IOException {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(file))).withCSVParser(parser).build();

		ArrayList<String> identifiers = new ArrayList<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			identifiers.add(nextLine[0]);

		}

		return identifiers;

	}

}
