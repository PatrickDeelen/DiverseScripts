/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.smallscripts2;

import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import com.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;
import umcg.genetica.collections.intervaltree.StringRange;

/**
 *
 * @author patri
 */
public class QueryGsa {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		final File queryFile = new File("C:\\UMCG\\Genetica\\Projects\\GSA\\parkinsonGenes.txt");
		final File gsaFile = new File("C:\\UMCG\\Genetica\\Projects\\GSA\\GSAMD-24v1-0_20011747_A4.csv");
		final File outputFile = new File("C:\\UMCG\\Genetica\\Projects\\GSA\\parkinsonGenesGsa.txt");

		PerChrIntervalTree<StringRange> queryTreeMap = loadQuery(queryFile);

		final CSVParser parser = new CSVParserBuilder().withSeparator(',').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(gsaFile))).withSkipLines(8).withCSVParser(parser).build();

		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		int c = 0;
		String[] outputLine = new String[5];
		outputLine[c++] = "chr";
		outputLine[c++] = "pos";
		outputLine[c++] = "id";
		outputLine[c++] = "alleles";
		outputLine[c++] = "genes";
		writer.writeNext(outputLine);

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			//System.out.println(Arrays.toString(nextLine));
			
			if(nextLine.length < 20){
				continue;
			}
			
			String name = nextLine[1];
			String alleles = nextLine[3];
			String chr = nextLine[9];
			int pos = Integer.parseInt(nextLine[10]);

			
			
			List<StringRange> overlappingGenes = queryTreeMap.searchPosition(chr, pos);

			if (!overlappingGenes.isEmpty()) {

				StringBuilder genes = new StringBuilder();
				boolean first = true;
				for (StringRange gene : overlappingGenes) {

					if (first) {
						first = false;
					} else {
						genes.append(';');
					}
					genes.append(gene.getValue());

				}

				c = 0;
				outputLine[c++] = chr;
				outputLine[c++] = String.valueOf(pos);
				outputLine[c++] = name;
				outputLine[c++] = alleles;
				outputLine[c++] = genes.toString();
				writer.writeNext(outputLine);

			}

		}
		
		writer.close();

	}

	private static PerChrIntervalTree<StringRange> loadQuery(File queryFile) throws Exception {

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new BufferedReader(new FileReader(queryFile))).withSkipLines(0).withCSVParser(parser).build();

		PerChrIntervalTree<StringRange> queryTreeMap = new PerChrIntervalTree<>(StringRange.class);

		HashMap<String, ArrayList<StringRange>> chrGenes = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			String gene = nextLine[0];
			String chr = nextLine[1];
			int start = Integer.parseInt(nextLine[2]);
			int stop = Integer.parseInt(nextLine[3]);

			StringRange range = new StringRange(start, stop, gene);

			ArrayList<StringRange> genes = chrGenes.get(chr);
			if (genes == null) {
				genes = new ArrayList<>();
				chrGenes.put(chr, genes);
			}
			genes.add(range);

		}

		for (Map.Entry<String, ArrayList<StringRange>> chrGenesEntry : chrGenes.entrySet()) {
			queryTreeMap.addChrElements(chrGenesEntry.getKey(), chrGenesEntry.getValue());
		}

		return queryTreeMap;

	}

}
