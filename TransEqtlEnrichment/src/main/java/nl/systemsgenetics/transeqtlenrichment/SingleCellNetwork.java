/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.transeqtlenrichment;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.IntStream;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.FastMath;
import umcg.genetica.math.matrix2.DoubleMatrixDataset;

/**
 *
 * @author patri
 */
public class SingleCellNetwork {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		DoubleMatrixDataset<String, String> allSingleCellData = DoubleMatrixDataset.loadDoubleTextData("/groups/umcg-wijmenga/scr02/eqtlGen/expression_all_cells.tsv", '\t');

		final CSVParser parser = new CSVParserBuilder().withSeparator('\t').withIgnoreQuotations(true).build();
		final CSVReader reader = new CSVReaderBuilder(new InputStreamReader(new FileInputStream("/groups/umcg-wijmenga/scr02/eqtlGen/tsne_rot_ident_all_cells.tsv"))).withCSVParser(parser).build();

		HashMap<String, ArrayList<String>> cellTypeToIdMapping = new HashMap<>();

		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			
			String id = nextLine[0];
			String cellType = nextLine[5];

			ArrayList<String> cellTypeIds = cellTypeToIdMapping.get(cellType);
			if (cellTypeIds == null) {
				System.out.println(cellType);
				cellTypeIds = new ArrayList<>();
				cellTypeToIdMapping.put(cellType, cellTypeIds);
			}
			cellTypeIds.add(id);

		}

		System.out.println("Done loading data");

		String cellTypeToUse = "CD4+ T";

		for (String cellId : cellTypeToIdMapping.get(cellTypeToUse)) {
			if(!allSingleCellData.containsCol(cellId)){
				System.out.println("could not find cell: " + cellId);
			}
		}

		DoubleMatrixDataset<String, String> currentCellTypeExpression = allSingleCellData.viewColSelection(cellTypeToIdMapping.get(cellTypeToUse));

		int numberOfCells = currentCellTypeExpression.columns();
		ArrayList<String> allGenes = currentCellTypeExpression.getRowObjects();
		
		System.out.println("Number of cells: " + numberOfCells);
		
		ArrayList<String> expressedGenes = new ArrayList<>();
		
		for(int row = 0 ; row < currentCellTypeExpression.rows() ; ++row){
			
			int nonZeroCount = currentCellTypeExpression.getRow(row).cardinality();
			
			if( (nonZeroCount / (double) numberOfCells) >= 0.05 ){
				expressedGenes.add(allGenes.get(row));
			}
			
		}
		
		final DoubleMatrixDataset<String, String> currentCellTypeExpression2 = currentCellTypeExpression.viewRowSelection(expressedGenes);
		
		System.out.println("Number of expressed genes: " + currentCellTypeExpression2.rows());
		
		NaturalRanking ranking = new NaturalRanking(NaNStrategy.FAILED, TiesStrategy.AVERAGE);

		IntStream.range(0, currentCellTypeExpression2.rows()).parallel().forEach(i -> {
			
			DoubleMatrix1D row = currentCellTypeExpression2.viewRow(i);

			double[] ranks = ranking.rank(row.toArray());

			row.assign(ranks);

			if (i % 1000 == 0) {
				System.out.println("Ranking: " + i + " / " + currentCellTypeExpression2.rows());
			}

		});

		System.out.println("Done ranking data");

		DoubleMatrixDataset<String, String> correlationRes = currentCellTypeExpression2.viewDice().calculateCorrelationMatrix();
		
		System.out.println("Done calculating correlations");
		
//		IntStream.range(0, correlationRes.rows()).parallel().forEach(row -> {
//			
//			TDistribution tDistribution = new TDistribution(currentCellTypeExpression2.columns() - 2);
//			
//			for (int col = row; col < correlationRes.columns(); ++col) {
//				double r = correlationRes.getElementQuick(row, col);
//				double t = FastMath.abs(r * FastMath.sqrt((numberOfCells - 2) / (1 - r * r)));
//				double p = 2 * tDistribution.cumulativeProbability(-t);
//				
//				correlationRes.setElementQuick(row, col, p);
//				correlationRes.setElementQuick(col, row, p);
//				
//			}
//			
//			if (row % 1000 == 0) {
//				System.out.println("Cor: " + row + " / " + correlationRes.rows());
//			}
//			
//		});
//
//		System.out.println("Done calculating p-values");
		


		correlationRes.save("/groups/umcg-wijmenga/scr02/eqtlGen/TcellNetwork.txt");

		System.out.println("Data saved");
	}

}
