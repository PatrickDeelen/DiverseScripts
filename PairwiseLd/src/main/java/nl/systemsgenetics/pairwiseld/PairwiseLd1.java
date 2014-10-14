/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.pairwiseld;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.MultiPartGenotypeData;
import org.molgenis.genotype.util.ChrPos;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.containers.Pair;

/**
 *
 * @author Patrick Deelen
 */
public class PairwiseLd1 {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, LdCalculatorException {

		if (args.length != 3) {
			System.out.println("Args:");
			System.out.println(" - SNPs to query file");
			System.out.println(" - Folder with vcf.gz files and vcg.gz.tbx files");
			System.out.println(" - Output file path");
			System.exit(0);
		}

		String querySnpsFilePath = args[0];
		String vcfFolder = args[1];
		String outputFilePath = args[2];

		System.out.println("SNPs to query file: " + querySnpsFilePath);
		System.out.println("Folder with vcf.gz files and vcg.gz.tbx files: " + vcfFolder);
		System.out.println("Output file path: " + outputFilePath);

		ArrayList<Pair<ChrPos, ChrPos>> querySnps = new ArrayList<Pair<ChrPos, ChrPos>>();
		
		CSVReader querySnpReader = new CSVReader(new FileReader(querySnpsFilePath), '\t');
		String[] nextLine;
		
		ArrayList<ChrPos> snps = new ArrayList<ChrPos>();
		
		while ((nextLine = querySnpReader.readNext()) != null) {
			snps.add(new ChrPos(nextLine[1], Integer.parseInt(nextLine[2])));
		}
		
		for(ChrPos snp1 : snps){
			for(ChrPos snp2 : snps){
				if(snp1 == snp2){
					//Calculate each pair once by stopping on same snp.
					//IE take lower triangle of matrix
					break;
				}
				querySnps.add(new Pair<ChrPos, ChrPos>(snp1, snp2));
			}
		}
		
		

		querySnpReader.close();

		System.out.println("Number of query SNP pairs loaded: " + querySnps.size());

		RandomAccessGenotypeData genotypeData = MultiPartGenotypeData.createFromVcfFolder(new File(vcfFolder), 100, 0);


		CSVWriter outputWriter = new CSVWriter(new FileWriter(outputFilePath), '\t', CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.NO_ESCAPE_CHARACTER);

		String[] outputHeaders = {"snp1", "snp2", "r2", "dPrime", "snp1Chr", "snp1Pos", "snp2Chr", "snp2Pos", "distance", "snp1Maf", "snp2Maf"};

		outputWriter.writeNext(outputHeaders);

		for (Pair<ChrPos, ChrPos> snpPair : querySnps) {

			GeneticVariant snp1 = genotypeData.getSnpVariantByPos(snpPair.getLeft().getChr(), snpPair.getLeft().getPos());
			GeneticVariant snp2 = genotypeData.getSnpVariantByPos(snpPair.getRight().getChr(), snpPair.getRight().getPos());

			if (snp1 == null) {
				System.err.println("Error snp1 not found: " + snpPair.getLeft().getChr() + ":" + snpPair.getLeft().getPos() );
				continue;
			}

			if (snp2 == null) {
				System.err.println("Error snp2 not found: " + snpPair.getRight().getChr() + ":" + snpPair.getRight().getPos() );
				continue;
			}
			
			String snp1Id = snp1.getPrimaryVariantId();
			String snp2Id = snp2.getPrimaryVariantId();
	
			String snp1Chr = snp1.getSequenceName();
			int snp1Pos = snp1.getStartPos();

			String snp2Chr = snp2.getSequenceName();
			int snp2Pos = snp2.getStartPos();

			Ld ld = snp1.calculateLd(snp2);

			String[] results = new String[11];
			results[0] = snp1Id;
			results[1] = snp2Id;
			results[2] = String.valueOf(ld.getR2());
			results[3] = String.valueOf(ld.getDPrime());
			results[4] = String.valueOf(snp1Chr);
			results[5] = String.valueOf(snp1Pos);
			results[6] = String.valueOf(snp2Chr);
			results[7] = String.valueOf(snp2Pos);
			results[8] = String.valueOf(Math.abs(snp1Pos - snp2Pos));
			results[9] = String.valueOf(snp1.getMinorAlleleFrequency());
			results[10] = String.valueOf(snp2.getMinorAlleleFrequency());
			outputWriter.writeNext(results);


		}

		outputWriter.close();

	}
}
