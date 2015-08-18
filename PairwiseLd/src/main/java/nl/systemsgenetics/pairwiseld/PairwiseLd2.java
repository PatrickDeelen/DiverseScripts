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
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.MultiPartGenotypeData;
import org.molgenis.genotype.util.ChrPos;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.util.LdCalculatorException;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class PairwiseLd2 {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, LdCalculatorException {

		if (args.length != 8) {
			System.out.println("Args:");
			System.out.println(" - SNPs to query file tab seperated");
			System.out.println(" - SNPs to query file chr column");
			System.out.println(" - SNPs to query file pos column");
			System.out.println(" - Folder with vcf.gz files and vcg.gz.tbx files");
			System.out.println(" - Query range (bases to left and right)");
			System.out.println(" - Min LD to report (r2)");
			System.out.println(" - Output file path");
			System.out.println(" - SNPs to query file name column");
			System.exit(0);
		}

		String querySnpsFilePath = args[0];
		int chrCol = Integer.parseInt(args[1]);
		int posCol = Integer.parseInt(args[2]);
		String vcfFolder = args[3];
		int queryRange = Integer.parseInt(args[4]);
		double minLd = Double.parseDouble(args[5]);
		String outputFilePath = args[6];
		int nameCol = Integer.parseInt(args[7]);

		System.out.println("SNPs to query file: " + querySnpsFilePath);
		System.out.println("Chr col: " + chrCol);
		System.out.println("Pos col: " + posCol);
		System.out.println("Folder with vcf.gz files and vcg.gz.tbx files: " + vcfFolder);
		System.out.println("Query range (bases to left and right): " + queryRange);
		System.out.println("Min LD to report (r2): " + minLd);	
		System.out.println("Output file path: " + outputFilePath);
		System.out.println("Name col: " + posCol);

		CSVReader querySnpReader = new CSVReader(new FileReader(querySnpsFilePath), '\t');
		String[] nextLine;

		

		//System.out.println("Number of query SNP to query: " + querySnps.size());

		RandomAccessGenotypeData genotypeData = MultiPartGenotypeData.createFromVcfFolder(new File(vcfFolder), 10000, 0);


		CSVWriter outputWriter = new CSVWriter(new FileWriter(outputFilePath), '\t', CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.NO_ESCAPE_CHARACTER);

		String[] outputHeaders = {"snp1", "snp2", "r2", "dPrime", "snp1Chr", "snp1Pos", "snp1Alleles", "snp2Chr", "snp2Pos", "snp2Alleles", "distance", "snp1Maf", "snp1Ma", "snp2Maf", "snp2Ma"};

		outputWriter.writeNext(outputHeaders);

		while ((nextLine = querySnpReader.readNext()) != null) {
			
			String querySnpChr = nextLine[chrCol];
			int querySnpPos = Integer.parseInt(nextLine[posCol]);
			String querySnpName = nextLine[nameCol];
			
			GeneticVariant snp1;
			try {
				snp1 = genotypeData.getSnpVariantByPos(querySnpChr, querySnpPos);
			} catch (GenotypeDataException ex) {
				System.err.println("Error reading SNP " + querySnpChr + ":" + querySnpPos);
				ex.printStackTrace();
				continue;
			}

			if (snp1 == null) {
				System.err.println("Error snp not found: " + querySnpChr + ":" + querySnpPos);
				continue;
			}
			
			if( !(snp1.getMinorAlleleFrequency() > 0) ){
				System.err.println("Error snp with MAF 0: " + querySnpChr + ":" + querySnpPos);
				continue;
			}

			String snp1Id = snp1.getPrimaryVariantId();
			String snp1Chr = snp1.getSequenceName();
			int snp1Pos = snp1.getStartPos();

			int queryStart = snp1Pos - queryRange;
			if(queryStart < 0){
				queryStart = 0;
			}
			
			targetSnps:
			for (GeneticVariant snp2 : genotypeData.getVariantsByRange(snp1Chr, queryStart, snp1Pos + queryRange)) {

				String snp2Id = snp2.getPrimaryVariantId();

				String snp2Chr = snp2.getSequenceName();
				int snp2Pos = snp2.getStartPos();
				
//				if(snp2Pos == snp1Pos){
//					continue targetSnps;
//				}
				
				if( !(snp2.getMinorAlleleFrequency() > 0) ){
					continue targetSnps;
				}

				Ld ld = snp1.calculateLd(snp2);
				
				if(Double.isNaN(ld.getR2())){
					continue targetSnps;
				}
				
				if(ld.getR2() < minLd){
					continue targetSnps;
				}

				String[] results = new String[15];
				int column = 0;
				results[column++] = querySnpName;
				results[column++] = snp2Id;
				results[column++] = String.valueOf(ld.getR2());
				results[column++] = String.valueOf(ld.getDPrime());
				results[column++] = String.valueOf(snp1Chr);
				results[column++] = String.valueOf(snp1Pos);
				results[column++] = snp1.getVariantAlleles().toString();
				results[column++] = String.valueOf(snp2Chr);
				results[column++] = String.valueOf(snp2Pos);
				results[column++] = snp2.getVariantAlleles().toString();
				results[column++] = String.valueOf(Math.abs(snp1Pos - snp2Pos));
				results[column++] = String.valueOf(snp1.getMinorAlleleFrequency());
				results[column++] = snp1.getMinorAllele().getAlleleAsString();
				results[column++] = String.valueOf(snp2.getMinorAlleleFrequency());
				results[column++] = snp2.getMinorAllele().getAlleleAsString();
				outputWriter.writeNext(results);
			}

		}

		outputWriter.close();

		System.out.println("LD done");
		querySnpReader.close();
	}
}

