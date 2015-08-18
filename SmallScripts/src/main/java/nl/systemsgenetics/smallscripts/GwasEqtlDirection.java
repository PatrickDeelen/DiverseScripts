package nl.systemsgenetics.smallscripts;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.multipart.MultiPartGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author Patrick Deelen
 */
public class GwasEqtlDirection {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException, Exception {


		if (args.length != 3) {
			System.out.println("Args:");
			System.out.println(" - Input");
			System.out.println(" - Folder with vcf.gz files and vcg.gz.tbx files");
			System.out.println(" - Output file path");
			System.exit(0);
		}

		String inputPath = args[0];
		String vcfFolder = args[1];
		String outputFilePath = args[2];


		System.out.println("Input: " + inputPath);
		System.out.println("Folder with vcf.gz files and vcg.gz.tbx files: " + vcfFolder);
		System.out.println("Output file path: " + outputFilePath);

		SimpleRegression regression = new SimpleRegression();

		RandomAccessGenotypeData genotypeData = MultiPartGenotypeData.createFromVcfFolder(new File(vcfFolder), 10000, 0);


		CSVWriter outputWriter = new CSVWriter(new FileWriter(outputFilePath), '\t', CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.NO_ESCAPE_CHARACTER);



		CSVReader inputReader = new CSVReader(new FileReader(inputPath), '\t');
		String[] nextLine = inputReader.readNext();//header

		String[] outputLine = new String[nextLine.length + 6];
		int outputCol = 0;
		for (String element : nextLine) {
			outputLine[outputCol++] = element;
		}
		outputLine[outputCol++] = "eSNP_MA";
		outputLine[outputCol++] = "eSNP_MAF";
		outputLine[outputCol++] = "gSNP_MA";
		outputLine[outputCol++] = "gSNP_MAF";
		outputLine[outputCol++] = "r";
		outputLine[outputCol++] = "eQtlDirectionRiskAllele";
		outputWriter.writeNext(outputLine);

		while ((nextLine = inputReader.readNext()) != null) {

			String eSnpChr = nextLine[7];
			int eSnpPos = Integer.parseInt(nextLine[8]);

			double eSnpZscore = Double.parseDouble(nextLine[3]);
			Allele eSnpAssessedAllele = Allele.create(nextLine[2]);

			String gSnpChr = nextLine[10];
			int gSnpPos = Integer.parseInt(nextLine[11]);

			Allele gSnpRiskAllele = Allele.create(nextLine[13]);

			GeneticVariant eSnp;
			try {
				eSnp = genotypeData.getSnpVariantByPos(eSnpChr, eSnpPos);
			} catch (GenotypeDataException ex) {
				System.err.println("Error reading eSNP " + eSnpChr + ":" + eSnpPos);
				ex.printStackTrace();
				continue;
			}

			if (eSnp == null) {
				System.err.println("Error eSNP not found: " + eSnpChr + ":" + eSnpPos);
				continue;
			}

			if (!(eSnp.getMinorAlleleFrequency() > 0)) {
				System.err.println("Error eSNP with MAF 0: " + eSnpChr + ":" + eSnpPos);
				continue;
			}

			GeneticVariant gSnp;
			try {
				gSnp = genotypeData.getSnpVariantByPos(gSnpChr, gSnpPos);
			} catch (GenotypeDataException ex) {
				System.err.println("Error reading gSNP " + gSnpChr + ":" + gSnpPos);
				ex.printStackTrace();
				continue;
			}

			if (gSnp == null) {
				System.err.println("Error gSNP not found: " + gSnpChr + ":" + gSnpPos);
				continue;
			}

			if (!(gSnp.getMinorAlleleFrequency() > 0)) {
				System.err.println("Error gSNP with MAF 0: " + gSnpChr + ":" + gSnpPos);
				continue;
			}

			if (!eSnp.getVariantAlleles().contains(eSnpAssessedAllele)) {
				throw new Exception("Assessed eQTL snp not found in variant.");
			}

			if(!gSnpRiskAllele.isSnpAllele()){
				System.err.println("Skipping, gwas risk allele not a SNP allle: " + gSnpChr + ":" + gSnpPos + " " + gSnpRiskAllele.getAlleleAsString());
				continue;
			}
			
			if (!gSnp.getVariantAlleles().contains(gSnpRiskAllele)) {
				gSnpRiskAllele = gSnpRiskAllele.getComplement();
			}

			Allele eSnpRef = eSnp.getRefAllele();
			Allele gSnpRef = gSnp.getRefAllele();

			float[] eSnpDosage = eSnp.getSampleDosages();
			float[] gSnpDosage = gSnp.getSampleDosages();

			regression.clear();
			for (int i = 0; i < eSnpDosage.length; i++) {
				regression.addData(eSnpDosage[i], gSnpDosage[i]);
			}
			double r = regression.getR();

			if (eSnpAssessedAllele != eSnpRef) {
				r = r * -1;
			}

			if (gSnpRiskAllele != gSnpRef) {
				r = r * -1;
			}
//
//			final boolean eSnpRefIncreaseExpression;
//
//			//If z = positive assessed allele increases expression
//			if (eSnpAssessedAllele == eSnpRef) {
//				if (eSnpZscore > 0) {
//					//eSnp ref increases expression
//					eSnpRefIncreaseExpression = true;
//				} else {
//					//eSnp ref decreases expression
//					eSnpRefIncreaseExpression = false;
//				}
//			} else {
//				if (eSnpZscore < 0) {
//					//eSnp ref increases expression
//					eSnpRefIncreaseExpression = true;
//				} else {
//					//eSnp ref decreases expression
//					eSnpRefIncreaseExpression = false;
//				}
//			}

			//final boolean gSnpRefIncreaseExpression = r >= 0 ? eSnpRefIncreaseExpression : !eSnpRefIncreaseExpression;

			//final boolean gSnpRiskIncreaseExpression = gSnpRiskAllele == gSnpRef ? gSnpRefIncreaseExpression : !gSnpRefIncreaseExpression;
			
			final boolean gSnpRiskIncreaseExpression = r >= 0 ? eSnpZscore > 0 : eSnpZscore < 0;


			//Test GC

			outputCol = 0;
			for (String element : nextLine) {
				outputLine[outputCol++] = element;
			}
			outputLine[outputCol++] = eSnp.getMinorAllele().getAlleleAsString();
			outputLine[outputCol++] = String.valueOf(eSnp.getMinorAlleleFrequency());
			outputLine[outputCol++] = gSnp.getMinorAllele().getAlleleAsString();
			outputLine[outputCol++] = String.valueOf(eSnp.getMinorAlleleFrequency());
			outputLine[outputCol++] = String.valueOf(r);
			outputLine[outputCol++] = gSnp.isAtOrGcSnp() ? "NA" : (gSnpRiskIncreaseExpression ? "increased" : "decreased");
			outputWriter.writeNext(outputLine);



		}

		outputWriter.close();
		inputReader.close();
		
		System.out.println("Done");


	}
}
