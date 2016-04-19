/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.collapsegenotypedsamples;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.log4j.Logger;
import org.molgenis.genotype.GenotypeDataException;
import org.molgenis.genotype.GenotypeInfo;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.multipart.IncompatibleMultiPartGenotypeDataException;
import org.molgenis.genotype.tabix.TabixFileNotFoundException;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author patri
 */
public class Main {

	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |      Collapse genotype samples        |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordication Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	private static final Options OPTIONS;
	private static Logger LOGGER;

	static {

		LOGGER = Logger.getLogger(GenotypeInfo.class);

		OPTIONS = new Options();

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The data with the samples to search for in the second dataset");
		OptionBuilder.withLongOpt("genotypes");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("genotypesFormat");
		OPTIONS.addOption(OptionBuilder.create("G"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Path to output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('o'));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The minimum posterior probability to call genotypes. default: " + 0.8);
		OptionBuilder.withLongOpt("posterior");
		OPTIONS.addOption(OptionBuilder.create("p"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The minimum dosage r2 to match sample. default: " + 0.9);
		OptionBuilder.withLongOpt("dosageR");
		OPTIONS.addOption(OptionBuilder.create("r"));
        
	}
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {
		
		System.out.println(HEADER);
		System.out.println();
		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		CommandLineParser parser = new PosixParser();
		final CommandLine commandLine;
		try {
			commandLine = parser.parse(OPTIONS, args, true);
		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: " + ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		final String[] genotypesBasePaths = commandLine.getOptionValues("g");
		final RandomAccessGenotypeDataReaderFormats genotypeDataType;

		try {
			if (commandLine.hasOption("G")) {
				genotypeDataType = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue("I1").toUpperCase());
			} else {
				if (genotypesBasePaths[0].endsWith(".vcf")) {
					System.err.println("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
					System.exit(1);
					return;
				}
				try {
					genotypeDataType = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(genotypesBasePaths[0]);
				} catch (GenotypeDataException e) {
					System.err.println("Unable to determine input 1 type based on specified path. Please specify --input1Type");
					System.exit(1);
					return;
				}
			}
		} catch (IllegalArgumentException e) {
			System.err.println("Error parsing --genotypesFormat \"" + commandLine.getOptionValue("G") + "\" is not a valid input data format");
			System.exit(1);
			return;
		}

		final double minimumPosteriorProbability;
		try {
			minimumPosteriorProbability = commandLine.hasOption("p") ? Double.parseDouble(commandLine.getOptionValue("p")) : 0.8;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing --posterior \"" + commandLine.getOptionValue("ip1") + "\" is not an double");
			System.exit(1);
			return;
		}

		final double minRToMatch;
		try {
			minRToMatch = commandLine.hasOption("r") ? Double.parseDouble(commandLine.getOptionValue("r")) : 0.9;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing -r \"" + commandLine.getOptionValue("r") + "\" is not an double");
			System.exit(1);
			return;
		}
                
		final String outputFilePath = commandLine.getOptionValue('o');

		StringBuilder input1Paths = new StringBuilder();
		for (String path : genotypesBasePaths) {
			input1Paths.append(path);
			input1Paths.append(' ');
		}

		
		System.out.println("Genotype base path: " + input1Paths);
		System.out.println("Genotype data type: " + genotypeDataType.getName());
		System.out.println("Min call probability: " + minimumPosteriorProbability);
		System.out.println("Minimum dosage r to match sampels: " + minRToMatch);
		System.out.println("Output path: " + outputFilePath);

		File outputFile = new File(outputFilePath);
		File outputFolder = outputFile.getParentFile();
		if (outputFolder != null) {
			outputFolder.mkdirs();
		}

		final RandomAccessGenotypeData genotypeData;

		try {
			genotypeData = genotypeDataType.createFilteredGenotypeData(genotypesBasePaths, 100, null, null, null, minimumPosteriorProbability);
		} catch (TabixFileNotFoundException e) {
			LOGGER.fatal("Tabix file not found for input data at: " + e.getPath() + "\n"
					+ "Please see README on how to create a tabix file");
			System.exit(1);
			return;
		} catch (IOException e) {
			LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
			System.exit(1);
			return;
		} catch (IncompatibleMultiPartGenotypeDataException e) {
			LOGGER.fatal("Error combining the impute genotype data files: " + e.getMessage(), e);
			System.exit(1);
			return;
		} catch (GenotypeDataException e) {
			LOGGER.fatal("Error reading input data: " + e.getMessage(), e);
			System.exit(1);
			return;
		}
		
		final String[] samples = genotypeData.getSampleNames();
		
		DenseObjectMatrix2D<SimpleRegression> regressionMatrix = new DenseObjectMatrix2D<>(samples.length, samples.length);

		for (int s1 = 0; s1 < samples.length; ++s1) {
			for (int s2 = 0; s2 < samples.length; ++s2) {
				regressionMatrix.setQuick(s1, s2, new SimpleRegression());
			}
		}
		
		int excludedNonSnpVariants = 0;
		int excludedNonBialelic = 0;
       
//        BufferedWriter logWriter = new BufferedWriter(new FileWriter(outputFile+"_log.txt"));
		for (GeneticVariant variantData : genotypeData) {

			if (variantData.getAlleleCount()>2) {
//                System.out.println(variantData1.getAlleleCount());
				++excludedNonSnpVariants;
				continue;
			}
			
            float[] dosages = null;
            
            dosages = variantData.getSampleDosages();
            

			for (int s1 = 0; s1 < samples.length; ++s1) {

				float variantData1SampleDosage = dosages[s1];

				if (variantData1SampleDosage < 0) {
					continue;
				}

				for (int s2 = 0; s2 < samples.length; ++s2) {

					float variantData2SampleDosage = dosages[s2];

					if (variantData2SampleDosage < 0) {
						continue;
					}
					regressionMatrix.getQuick(s1, s2).addData(variantData1SampleDosage, variantData2SampleDosage);

				}
			}
		}
		
		final DoubleMatrix2D r2 = new DenseDoubleMatrix2D(samples.length, samples.length);
		for (int s1 = 0; s1 < samples.length; ++s1) {
			for (int s2 = 0; s2 < samples.length; ++s2) {
				r2.setQuick(s1, s2, regressionMatrix.getQuick(s1, s2).getRSquare());
			}
		}		
		
		Collapse[] sampleCollapses = new Collapse[samples.length];
		HashSet<Collapse> collapsedSamples = new HashSet<>();
		
		for(int s1 = 0 ; s1 < samples.length ; s1++){
			
			Collapse sampleCollapse = sampleCollapses[s1];
			
			if(sampleCollapse == null){
				sampleCollapse = new Collapse(samples[s1]);
				collapsedSamples.add(sampleCollapse);
			}
			
			for(int s2 = s1 + 1 ; s2 < samples.length ; s2++){
				
				if(r2.getQuick(s1, s2) >= minRToMatch){
					sampleCollapse.addSample(samples[s2]);
					
					if(sampleCollapses[s2] != null){
						throw new RuntimeException("This is a bug");
					}
					
					sampleCollapses[s2] = sampleCollapse;
					
				}
			}
			
			
		}
		
		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(outputFile));
		
		StringBuilder outputLine = new StringBuilder();
		for(Collapse collapse : collapsedSamples){
			
			outputLine.setLength(0);
			
			for(String sample : collapse.getSamples()){
				if(outputLine.length() > 0){
					outputLine.append(';');
				}
				outputLine.append(sample);
			}
			
			outputWriter.write(outputLine.toString());
			
		}
		
		outputWriter.close();
		
		System.out.println("Number of samples after collapsing: " + collapsedSamples.size());
		
		
	}
	
}
