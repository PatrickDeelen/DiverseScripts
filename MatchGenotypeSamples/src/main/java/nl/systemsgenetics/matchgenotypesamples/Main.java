package nl.systemsgenetics.matchgenotypesamples;


import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.procedure.TIntProcedure;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
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
import umcg.genetica.containers.DenseGenericObjectMatrix2D;

/**
 * Hello world!
 *
 */
public class Main {

	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |        Match Genotype Samples         |\n"
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
		OptionBuilder.withLongOpt("input1");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("i1"));

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
		OptionBuilder.withLongOpt("input1Type");
		OPTIONS.addOption(OptionBuilder.create("I1"));

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The data with the samples to search for in the second dataset");
		OptionBuilder.withLongOpt("input2");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("i2"));

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
		OptionBuilder.withLongOpt("input2Type");
		OPTIONS.addOption(OptionBuilder.create("I2"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Path to output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('o'));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The minimum posterior probability to call genotypes in the input 1 data. default: " + 0.8);
		OptionBuilder.withLongOpt("input1Prob");
		OPTIONS.addOption(OptionBuilder.create("ip1"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The minimum posterior probability to call genotypes in the input 2 data. default: " + 0.8);
		OptionBuilder.withLongOpt("input2Prob");
		OPTIONS.addOption(OptionBuilder.create("ip2"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The minimum dosage r to match sample. default: " + 0.95);
		OptionBuilder.withLongOpt("dosageR");
		OPTIONS.addOption(OptionBuilder.create("r"));


	}

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

		final String[] input1BasePaths = commandLine.getOptionValues("i1");
		final RandomAccessGenotypeDataReaderFormats input1Type;

		try {
			if (commandLine.hasOption("I1")) {
				input1Type = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue("I1").toUpperCase());
			} else {
				if (input1BasePaths[0].endsWith(".vcf")) {
					System.err.println("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
					System.exit(1);
					return;
				}
				try {
					input1Type = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(input1BasePaths[0]);
				} catch (GenotypeDataException e) {
					System.err.println("Unable to determine input 1 type based on specified path. Please specify --input1Type");
					System.exit(1);
					return;
				}
			}
		} catch (IllegalArgumentException e) {
			System.err.println("Error parsing --input1Type \"" + commandLine.getOptionValue("I1") + "\" is not a valid input data format");
			System.exit(1);
			return;
		}

		final String[] input2BasePaths = commandLine.getOptionValues("i2");
		final RandomAccessGenotypeDataReaderFormats input2Type;

		try {
			if (commandLine.hasOption("I2")) {
				input2Type = RandomAccessGenotypeDataReaderFormats.valueOf(commandLine.getOptionValue("I2").toUpperCase());
			} else {
				if (input2BasePaths[0].endsWith(".vcf")) {
					System.err.println("Only vcf.gz is supported. Please see manual on how to do create a vcf.gz file.");
					System.exit(1);
					return;
				}
				try {
					input2Type = RandomAccessGenotypeDataReaderFormats.matchFormatToPath(input2BasePaths[0]);
				} catch (GenotypeDataException e) {
					System.err.println("Unable to determine input2 type based on specified path. Please specify --input2Type");
					System.exit(1);
					return;
				}
			}
		} catch (IllegalArgumentException e) {
			System.err.println("Error parsing --input2Type \"" + commandLine.getOptionValue("I2") + "\" is not a valid input data format");
			System.exit(1);
			return;
		}

		final double input1MinimumPosteriorProbability;
		try {
			input1MinimumPosteriorProbability = commandLine.hasOption("ip1") ? Double.parseDouble(commandLine.getOptionValue("ip1")) : 0.8;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing --input1Prob \"" + commandLine.getOptionValue("ip1") + "\" is not an double");
			System.exit(1);
			return;
		}

		final double input2MinimumPosteriorProbability;
		try {
			input2MinimumPosteriorProbability = commandLine.hasOption("ip2") ? Double.parseDouble(commandLine.getOptionValue("ip2")) : 0.8;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing --input2Prob \"" + commandLine.getOptionValue("ip2") + "\" is not an double");
			System.exit(1);
			return;
		}

		final double minRToMatch;
		try {
			minRToMatch = commandLine.hasOption("r") ? Double.parseDouble(commandLine.getOptionValue("r")) : 0.95;
		} catch (NumberFormatException e) {
			System.err.println("Error parsing -r \"" + commandLine.getOptionValue("r") + "\" is not an double");
			System.exit(1);
			return;
		}

		final String outputFilePath = commandLine.getOptionValue('o');

		StringBuilder input1Paths = new StringBuilder();
		for (String path : input1BasePaths) {
			input1Paths.append(path);
			input1Paths.append(' ');
		}

		StringBuilder input2Paths = new StringBuilder();
		for (String path : input2BasePaths) {
			input2Paths.append(path);
			input2Paths.append(' ');
		}

		System.out.println("Input 1 base path: " + input1Paths);
		System.out.println("Input 1 data type: " + input1Type.getName());
		System.out.println("Input 1 min call probability: " + input1MinimumPosteriorProbability);
		System.out.println("Input 2 base path: " + input2Paths);
		System.out.println("Input 2 data type: " + input2Type.getName());
		System.out.println("Input 2 min call probability: " + input2MinimumPosteriorProbability);
		System.out.println("Minimum dosage r to match sampels: " + minRToMatch);
		System.out.println("Output path: " + outputFilePath);

		File outputFile = new File(outputFilePath);
		File outputFolder = outputFile.getParentFile();
		if (outputFolder != null) {
			outputFolder.mkdirs();
		}


		final RandomAccessGenotypeData data1;
		final RandomAccessGenotypeData data2;

		try {
                    //System.out.println("Creating Genotype Data input 1: " + " time: " + System.currentTimeMillis());
			data1 = input1Type.createFilteredGenotypeData(input1BasePaths, 100, null, null, null, input1MinimumPosteriorProbability);
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

		try {
                      //  System.out.println("Creating Genotype Data input 2" + " time: " + System.currentTimeMillis());
			data2 = input2Type.createFilteredGenotypeData(input2BasePaths, 100, null, null, null, input2MinimumPosteriorProbability);
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

		String[] data1Samles = data1.getSampleNames();
		String[] data2Samles = data2.getSampleNames();

		DenseRegressionMatrix2D regressionMatrix = new DenseRegressionMatrix2D(data1Samles.length, data2Samles.length);
		//System.out.println("Loading in Samples now: " + " Time: " + System.currentTimeMillis());
                for (int s1 = 0; s1 < data1Samles.length; ++s1) {
			for (int s2 = 0; s2 < data2Samles.length; ++s2) {
				regressionMatrix.setQuick(s1, s2, new SimpleRegression());
			}
		}

		int excludedNonSnpVariants = 0;
		int excludedNonBialelic = 0;
		int excludedNotInInput2 = 0;
		int excludedDifferentAllelesInput2 = 0;
		int variantsFoundInBoth = 0;
                
		for (GeneticVariant variantData1 : data1) {                      
                       
			if (!variantData1.isSnp()) {
				++excludedNonSnpVariants;
				continue;
			}

			if (!variantData1.isBiallelic()) {
				++excludedNonBialelic;
				continue;
			}

			String chr = variantData1.getSequenceName();
			int pos = variantData1.getStartPos();

			GeneticVariant variantData2 = data2.getSnpVariantByPos(chr, pos);

			if (variantData2 == null) {
				++excludedNotInInput2;
				continue;
			}

			if (!variantData1.getVariantAlleles().containsAll(variantData2.getVariantAlleles())) {
				++excludedDifferentAllelesInput2;
				continue;
			}

			++variantsFoundInBoth;

			float[] variantData1SampleDosages = variantData1.getSampleDosages();
			float[] variantData2SampleDosages = variantData2.getSampleDosages();


                        //System.out.println("Entering sample loop NEW JAR: " + " Time: " + System.currentTimeMillis());
                        
			for (int s1 = 0; s1 < data1Samles.length; ++s1) {
                            
//                            if (j % 500 == 0) {
//                                    System.out.println("Sample: " + j + " Time: " + System.currentTimeMillis());
//                            }
//                            ++j; 
                            
				float variantData1SampleDosage = variantData1SampleDosages[s1];

				if (variantData1SampleDosage < 0) {
					continue;
				}

				for (int s2 = 0; s2 < data2Samles.length; ++s2) {
					float variantData2SampleDosage = variantData2SampleDosages[s2];

					if (variantData2SampleDosage < 0) {
						continue;
					}
                                        //The slow step...?
					regressionMatrix.getQuick(s1, s2).addData(variantData1SampleDosage, variantData2SampleDosage);

				}
			}

		}

		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(outputFile));

		int lowVarWarning = 0;
		int sampleMatched = 0;
		int sampleMultipleMatched = 0;
		TObjectIntMap<String> matchedSamples = new TObjectIntHashMap<String>(data1Samles.length, 0.5f, 0);

		for (int s1 = 0; s1 < data1Samles.length; ++s1) {

			boolean matchFound = false;
			boolean multipleMatched = false;

			outputWriter.append(data1Samles[s1]);
			outputWriter.append('\t');
			
			StringBuilder samplesMatchedString = new StringBuilder();
			StringBuilder samplesMatchedRString = new StringBuilder();

			for (int s2 = 0; s2 < data2Samles.length; ++s2) {

				SimpleRegression regression = regressionMatrix.getQuick(s1, s2);
				if (regression.getN() < 5000) {
					++lowVarWarning;
				}

				if (regression.getR() >= minRToMatch) {
					if (matchFound) {
						samplesMatchedString.append(',');
						samplesMatchedRString.append(',');
						multipleMatched = true;
					}
					matchFound = true;
					String data2Sample = data2Samles[s2];
					samplesMatchedString.append(data2Sample);
					samplesMatchedRString.append(Double.toString(regression.getR()));
					matchedSamples.adjustOrPutValue(data2Sample, 1, 1);
				}

			}
			outputWriter.append(samplesMatchedString);
			outputWriter.append('\t');
			outputWriter.append(samplesMatchedRString);
			outputWriter.append('\n');

			if (matchFound) {
				++sampleMatched;
			}

			if (multipleMatched) {
				++sampleMultipleMatched;
			}
		}

		outputWriter.close();
		
		CountBiggerThanOne matchedToMultipleCounter = new CountBiggerThanOne();
		matchedSamples.forEachValue(matchedToMultipleCounter);
		
			
		System.out.println("Matching procedure completed");
		System.out.println();
		System.out.println("Variant excluded from input 1:");
		System.out.println(" - Not a SNP: " + excludedNonSnpVariants);
		System.out.println(" - Not bi-allelic: " + excludedNonBialelic);
		System.out.println(" - Not found in input 2: " + excludedNotInInput2);
		System.out.println(" - Different alleles in input 2: " + excludedDifferentAllelesInput2 + " (might be strand issues)");
		System.out.println("Variants matched between input 1 and 2: " + variantsFoundInBoth);
		System.out.println();
		System.out.println("Comparisons using < 5000 variants: " + lowVarWarning);
		System.out.println("Samples with match in input 2: " + sampleMatched + " out of " + data1Samles.length + " of which " + sampleMultipleMatched + " match to multiple samples in input 2");
		System.out.println("Samples in input 2 matched multiple times: " + matchedToMultipleCounter.getCount());

	}

	private static class CountBiggerThanOne implements TIntProcedure {

		private int count = 0;
		
		public boolean execute(int value) {
			if(value > 1){
				++count;	
			} 
			return true;
			
		}

		public int getCount() {
			return count;
		}
		
		
	}
}
