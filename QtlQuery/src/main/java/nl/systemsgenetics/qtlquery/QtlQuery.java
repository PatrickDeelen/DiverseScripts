package nl.systemsgenetics.qtlquery;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.NavigableMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.StringUtils;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.collections.ChrPosMap;
import umcg.genetica.collections.ChrPosTreeMap;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

public class QtlQuery {

	private static final Options OPTIONS;
	private static final int DEFAULT_WINDOW = 250000;
	private static final double DEFAULT_R2 = 0.8;

	private static final String HEADER
			= "  /---------------------------------------\\\n"
			+ "  |             QTL Query                 |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordication Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";

	static {

		OPTIONS = new Options();

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("QTL file");
		OptionBuilder.withLongOpt("qtlres");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('q'));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output file.");
		OptionBuilder.withLongOpt("output");
		OPTIONS.addOption(OptionBuilder.create('o'));

		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("The path to the reference genotypes. These genotypes will be used for LD calculations");
		OptionBuilder.withLongOpt("genotypes");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The input data type.\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("genotypesType");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("G"));

		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Search window for LD SNPs. Default: " + DEFAULT_WINDOW);
		OptionBuilder.withLongOpt("window");
		OPTIONS.addOption(OptionBuilder.create("w"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum r2 to include annotation. Default: " + DEFAULT_R2);
		OPTIONS.addOption(OptionBuilder.create("r2"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Search file. chr<tab>pos<tab>.... use headers");
		OptionBuilder.withLongOpt("search");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('s'));

	}

	@SuppressWarnings("RedundantStringConstructorCall")
	public static void main(String[] args) throws IOException, Exception {

		System.out.println(HEADER);
		System.out.println();
		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}

		final File qtlFile;
		final File outputFile;
		final File queryFile;
		final String genotypeDataType;
		final String[] genotypeDataPaths;
		final double minR2;
		final int window;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			qtlFile = new File(commandLine.getOptionValue('q'));
			queryFile = new File(commandLine.getOptionValue('s'));

			if (commandLine.hasOption("o")) {
				outputFile = new File(commandLine.getOptionValue('o'));
			} else {
				outputFile = new File(queryFile.getAbsolutePath() + "_annotated.txt");
			}

			genotypeDataType = commandLine.getOptionValue('G');
			genotypeDataPaths = commandLine.getOptionValues('g');

			if (commandLine.hasOption("r2")) {
				minR2 = Double.parseDouble(commandLine.getOptionValue("r2"));
			} else {
				minR2 = DEFAULT_R2;
			}

			if (commandLine.hasOption("w")) {
				window = Integer.parseInt(commandLine.getOptionValue("w"));
			} else {
				window = DEFAULT_WINDOW;
			}

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		System.out.println("QTL file: " + qtlFile.getAbsolutePath());
		System.out.println("Annotion file: " + queryFile.getAbsolutePath());
		System.out.println("Output file: " + outputFile.getAbsolutePath());
		System.out.println("Genotype data: " + Arrays.toString(genotypeDataPaths));
		System.out.println("Genotype data type: " + genotypeDataType);
		System.out.println("Windows: " + window);
		System.out.println("r2: " + minR2);

		if (!qtlFile.canRead()) {
			System.err.println("Failed to read QTL file.");
			System.exit(1);
			return;
		}

		File outputParent = outputFile.getParentFile();
		if (outputParent != null) {
			if (!outputParent.exists()) {
				if (!outputParent.mkdirs()) {
					System.err.println("Failed to create directory for output file at: " + outputParent.getAbsolutePath());
				}
			}
		}

		RandomAccessGenotypeData genotypeData = RandomAccessGenotypeDataReaderFormats.valueOfSmart(genotypeDataType.toUpperCase()).createGenotypeData(genotypeDataPaths, 10000);

		QTLTextFile eQTLsTextFile = new QTLTextFile(qtlFile.getAbsolutePath(), false);
		ChrPosTreeMap<ArrayList<EQTL>> eQtls = eQTLsTextFile.readQtlsAsTreeMap();

		CSVReader reader = new CSVReader(new FileReader(queryFile), '\t', '\0');
		String[] queryHeader = reader.readNext();
		String[] queryLine;

		CSVWriter outputWriter = new CSVWriter(new FileWriter(outputFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] outputLine = new String[queryHeader.length + 8];
		int c = 0;

		for (int i = 0; i < queryHeader.length; ++i) {
			outputLine[c++] = queryHeader[i];
		}

		outputLine[c++] = "Query SNP allele ";
		outputLine[c++] = "eQTL_SNP";
		outputLine[c++] = "eQTL_Assed_allele";
		outputLine[c++] = "eQTL_Z-score";
		outputLine[c++] = "eQTL_P-value";
		outputLine[c++] = "eQTL_Gene";
		outputLine[c++] = "R2";
		outputLine[c++] = "D'";

		outputWriter.writeNext(outputLine);

		while ((queryLine = reader.readNext()) != null) {

			String queryChr = queryLine[0];
			int queryPos = Integer.parseInt(queryLine[1]);

			GeneticVariant queryVariant = genotypeData.getSnpVariantByPos(queryChr, queryPos);

			if (queryVariant == null) {
				System.out.println(queryChr + ":" + queryPos + " not found in genotype data");
				continue;
			}

			for (ArrayList<EQTL> posEQtl : eQtls.getChrRange(queryChr, queryPos - window, true, queryPos + window, true).values()) {
				for (EQTL eQtl : posEQtl) {

					GeneticVariant qtlVariant = genotypeData.getSnpVariantByPos(eQtl.getRsChr() + "", eQtl.getRsChrPos());

					if (qtlVariant == null) {
						System.out.println("Missing QTL variant in genotype data");
						continue;
					}

					Ld ld = queryVariant.calculateLd(qtlVariant);

					if (ld.getR2() > minR2) {

						c = 0;

						for (int i = 0; i < queryHeader.length; ++i) {
							outputLine[c++] = queryLine[i];
						}

						String commonHap = null;
						double commonHapFreq = -1;
						for (Map.Entry<String, Double> hapFreq : ld.getHaplotypesFreq().entrySet()) {

							double f = hapFreq.getValue();

							if (f > commonHapFreq) {
								commonHapFreq = f;
								commonHap = hapFreq.getKey();
							}

						}

						String[] commonHapAlleles = StringUtils.split(commonHap, '/');
						
						double zscore = eQtl.getZscore();
						
						if(!commonHapAlleles[1].equals(eQtl.getAlleleAssessed())){
							zscore = zscore * -1;
						}

						
						outputLine[c++]	= commonHapAlleles[0];
						outputLine[c++] = eQtl.getRsName();
						outputLine[c++] = commonHapAlleles[1];
						outputLine[c++] = String.valueOf(zscore);
						outputLine[c++] = String.valueOf(eQtl.getPvalue());
						outputLine[c++] = eQtl.getProbe();
						outputLine[c++] = String.valueOf(ld.getR2());
						outputLine[c++] = String.valueOf(ld.getDPrime());

						outputWriter.writeNext(outputLine);

					}

				}
			}

		}
		reader.close();
		outputWriter.close();

		System.out.println("Query done");

	}
}
