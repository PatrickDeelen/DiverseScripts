package nl.systemsgenetics.qtlannotator;

import au.com.bytecode.opencsv.CSVReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.String;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.collections.ChrPosMap;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 * Hello world!
 *
 */
public class QtlAnnotator {

	private static final Options OPTIONS;
	private static final int DEFAULT_WINDOW = 250000;
	private static final double DEFAULT_R2 = 0.9;

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
		OptionBuilder.withDescription("Output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
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
		OptionBuilder.withDescription("Search window for LD SNPs");
		OptionBuilder.withLongOpt("window");
		OPTIONS.addOption(OptionBuilder.create("w"));

		OptionBuilder.withArgName("double");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum r2 to include annotation");
		OPTIONS.addOption(OptionBuilder.create("r2"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Annotation file. chr<tab>pos<tab>annotation1<tab>annotation2 etc, use header to name the annotations");
		OptionBuilder.withLongOpt("annotation");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('a'));

	}

	@SuppressWarnings("RedundantStringConstructorCall")
	public static void main(String[] args) throws IOException, Exception {

		final File qtlFile;
		final File outputFile;
		final File annotationFile;
		final String genotypeDataType;
		final String[] genotypeDataPaths;
		final double minR2;
		final int window;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			qtlFile = new File(commandLine.getOptionValue('q'));
			outputFile = new File(commandLine.getOptionValue('o'));
			annotationFile = new File(commandLine.getOptionValue('a'));
			genotypeDataType = commandLine.getOptionValue('G');
			genotypeDataPaths = commandLine.getOptionValues('g');

			if (commandLine.hasOption("r2")) {
				minR2 = Integer.parseInt(commandLine.getOptionValue("r2"));
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
		
		System.out.println("QTL file: " + qtlFile);
		System.out.println("Annotion file: " + annotationFile);
		System.out.println("Output file: " + outputFile);
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

		
		
		
		ChrPosMap<String[]> annotationMap = new ChrPosMap<String[]>();
		CSVReader reader = new CSVReader(new FileReader(annotationFile));
		String[] nextLine = reader.readNext();
		
		if(nextLine.length < 3){
			System.out.println("No annotations found in header");
		}
		final String[] annotatationNames = new String[nextLine.length - 2];
		for(int i = 2 ; i < nextLine.length ; ++i){
			annotatationNames[i-2] = nextLine[i];
		}
		
		while ((nextLine = reader.readNext()) != null) {
			if(nextLine.length != annotatationNames.length + 2){
				throw new Exception("Error in annotation different number of columns then headers on line: " + Arrays.toString(nextLine));
			}
			annotationMap.put(nextLine[1], Integer.valueOf(nextLine[2]), Arrays.copyOfRange(nextLine, 2, nextLine.length));
		}
		reader.close();
		
		RandomAccessGenotypeData genotypeData = RandomAccessGenotypeDataReaderFormats.valueOf(genotypeDataType.toUpperCase()).createGenotypeData(genotypeDataPaths, 10000);

		QTLTextFile eQTLsTextFile = new QTLTextFile(qtlFile.getAbsolutePath(), false);
		qtls:
		for (Iterator<EQTL> eQtlIt = eQTLsTextFile.getEQtlIterator(); eQtlIt.hasNext();) {

			EQTL eQtl = eQtlIt.next();
			String chr = eQtl.getRsChr() + "";
			int pos = eQtl.getRsChrPos();
			
			GeneticVariant qtlVariant = genotypeData.getSnpVariantByPos(chr, pos);
			
			if(qtlVariant == null){
				System.err.println("ERROR: qtl variant " + eQtl.getRsName() + " not found in reference. Skipping variant.");
				continue qtls;
			}
			
			otherVariants:
			for(GeneticVariant variant : genotypeData.getVariantsByRange(chr, pos - window, pos + window)){
				
				if(variant == qtlVariant){
					continue otherVariants;
				}
				
				
				
			}

		}

	}
}
