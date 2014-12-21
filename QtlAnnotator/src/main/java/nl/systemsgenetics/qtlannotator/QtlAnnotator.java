package nl.systemsgenetics.qtlannotator;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.util.Ld;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.collections.ChrPosMap;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

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




		ChrPosMap<HashSet<String>[]> annotationMap = new ChrPosMap<HashSet<String>[]>();
		CSVReader reader = new CSVReader(new FileReader(annotationFile));
		String[] nextLine = reader.readNext();

		if (nextLine.length < 3) {
			System.out.println("No annotations found in header");
		}
		final String[] annotatationNames = new String[nextLine.length - 2];
		for (int i = 2; i < nextLine.length; ++i) {
			annotatationNames[i - 2] = nextLine[i];
		}

		while ((nextLine = reader.readNext()) != null) {
			if (nextLine.length != annotatationNames.length + 2) {
				throw new Exception("Error in annotation different number of columns then headers on line: " + Arrays.toString(nextLine));
			}

			int pos = Integer.valueOf(nextLine[2]);
			HashSet<String>[] currentAnnotationPos = annotationMap.get(nextLine[1], pos);
			if (currentAnnotationPos == null) {
				currentAnnotationPos = new HashSet[annotatationNames.length];
				for (int i = 0; i < currentAnnotationPos.length; ++i) {
					currentAnnotationPos[i] = new HashSet<String>(1);
				}
				annotationMap.put(nextLine[1], pos, currentAnnotationPos);
			}
			for (int i = 0; i < currentAnnotationPos.length; ++i) {
				currentAnnotationPos[i].add(nextLine[i + 2]);
			}


		}
		reader.close();

		RandomAccessGenotypeData genotypeData = RandomAccessGenotypeDataReaderFormats.valueOf(genotypeDataType.toUpperCase()).createGenotypeData(genotypeDataPaths, 10000);

		QTLTextFile eQTLsTextFile = new QTLTextFile(qtlFile.getAbsolutePath(), false);



		CSVWriter mappingReportWriter = new CSVWriter(new FileWriter(outputFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] qtlAnnotatedOutput = new String[5 + annotatationNames.length];
		int c = 0;
		qtlAnnotatedOutput[c++] = "p-value";
		qtlAnnotatedOutput[c++] = "variant";
		qtlAnnotatedOutput[c++] = "chr";
		qtlAnnotatedOutput[c++] = "pos";

		for (String annotationName : annotatationNames) {
			qtlAnnotatedOutput[c++] = annotationName;
		}
		mappingReportWriter.writeNext(qtlAnnotatedOutput);

		final StringBuilder[] qtlAnnotations = new StringBuilder[annotatationNames.length];

		qtls:
		for (Iterator<EQTL> eQtlIt = eQTLsTextFile.getEQtlIterator(); eQtlIt.hasNext();) {

			EQTL eQtl = eQtlIt.next();
			String chr = eQtl.getRsChr() + "";
			int pos = eQtl.getRsChrPos();

			for (int i = 0; i < annotatationNames.length; ++i) {
				HashSet<String>[] annotation = annotationMap.get(chr, pos);
				if (annotation == null) {
					qtlAnnotations[i] = new StringBuilder();
				} else {
					qtlAnnotations[i] = annotationsToString(annotation[i]);
				}

			}

			GeneticVariant qtlVariant = genotypeData.getSnpVariantByPos(chr, pos);

			if (qtlVariant == null) {
				System.err.println("ERROR: QTL variant " + eQtl.getRsName() + " not found in reference. Unable to search for variants in LD for this QTL.");
			} else {

				otherVariants:
				for (GeneticVariant variant : genotypeData.getVariantsByRange(chr, pos - window, pos + window)) {

					if (variant == qtlVariant) {
						continue otherVariants;
					}

					Ld ld = variant.calculateLd(qtlVariant);

					if (ld.getR2() < minR2) {
						continue otherVariants;
					}

					HashSet<String>[] annotation = annotationMap.get(variant.getSequenceName(), variant.getStartPos());

					if (annotation != null) {
						for (int i = 0; i < annotatationNames.length; ++i) {
							if (qtlAnnotations[i].length() != 0) {
								qtlAnnotations[i].append(';');
							}
							qtlAnnotations[i].append(annotationsToString(annotation[i]));
						}
					}

				}

			}

			c = 0;
			qtlAnnotatedOutput[c++] = String.valueOf(eQtl.getPvalue());
			qtlAnnotatedOutput[c++] = String.valueOf(eQtl.getRsName());
			qtlAnnotatedOutput[c++] = chr;
			qtlAnnotatedOutput[c++] = String.valueOf(pos);

			for (StringBuilder annotation : qtlAnnotations) {
				qtlAnnotatedOutput[c++] = annotation.toString();
			}

			mappingReportWriter.writeNext(qtlAnnotatedOutput);

		}

		mappingReportWriter.close();

		System.out.println("Annotation done");

	}
	
	private static StringBuilder annotationsToString(HashSet<String> annotations){
		
		StringBuilder result = new StringBuilder();
		
		boolean first = true;
		for(String annotation : annotations){
			if(first){
				first = false;
			} else {
				result.append(',');
			}
			result.append(annotation);
		}
		
		return result;
		
	}
}
