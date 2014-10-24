package nl.systemsgenetics.qtlsnpbedannotator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import umcg.genetica.io.bed.BedEntry;
import umcg.genetica.io.bed.BedFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 * Hello world!
 *
 */
public class Main {

	private static final Options OPTIONS;

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
		OptionBuilder.withDescription("Directory with BED files. Track name should be specified in header of file.");
		OptionBuilder.withLongOpt("bedDir");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('b'));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output file");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create('o'));


	}

	public static void main(String[] args) throws IOException, Exception {

		final File qtlFile;
		final File bedDirectory;
		final File outputFile;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			qtlFile = new File(commandLine.getOptionValue('q'));
			bedDirectory = new File(commandLine.getOptionValue('b'));
			outputFile = new File(commandLine.getOptionValue('o'));

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		System.out.println("QTL file: " + qtlFile.getAbsolutePath());
		System.out.println("BED directory: " + bedDirectory.getAbsolutePath());
		System.out.println("Ouput file: " + outputFile.getAbsolutePath());

		if (!qtlFile.canRead()) {
			System.err.println("Failed to read QTL file.");
			System.exit(1);
			return;
		}

		if (!bedDirectory.isDirectory()) {
			if (bedDirectory.exists()) {
				System.err.println("Specified BED directory is a file not a folder");
			} else {
				System.err.println("BED directory does not exist");
			}
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

		File[] bedFiles = bedDirectory.listFiles(new bedFileFilter());
		BedAnnotationSource[] bedAnnotationSources = new BedAnnotationSource[bedFiles.length];

		int i = 0;
		for (File bedFile : bedFiles) {
			BedFile bed = new BedFile(bedFile, true, true);

			String bedTrackName = bed.getTrackInfo().get("name");
			if (bedTrackName == null) {
				System.err.println("Error: bed files from directory must have track name specified in header");
				System.exit(1);
			}

			bedAnnotationSources[i++] = new BedAnnotationSource(bedTrackName, bed.createIntervalTree());

		}

		BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
		writer.append("SNP\tChr\tPos");
		for (BedAnnotationSource bedAnnotationSource : bedAnnotationSources) {
			writer.append('\t');
			writer.append(bedAnnotationSource.getName());
		}
		writer.append('\n');

		eQTLTextFile eQTLsTextFile = new eQTLTextFile(qtlFile.getAbsolutePath(), false);

		for (Iterator<EQTL> eQtlIt = eQTLsTextFile.getEQtlIterator(); eQtlIt.hasNext();) {

			EQTL eQtl = eQtlIt.next();
			String chr = eQtl.getRsChr() + "";
			int pos = eQtl.getRsChrPos();

			writer.append(eQtl.getRsName());
			writer.append('\t');
			writer.append(chr);
			writer.append('\t');
			writer.append(String.valueOf(pos));

			for (BedAnnotationSource bedAnnotationSource : bedAnnotationSources) {
				writer.append('\t');
				List<BedEntry> hits = bedAnnotationSource.searchPosition(chr, pos);
				writer.append(String.valueOf(hits.size()));
			}
			
			writer.append('\n');

		}

		writer.close();
		System.out.println("Done");

	}

	private static class bedFileFilter implements FileFilter {

		public boolean accept(File pathname) {
			return pathname.getName().endsWith(".bed");
		}
	}
}
