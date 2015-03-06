package nl.systemsgenetics.annotatecsvensgenes;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import gnu.trove.set.hash.TIntHashSet;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

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
		OptionBuilder.withDescription("");
		OptionBuilder.withLongOpt("input");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("i"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));

		OptionBuilder.withArgName("path");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("First col must be ENSG ID");
		OptionBuilder.withLongOpt("geneRef");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("g"));

		OptionBuilder.withArgName("ints");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("colums with ens ids to annotate (zero based)");
		OptionBuilder.withLongOpt("colums");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("c"));

		OptionBuilder.withArgName("ints");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("reference file colums to use as annotation (zero based)");
		OptionBuilder.withLongOpt("refColums");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("rc"));





	}

	public static void main(String[] args) throws FileNotFoundException, IOException {

		final File inputFile;
		final File outputFile;
		final File referenceFile;
		final TIntHashSet columnsToAnnotate;
		final int[] refColumnsToUse;

		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			inputFile = new File(commandLine.getOptionValue("i"));
			outputFile = new File(commandLine.getOptionValue("o"));
			referenceFile = new File(commandLine.getOptionValue("g"));
			columnsToAnnotate = new TIntHashSet(StringArrayToIntArray(commandLine.getOptionValues("c")));
			refColumnsToUse = StringArrayToIntArray(commandLine.getOptionValues("rc"));

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		final HashMap<String, String[]> ensgAnnotations = new HashMap<String, String[]>();


		CSVReader refReader = new CSVReader(new FileReader(referenceFile), '\t', '\0', '\0');
		final String[] headerRef = refReader.readNext();
		String[] nextLine;
		while ((nextLine = refReader.readNext()) != null) {
			ensgAnnotations.put(nextLine[0], nextLine);
		}

		refReader.close();

		CSVReader inputReader = new CSVReader(new FileReader(inputFile), '\t', '\0', '\0');
		CSVWriter outputWriter = new CSVWriter(new FileWriter(outputFile), '\t', '\0', '\0', "\n");

		String[] inputHeader = inputReader.readNext();
		String[] output = new String[inputHeader.length + (columnsToAnnotate.size() * refColumnsToUse.length)];

		int outputI = 0;
		for (int i = 0; i < inputHeader.length; ++i) {
			output[outputI++] = inputHeader[i];
			if (columnsToAnnotate.contains(i)) {
				for (int j : refColumnsToUse) {
					output[outputI++] = headerRef[j];
				}
			}
		}

		outputWriter.writeNext(output);

		while ((nextLine = inputReader.readNext()) != null) {

			outputI = 0;
			for (int i = 0; i < nextLine.length; ++i) {

				output[outputI++] = nextLine[i];

				if (columnsToAnnotate.contains(i)) {

					String[] annotation = ensgAnnotations.get(nextLine[i]);
					if (annotation == null) {
						for (int j = 0; j < refColumnsToUse.length; ++j) {
							output[outputI++] = "NA";
						}
					} else {
						for (int j : refColumnsToUse) {
							output[outputI++] = annotation[j];
						}
					}


				}


			}

			outputWriter.writeNext(output);

		}

		outputWriter.close();
		inputReader.close();

	}

	private static int[] StringArrayToIntArray(String[] stringArray) {

		int[] intArray = new int[stringArray.length];
		for (int i = 0; i < stringArray.length; ++i) {
			intArray[i] = Integer.parseInt(stringArray[i]);
		}
		return intArray;

	}
}
