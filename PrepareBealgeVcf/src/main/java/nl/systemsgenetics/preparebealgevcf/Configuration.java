/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.preparebealgevcf;

import java.io.File;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author Patrick Deelen
 */
public class Configuration {

	private final String pathToInputVcf;
	private final String referenceVariantList;
	private final String pathToOuputVcf;
	private final String pathToReferenceExcludeList;
	private final int chunkSize;
	private static final Options OPTIONS;
	private static final HelpFormatter HELP_FORMATTER = new HelpFormatter();

	static {

		OPTIONS = new Options();

		Option option;

		option = OptionBuilder.withArgName("file")
				.hasArgs()
				.withDescription("Study VCF, single chromosome")
				.withLongOpt("studyVcf")
				.isRequired()
				.create("s");
		OPTIONS.addOption(option);
		
		option = OptionBuilder.withArgName("file")
				.hasArgs()
				.withDescription("Ref variants for chromosome matching study VCF. chr<tab>pos<tab>ref<tab>alt")
				.withLongOpt("refVariants")
				.isRequired()
				.create("r");
		OPTIONS.addOption(option);
		
		option = OptionBuilder.withArgName("path")
				.hasArgs()
				.withDescription("Filtered study VCF")
				.withLongOpt("outputVcf")
				.isRequired()
				.create("o");
		OPTIONS.addOption(option);
		
		option = OptionBuilder.withArgName("int")
				.hasArgs()
				.withDescription("Chunk size when imputing")
				.withLongOpt("chunkSize")
				.isRequired()
				.create("c");
		OPTIONS.addOption(option);
		
		option = OptionBuilder.withArgName("path")
				.hasArgs()
				.withDescription("List of variants to exclude using bealge option excludemarkers")
				.withLongOpt("excludedMarkers")
				.isRequired()
				.create("e");
		OPTIONS.addOption(option);

	}

	public Configuration(String[] args) throws ParseException {
		
		CommandLineParser parser = new PosixParser();
		CommandLine commandLine = parser.parse(OPTIONS, args, true);
	
		pathToInputVcf = commandLine.getOptionValue('s');
		referenceVariantList = commandLine.getOptionValue('r');
		pathToOuputVcf = commandLine.getOptionValue('o');
		pathToReferenceExcludeList = commandLine.getOptionValue('e');
		chunkSize = Integer.parseInt(commandLine.getOptionValue('c'));
		
	}
	
	public static void printHelp(){
		HELP_FORMATTER.printHelp(" ", OPTIONS);
	}

	public String getPathToInputVcf() {
		return pathToInputVcf;
	}

	public String getReferenceVariantList() {
		return referenceVariantList;
	}

	public String getPathToOuputVcf() {
		return pathToOuputVcf;
	}

	public String getPathToReferenceExcludeList() {
		return pathToReferenceExcludeList;
	}

	public int getChunkSize() {
		return chunkSize;
	}
	
	
	
}
