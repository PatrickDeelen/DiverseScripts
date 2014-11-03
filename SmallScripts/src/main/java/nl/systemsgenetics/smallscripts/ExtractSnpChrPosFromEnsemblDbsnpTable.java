package nl.systemsgenetics.smallscripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedHashSet;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author Patrick Deelen
 */
public class ExtractSnpChrPosFromEnsemblDbsnpTable {

	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		if (args.length != 3) {
			System.out.println("file with rs on rows");
			System.out.println("dbSNP table");
			System.out.println("output file");
		}

		final File inputRsFile = new File(args[0]);
		final File dbSnpFile = new File(args[1]);
		final File outputFile = new File(args[2]);

		System.out.println("Input: " + inputRsFile.getAbsolutePath());
		System.out.println("dbSNP table: " + dbSnpFile.getAbsolutePath());
		System.out.println("Output: " + outputFile.getAbsolutePath());
		
		LinkedHashSet<String> querySnps = new LinkedHashSet<String>();
		
		BufferedReader inputSnpReader = new BufferedReader(new InputStreamReader(new FileInputStream(inputRsFile), "UTF-8"));
		
		String line;
		
		while ((line = inputSnpReader.readLine()) != null) {
			querySnps.add(line);
		}
		
		inputSnpReader.close();
		
		System.out.println("Unique input SNPs: " + querySnps.size());
		
		BufferedReader dbSnpTableReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(dbSnpFile)), "UTF-8"));
		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(outputFile));
		
		while ((line = dbSnpTableReader.readLine()) != null) {
			
			String[] elements = TAB_PATTERN.split(line);
			if(elements.length < 3){
				System.err.println("Error skipping line: " + line);
				continue;
			}
			
			if(querySnps.remove(elements[0]) || ( elements.length == 4 && querySnps.remove(elements[3] ) )){
				outputWriter.append(elements[0]);
				outputWriter.append('\t');
				outputWriter.append(elements[1]);
				outputWriter.append('\t');
				outputWriter.append(elements[2]);
				outputWriter.append('\n');
			}
			
		}
		
		System.out.println("SNPs not found: " + querySnps.size());
		
		outputWriter.close();
		dbSnpTableReader.close();

	}
}
