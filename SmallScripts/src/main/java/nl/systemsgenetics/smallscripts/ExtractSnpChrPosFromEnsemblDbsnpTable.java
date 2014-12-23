package nl.systemsgenetics.smallscripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author Patrick Deelen
 */
public class ExtractSnpChrPosFromEnsemblDbsnpTable {

	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private static final Pattern INPUT_LINE_PATTERN = Pattern.compile("([^\\t]+)(\\t?.*)");

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException, Exception {

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

		LinkedHashMap<String, ArrayList<String>> querySnps = new LinkedHashMap<String, ArrayList<String>>();

		BufferedReader inputSnpReader = new BufferedReader(new InputStreamReader(new FileInputStream(inputRsFile), "UTF-8"));

		String line;

		while ((line = inputSnpReader.readLine()) != null) {

			Matcher m = INPUT_LINE_PATTERN.matcher(line);
			if (!m.matches()) {
				throw new Exception("Line: " + line);
			}
			String snp = m.group(1);
			ArrayList<String> snpStrings = querySnps.get(snp);
			if (snpStrings == null) {
				snpStrings = new ArrayList<String>(1);
				querySnps.put(snp, snpStrings);
			}
			snpStrings.add(m.group(2));

		}

		inputSnpReader.close();

		System.out.println("Unique input SNPs: " + querySnps.size());

		BufferedReader dbSnpTableReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(dbSnpFile)), "UTF-8"));
		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(outputFile));

		while ((line = dbSnpTableReader.readLine()) != null) {

			String[] elements = TAB_PATTERN.split(line);
			if (elements.length < 3) {
				System.err.println("Error skipping line: " + line);
				continue;
			}
			ArrayList<String> snpLines = querySnps.remove(elements[0]);

			if (snpLines != null) {
				for (String snpLine : snpLines) {
					outputWriter.append(elements[1]);
					outputWriter.append('\t');
					outputWriter.append(elements[2]);
					outputWriter.append('\t');
					outputWriter.append(elements[0]);
					outputWriter.append(snpLine);
					outputWriter.append('\n');
				}
			}

			if (elements.length == 4) {
				snpLines = querySnps.remove(elements[3]);
				if (snpLines != null) {
					for (String snpLine : snpLines) {
						outputWriter.append(elements[1]);
						outputWriter.append('\t');
						outputWriter.append(elements[2]);
						outputWriter.append('\t');
						outputWriter.append(elements[0]);
						outputWriter.append(snpLine);
						outputWriter.append('\n');
					}
				}

			}



		}

		System.out.println("SNPs not found: " + querySnps.size());

		outputWriter.close();
		dbSnpTableReader.close();

	}
}
