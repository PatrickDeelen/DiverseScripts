package nl.systemsgenetics.smallscripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Patrick Deelen
 */
public class ExtractSnpChrPosFromUcscDbsnpTable {
	
	private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	private static final Pattern CHR_PATTERN = Pattern.compile("^chr(.*)$", Pattern.CASE_INSENSITIVE);
	
	public static void main(String[] args) throws Exception {
		
		if(args.length != 3){
			System.out.println("file with rs on rows");
			System.out.println("uscs dbSNP table");
			System.out.println("output file");
		}
		
		final File inputRsFile = new File(args[0]);
		final File dbSnpFile = new File(args[1]);
		final File outputFile = new File(args[2]);
		
		System.out.println("Input: " + inputRsFile.getAbsolutePath());
		System.out.println("dbSNP table: " + dbSnpFile.getAbsolutePath());
		System.out.println("Output: " + outputFile.getAbsolutePath());
		
		final LinkedHashSet<String> rsIdsToQuery = new LinkedHashSet<String>();
		
		
		BufferedReader inputRsreader = new BufferedReader(new InputStreamReader(new FileInputStream(inputRsFile), "UTF-8"));
		
		String line;
		while ((line = inputRsreader.readLine()) != null) {
			rsIdsToQuery.add(line);
		}
		
		System.out.println("Found " + rsIdsToQuery.size() + " RS IDs to query");
		
		BufferedWriter outputWriter = new BufferedWriter(new FileWriter(outputFile));
		
		BufferedReader dbSnpReader = new BufferedReader(new InputStreamReader(new FileInputStream(dbSnpFile), "UTF-8"));
		while ((line = dbSnpReader.readLine()) != null) {
			if(line.charAt(0) == '#'){
				continue;
			}
			String[] elements = TAB_PATTERN.split(line);
			if(elements.length < 3){
				System.err.println("Error skipping line: " + line);
				continue;
			}
			String rs = elements[2];
			
			if(rsIdsToQuery.contains(rs)){
				rsIdsToQuery.remove(rs);
				
				outputWriter.append(rs);
				outputWriter.append('\t');
				outputWriter.append(removeChr(elements[0]));
				outputWriter.append('\t');
				outputWriter.append(elements[1]);
				outputWriter.append('\n');
			}
			
		}
		outputWriter.close();
		System.out.println("Rs without mapping: " + rsIdsToQuery.size());
	}
	
	private static String removeChr(String chromosome) {

		Matcher chrMatcher = CHR_PATTERN.matcher(chromosome);
		if (chrMatcher.find()) {
			return chrMatcher.group(1);
		} else {
			return chromosome;
		}

	}

}


