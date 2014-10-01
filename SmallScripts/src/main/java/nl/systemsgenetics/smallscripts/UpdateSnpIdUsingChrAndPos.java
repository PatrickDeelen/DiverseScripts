package nl.systemsgenetics.smallscripts;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import umcg.genetica.collections.ChrPosMap;

/**
 *
 * @author Patrick Deelen
 */
public class UpdateSnpIdUsingChrAndPos {

	private static final Pattern CHR_PATTERN = Pattern.compile("^chr(.*)$", Pattern.CASE_INSENSITIVE);
	
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		String dbSnpFile = args[0];

		ChrPosMap<String> snpIds = new ChrPosMap<String>();
		BufferedReader dbSnpReader = new BufferedReader(new InputStreamReader(new FileInputStream(dbSnpFile), "UTF-8"));

		String line;

		int i = 0;
		while ((line = dbSnpReader.readLine()) != null) {

			
			String[] elements = StringUtils.split(line, '\t');

			if (!elements[11].equals("single")) {
				continue;
			}
			snpIds.put(removeChr(elements[1]), Integer.parseInt(elements[2])+1, new String(elements[4]));

			++i;
			if(i %1000000 == 0){
				System.out.println("Parsed " + i);
			}
			

		}


		String inputFile = args[1];
		String outputFile = args[2];
		char inputFileSep = '\t';
		int inputFileChrCol = 5;
		int inputFilePosCol = 6;
		int inputFileSnpColToOverwrite = 7;

		CSVReader reader = new CSVReader(new FileReader(inputFile), inputFileSep);
		CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);
		
		writer.writeNext(reader.readNext());
		
		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {

			String newId = snpIds.get(nextLine[inputFileChrCol], Integer.parseInt(nextLine[inputFilePosCol]));
			if(newId != null){
				nextLine[inputFileSnpColToOverwrite] = newId;
			}
			writer.writeNext(nextLine);
			
		}
		
		writer.close();

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
