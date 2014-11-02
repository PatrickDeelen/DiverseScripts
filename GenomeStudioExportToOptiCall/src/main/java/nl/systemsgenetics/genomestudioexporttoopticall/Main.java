package nl.systemsgenetics.genomestudioexporttoopticall;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;

/**
 * Hello world!
 *
 */
public class Main {

	private static final Pattern SAMPLE_COUNT_PATTERN = Pattern.compile("Total Samples\\s+(\\d+)");
	private static final Pattern SNP_COUNT_PATTERN = Pattern.compile("Total SNPs\\s+(\\d+)");

	@SuppressWarnings("RedundantStringConstructorCall")
	public static void main(String[] args) throws Exception {
		
		int rsColumn = 0;
		int sampleColumn = 1;
		int chrColumn = 18;
		int chrPosColumn = 19;
		int allelesColumn = 22;
		int aSignalColumn =  29;
		int bSignalColumn =  30;
				
		System.out.println("Input file: " + args[0]);
		System.out.println("Output folder: " + args[1]);
		
		File outputFolder = new File(args[1]);
		if(!outputFolder.isDirectory()){
			if(!outputFolder.mkdirs()){
				throw new Exception("Cannot create output dir");
			}
		}
		
		if(!outputFolder.canWrite()){
			throw new Exception("Cannot write to output folder");
		}
		
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(args[0]), "UTF-8"));

		String line;

		line = reader.readLine();

		if (!line.equals("[Header]")) {
			throw new Exception("Expected: [Header]");
		}

		boolean inHeader = true;
		int sampleCount = 0;
		int snpCount = 0;
		int stringBuilderSize = 0;
		HashMap<String, TreeMap<Integer, StringBuilder>> outputData = new HashMap<String, TreeMap<Integer, StringBuilder>>(25);
		LinkedHashSet<String> samples = null;
		
		int excludedNaSnps = 0;
		int excludedInvalidSnps = 0;
		
		String currentSample = "";

		while ((line = reader.readLine()) != null) {

			if (inHeader) {

				Matcher sampleCountMatcher = SAMPLE_COUNT_PATTERN.matcher(line);

				if (sampleCountMatcher.matches()) {
					sampleCount = Integer.parseInt(sampleCountMatcher.group(1));
					continue;
				}
				
				Matcher snpCountMatcher = SNP_COUNT_PATTERN.matcher(line);
				
				if (snpCountMatcher.matches()) {
					snpCount = Integer.parseInt(snpCountMatcher.group(1));
					continue;
				}

				if (line.equals("[Data]")) {
					
					if(sampleCount == 0){
						throw new Exception("Did not find sample count in header");
					}
					if(snpCount == 0){
						throw new Exception("Did not find SNP count in header");
					}
					
					samples = new LinkedHashSet<String>(sampleCount);
					
					stringBuilderSize = 100 + 12 * sampleCount;
					
					
					//skip header of data
					reader.readLine();
					
					inHeader = false;
					
					
				}


			} else {
				
				String[] elements = StringUtils.splitPreserveAllTokens(line, '\t');
				
				if(elements[allelesColumn].equals("[N/A]")){
					++excludedNaSnps;
					continue;
				}
				
				if(elements[allelesColumn].length() != 5) {
					++excludedInvalidSnps;
					System.err.println("Illegal SNP alleles: " + elements[allelesColumn]);
					continue;
				}
				
				if(!elements[sampleColumn].equals(currentSample)){
					if(samples.contains(elements[sampleColumn])){
						throw new Exception("Sample order problem");
					}
					samples.add(new String(elements[sampleColumn]));
					currentSample = elements[sampleColumn];
				}
				
				TreeMap<Integer, StringBuilder> outputDataChr;
				if(outputData.containsKey(elements[chrColumn])){
					outputDataChr = outputData.get(elements[chrColumn]);
				} else {
					outputDataChr = new TreeMap<Integer, StringBuilder>();
					outputData.put(new String(elements[chrColumn]), outputDataChr);
				}
				
				StringBuilder snpOutput;
				if(!outputDataChr.containsKey(Integer.valueOf(elements[chrPosColumn]))){
					snpOutput = new StringBuilder(stringBuilderSize);
					snpOutput.append(elements[rsColumn]);
					snpOutput.append('\t');
					snpOutput.append(elements[chrPosColumn]);
					snpOutput.append('\t');
					snpOutput.append(elements[allelesColumn].charAt(1));
					snpOutput.append(elements[allelesColumn].charAt(3));
					outputDataChr.put(Integer.valueOf(elements[chrPosColumn]), snpOutput);
				} else {
					snpOutput = outputDataChr.get(Integer.valueOf(elements[chrPosColumn]));
				}
				snpOutput.append('\t');
				snpOutput.append(elements[aSignalColumn]);
				snpOutput.append('\t');
				snpOutput.append(elements[bSignalColumn]);
				
			}

		}
		
		System.out.println("Excluded N/A SNPs: " + excludedNaSnps);
		System.out.println("Excluded invalid alleles SNPs: " + excludedInvalidSnps);
		
		StringBuilder header = new StringBuilder("SNP\tCoor\tAlleles");
		for(String sample : samples){
			header.append('\t');
			header.append(sample);
			header.append('A');
			header.append('\t');
			header.append(sample);
			header.append('B');
		}
		header.append('\n');
		
		String headerString = header.toString();
		
		
		for(Map.Entry<String, TreeMap<Integer, StringBuilder>> chrEntry : outputData.entrySet()){
			
			String chr = chrEntry.getKey();
			TreeMap<Integer, StringBuilder> outputDataChr = chrEntry.getValue();
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(outputFolder, "chr_" + chr + "_intensities.tsv")));
			
			writer.append(headerString);
			
			for(StringBuilder snpOutput : outputDataChr.values()){
				snpOutput.append('\n');
				writer.append(snpOutput);
			}
			
			writer.close();
			
		}
	
	}
}
