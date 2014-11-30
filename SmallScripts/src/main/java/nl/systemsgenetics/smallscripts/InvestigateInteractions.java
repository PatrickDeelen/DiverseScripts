package nl.systemsgenetics.smallscripts;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import gnu.trove.list.array.TDoubleArrayList;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Patrick Deelen
 */
public class InvestigateInteractions {

	public static void main(String[] args) throws Exception {
		
	
		final int tfColumn = 2;
		final int metaZinteractionColumn = 17;
		
		//final HashMap<String, TDoubleArrayList> overallZ = new HashMap<String, TDoubleArrayList>();
		final HashMap<String, TDoubleArrayList> interactionZs = new HashMap<String, TDoubleArrayList>();
	
		
		CSVReader interactionReader = new CSVReader(new FileReader("D:\\UMCG\\Genetica\\Projects\\rp3\\oppositeEffects\\interactionWithTFs\\LL+LLS+RS+CODAM_meta-analysis.txt"), '\t', '\0', 1);
		
		String[] lineElements;
		while ((lineElements = interactionReader.readNext()) != null) {
			
			String tf = lineElements[tfColumn];
			double interactionMetaZ = Double.parseDouble(lineElements[metaZinteractionColumn]);
			
			TDoubleArrayList tfInteractionZs = interactionZs.get(tf);
			if(tfInteractionZs == null){
				tfInteractionZs = new TDoubleArrayList();
				interactionZs.put(tf, tfInteractionZs);
			}
			
			tfInteractionZs.add(interactionMetaZ);
			
		}
		
		int numberOfQtlsPerTf = -1;		
		for(TDoubleArrayList tfInteractionZs : interactionZs.values()){
			if(numberOfQtlsPerTf == -1){
				numberOfQtlsPerTf = tfInteractionZs.size();
			} else if(numberOfQtlsPerTf != numberOfQtlsPerTf) {
				throw new Exception();
			}
		}
		
		
		CSVWriter writer = new CSVWriter(new FileWriter("D:\\UMCG\\Genetica\\Projects\\rp3\\oppositeEffects\\interactionWithTFs\\interactionMetaZ.txt"), '\t', '\0');
		
		final String[] entries = new String[numberOfQtlsPerTf + 1];
		
		
		for(Map.Entry<String, TDoubleArrayList> tfInteractionZsEntry : interactionZs.entrySet()){
			
			TDoubleArrayList tfInteractionZs = tfInteractionZsEntry.getValue();
			entries[0] = tfInteractionZsEntry.getKey();
			for(int i = 1 ; i <= numberOfQtlsPerTf ; ++i){
				entries[i] = String.valueOf(tfInteractionZs.getQuick(i - 1));
			}
			writer.writeNext(entries);
			
		}
		
		writer.close();
		
	}
	
	
}
