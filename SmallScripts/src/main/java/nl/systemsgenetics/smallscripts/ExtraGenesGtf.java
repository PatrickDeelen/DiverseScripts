package nl.systemsgenetics.smallscripts;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.HashSet;
import umcg.genetica.io.gtf.GffElement;
import umcg.genetica.io.gtf.GtfReader;

/**
 *
 * @author Patrick Deelen
 */
public class ExtraGenesGtf {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
        String gtf = "D:\\UMCG\\Genetica\\Projects\\RnaSeqEqtl\\Homo_sapiens.GRCh37.71.cut.sorted.gtf";
		String genesPath = "D:\\UMCG\\Genetica\\Projects\\RnaSeqEqtl\\genes.txt";
		
		GtfReader gftReader = new GtfReader(new File(gtf));
		
		BufferedWriter genesWriter = new BufferedWriter(new FileWriter(genesPath));
		genesWriter.append("Ensembl_ID\tGene_symbol\tBiotype\n");
		
		HashSet<String> foundGenes = new HashSet<String>();
		
		for(GffElement element : gftReader){
			
			String ensId = element.getAttributeValue("gene_id");
			
			if(foundGenes.contains(ensId)){
				continue;
			}
			
			foundGenes.add(ensId);
			
			String symbol = element.getAttributeValue("gene_name");
			String type = element.getAttributeValue("gene_biotype");
			genesWriter.append(ensId);
			genesWriter.append('\t');
			genesWriter.append(symbol);
			genesWriter.append('\t');
			genesWriter.append(type);
			genesWriter.append('\n');
			
		}
		
		genesWriter.close();
    }

}
