package nl.systemsgenetics.smallscripts;

import au.com.bytecode.opencsv.CSVWriter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import umcg.genetica.io.gtf.GffElement;
import umcg.genetica.io.gtf.GtfReader;

/**
 *
 * @author Patrick Deelen
 */
public class GtfToAnnotation {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		GtfReader gtfReader = new GtfReader(new File(args[0]));

		CSVWriter annotationWriter = new CSVWriter(new FileWriter(args[1]), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] outRow = new String[48];
		int i = 0;
		outRow[i++] = "platform";
		outRow[i++] = "exon";
		outRow[i++] = "gene";
		outRow[i++] = "chr";
		outRow[i++] = "begin";
		outRow[i++] = "end";

		annotationWriter.writeNext(outRow);

		outRow[0] = "meta-exons_v71";

		for (GffElement gtfElement : gtfReader) {

			i = 1;
			outRow[i++] = gtfElement.getAttributeValue("meta-exon_id");
			outRow[i++] = gtfElement.getAttributeValue("gene_id");
			outRow[i++] = gtfElement.getSeqname();
			outRow[i++] = String.valueOf(gtfElement.getStart());
			outRow[i++] = String.valueOf(gtfElement.getEnd());
			annotationWriter.writeNext(outRow);

		}
		
		annotationWriter.close();
	}
}
