package nl.systemsgenetics.smallscripts;

import au.com.bytecode.opencsv.CSVWriter;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.FileWriter;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 *
 * @author Patrick Deelen
 */
public class IlluminaProbeRemapping {

	public static void main(String[] args) throws Exception {

		File samFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\ImmunoProbes.sam");
		String outputPrefix = "D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\ImmunoProbes";
		File mappingReportFile = new File(outputPrefix + "MappingReport.txt");
		File snpReport = new File(outputPrefix + "SnpReport.txt");

		
		SAMFileReader inputSam = new SAMFileReader(samFile);
		TObjectIntHashMap probeMatchedCounter = new TObjectIntHashMap();
		for (SAMRecord record : inputSam) {
			probeMatchedCounter.adjustOrPutValue(record.getReadName(), 1, 1);
		}
		
		CSVWriter mappingReportWriter = new CSVWriter(new FileWriter(mappingReportFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] mappingReportEntry = new String[14];
		int i = 0;
		mappingReportEntry[i++] = "SNP";
		mappingReportEntry[i++] = "Chr";
		mappingReportEntry[i++] = "SnpPos";
		mappingReportEntry[i++] = "LeftProbePos";
		mappingReportEntry[i++] = "Strand";
		mappingReportEntry[i++] = "ProbeLength";
		mappingReportEntry[i++] = "Cigar";
		mappingReportEntry[i++] = "UnexpectedCigarOpperator";
		mappingReportEntry[i++] = "ClippingLeft";
		mappingReportEntry[i++] = "ClippingRight";
		mappingReportEntry[i++] = "EditDistance";
		mappingReportEntry[i++] = "EditDistanceAndClipping";
		mappingReportEntry[i++] = "ClippingAtSnp";
		mappingReportEntry[i++] = "ProbeMatchedCount";
		mappingReportWriter.writeNext(mappingReportEntry);


		

		
		
		
		inputSam = new SAMFileReader(samFile);
		for (SAMRecord record : inputSam) {

			String snpName = record.getReadName();

			boolean reverseStrand = record.getReadNegativeStrandFlag();

			Cigar cigar = record.getCigar();

			final int clippingLeft;
			final int clippingRight;

			if (cigar.numCigarElements() > 1) {
				CigarElement firstCigarElement = cigar.getCigarElement(0);
				CigarElement lastCigarElement = cigar.getCigarElement(cigar.numCigarElements() - 1);

				if (firstCigarElement.getOperator() == CigarOperator.S) {
					clippingLeft = firstCigarElement.getLength();
				} else {
					clippingLeft = 0;
				}

				if (lastCigarElement.getOperator() == CigarOperator.S) {
					clippingRight = lastCigarElement.getLength();
				} else {
					clippingRight = 0;
				}

			} else {
				clippingLeft = 0;
				clippingRight = 0;
			}

			boolean unexpectedCigarOperator = false;
			for (CigarElement cigarElement : cigar.getCigarElements()) {
				if (cigarElement.getOperator() != CigarOperator.M && cigarElement.getOperator() != CigarOperator.S) {
					unexpectedCigarOperator = true;
				}
			}

			boolean clippingAtSnp = (clippingLeft > 0 && reverseStrand) || (clippingRight > 0 && !reverseStrand);

			Integer nm = record.getIntegerAttribute("NM");

			int leftProbePos = record.getAlignmentStart() - clippingLeft;
			int readLength = record.getReadLength();
			
			int snpPos = reverseStrand ? leftProbePos - 1 : leftProbePos + readLength;

			i = 0;
			mappingReportEntry[i++] = snpName;
			mappingReportEntry[i++] = record.getReferenceName();
			mappingReportEntry[i++] = Integer.toString(snpPos);
			mappingReportEntry[i++] = Integer.toString(leftProbePos);
			mappingReportEntry[i++] = reverseStrand ? "-" : "+";
			mappingReportEntry[i++] = Integer.toString(readLength);
			mappingReportEntry[i++] = record.getCigarString();
			mappingReportEntry[i++] = Boolean.toString(unexpectedCigarOperator);
			mappingReportEntry[i++] = Integer.toString(clippingLeft);
			mappingReportEntry[i++] = Integer.toString(clippingRight);
			mappingReportEntry[i++] = Integer.toString(clippingLeft);
			mappingReportEntry[i++] = Integer.toString(clippingLeft);
			mappingReportEntry[i++] = Boolean.toString(clippingAtSnp);
			mappingReportEntry[i++] = Integer.toString(probeMatchedCounter.get(snpName));
			
			mappingReportWriter.writeNext(mappingReportEntry);



		}

		mappingReportWriter.close();
		System.out.println("Done");

	}
}
