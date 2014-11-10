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
import umcg.genetica.io.fasta.ReferenceGenomeFasta;

/**
 *
 * @author Patrick Deelen
 */
public class IlluminaProbeRemapping {

	public static void main(String[] args) throws Exception {

		File samFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\ImmunoProbes.sam");
		String outputPrefix = "D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\ImmunoProbes";
		File mappingReportFile = new File(outputPrefix + "MappingReport.txt");
		File referenceGenomeFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\human_g1k_v37.fasta");
		int maxTotalEditDistance = 5; //Including clipping		


		ReferenceGenomeFasta referenceGenome = new ReferenceGenomeFasta(referenceGenomeFile, ReferenceGenomeFasta.HUMAN_NORMAL_CHR);
		for (String chr : referenceGenome.getChromosomes()) {
			System.out.println(chr);
		}


		SAMFileReader inputSam = new SAMFileReader(samFile);
		TObjectIntHashMap probeMatchedCounter = new TObjectIntHashMap();
		TObjectIntHashMap probeWeakMatchedCounter = new TObjectIntHashMap();
		for (SAMRecord record : inputSam) {
		
			if(getClipping(record.getCigar()).getClippingTotal() + record.getIntegerAttribute("NM") > maxTotalEditDistance){
				probeWeakMatchedCounter.adjustOrPutValue(record.getReadName(), 1, 1);
				continue;
			}
			
			probeMatchedCounter.adjustOrPutValue(record.getReadName(), 1, 1);
		}

		CSVWriter mappingReportWriter = new CSVWriter(new FileWriter(mappingReportFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] mappingReportEntry = new String[16];
		int i = 0;
		mappingReportEntry[i++] = "SNP";
		mappingReportEntry[i++] = "Chr";
		mappingReportEntry[i++] = "SnpPos";
		mappingReportEntry[i++] = "RefAllele";
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
		mappingReportEntry[i++] = "ProbeWeakMatchedCount";
		mappingReportWriter.writeNext(mappingReportEntry);







		inputSam = new SAMFileReader(samFile);
		for (SAMRecord record : inputSam) {

			String snpName = record.getReadName();

			boolean reverseStrand = record.getReadNegativeStrandFlag();

			Cigar cigar = record.getCigar();
			
			Clipping clipping = getClipping(cigar);

			boolean unexpectedCigarOperator = false;
			for (CigarElement cigarElement : cigar.getCigarElements()) {
				if (cigarElement.getOperator() != CigarOperator.M && cigarElement.getOperator() != CigarOperator.S) {
					unexpectedCigarOperator = true;
				}
			}

			boolean clippingAtSnp = (clipping.getClippingLeft() > 0 && reverseStrand) || (clipping.getClippingRight() > 0 && !reverseStrand);

			Integer nm = record.getIntegerAttribute("NM");
			
			if(nm + clipping.getClippingTotal() > maxTotalEditDistance){
				continue;
			}

			int leftProbePos = record.getAlignmentStart() - clipping.getClippingLeft();
			int readLength = record.getReadLength();

			int snpPos = reverseStrand ? leftProbePos - 1 : leftProbePos + readLength;

			String chr = record.getReferenceName();

			String refAllele = referenceGenome.loadedChr(chr) ? String.valueOf(referenceGenome.getNucleotide(chr, snpPos)) : "";

			i = 0;
			mappingReportEntry[i++] = snpName;
			mappingReportEntry[i++] = chr;
			mappingReportEntry[i++] = Integer.toString(snpPos);
			mappingReportEntry[i++] = refAllele;
			mappingReportEntry[i++] = Integer.toString(leftProbePos);
			mappingReportEntry[i++] = reverseStrand ? "-" : "+";
			mappingReportEntry[i++] = Integer.toString(readLength);
			mappingReportEntry[i++] = record.getCigarString();
			mappingReportEntry[i++] = Boolean.toString(unexpectedCigarOperator);
			mappingReportEntry[i++] = Integer.toString(clipping.getClippingLeft());
			mappingReportEntry[i++] = Integer.toString(clipping.getClippingRight());
			mappingReportEntry[i++] = Integer.toString(nm);
			mappingReportEntry[i++] = Integer.toString(nm + clipping.getClippingTotal());
			mappingReportEntry[i++] = Boolean.toString(clippingAtSnp);
			mappingReportEntry[i++] = Integer.toString(probeMatchedCounter.get(snpName));
			mappingReportEntry[i++] = Integer.toString(probeWeakMatchedCounter.get(snpName));

			mappingReportWriter.writeNext(mappingReportEntry);



		}

		mappingReportWriter.close();
		System.out.println("Done");

	}

	private static Clipping getClipping(Cigar cigar) {

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

		return new Clipping(clippingLeft, clippingRight);

	}

	private static class Clipping {

		final int clippingLeft;
		final int clippingRight;

		public Clipping(int clippingLeft, int clippingRight) {
			this.clippingLeft = clippingLeft;
			this.clippingRight = clippingRight;
		}

		public int getClippingLeft() {
			return clippingLeft;
		}

		public int getClippingRight() {
			return clippingRight;
		}
		
		public int getClippingTotal() {
			return clippingLeft + clippingRight;
		}
	}
}
