package nl.systemsgenetics.smallscripts;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.variant.GeneticVariant;
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
		File illuminaMappingFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\Immuno_BeadChip_11419691_B.csv");
		int maxTotalEditDistance = 5; //Including clipping	
		int probeBasesToCheckForMismatch = 10; // Bases nearest to SNP to check for mismatch
		double mafCutoffOverlappingVariants = 0.01;


		String vcf1000gPath = "D:\\UMCG\\Genetica\\Data\\ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites";
		String vcfFolderGoNLPath = "D:\\UMCG\\Genetica\\Data\\GoNL\\release5.2\\01_public_SNVs_and_InDels";

		RandomAccessGenotypeData g1000 = RandomAccessGenotypeDataReaderFormats.VCF.createGenotypeData(vcf1000gPath);
		RandomAccessGenotypeData gonl = RandomAccessGenotypeDataReaderFormats.VCF_FOLDER.createGenotypeData(vcfFolderGoNLPath);

		HashMap<String, IlluminaProbeInfo> illuminaProbeInfo = new HashMap<String, IlluminaProbeInfo>();

		CSVReader reader = new CSVReader(new FileReader(illuminaMappingFile), ',', '\0');
		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			if (nextLine.length == 19 && !nextLine[0].equals("IlmnID")) {
				illuminaProbeInfo.put(new String(nextLine[1]), new IlluminaProbeInfo(nextLine[3], nextLine[16]));
			}
		}










		ReferenceGenomeFasta referenceGenome = new ReferenceGenomeFasta(referenceGenomeFile, ReferenceGenomeFasta.HUMAN_NORMAL_CHR);
		for (String chr : referenceGenome.getChromosomes()) {
			System.out.println(chr);
		}


		SAMFileReader inputSam = new SAMFileReader(samFile);
		TObjectIntHashMap probeMatchedCounter = new TObjectIntHashMap();
		TObjectIntHashMap probeWeakMatchedCounter = new TObjectIntHashMap();
		for (SAMRecord record : inputSam) {

			if (getClipping(record.getCigar()).getClippingTotal() + record.getIntegerAttribute("NM") > maxTotalEditDistance) {
				probeWeakMatchedCounter.adjustOrPutValue(record.getReadName(), 1, 1);
				continue;
			}

			probeMatchedCounter.adjustOrPutValue(record.getReadName(), 1, 1);
		}

		CSVWriter mappingReportWriter = new CSVWriter(new FileWriter(mappingReportFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] mappingReportEntry = new String[38];
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
		mappingReportEntry[i++] = "MismatchProbeEnd";
		mappingReportEntry[i++] = "ClippingAtSnp";
		mappingReportEntry[i++] = "ProbeMatchedCount";
		mappingReportEntry[i++] = "ProbeWeakMatchedCount";
		mappingReportEntry[i++] = "Illumina_a1";
		mappingReportEntry[i++] = "Illumina_a2";
		mappingReportEntry[i++] = "Illumina_isSnp";
		mappingReportEntry[i++] = "Illumina_source_a1";
		mappingReportEntry[i++] = "Illumina_source_a2";
		mappingReportEntry[i++] = "1000G_Ref";
		mappingReportEntry[i++] = "1000G_Alt";
		mappingReportEntry[i++] = "1000G_AF";
		mappingReportEntry[i++] = "1000G_AlleleCount";
		mappingReportEntry[i++] = "1000G_Snp";
		mappingReportEntry[i++] = "1000G_MultiVarAtPos";
		mappingReportEntry[i++] = "1000G_InDelOverLap";
		mappingReportEntry[i++] = "1000G_VarInProbeEnd";
		mappingReportEntry[i++] = "GoNL_Ref";
		mappingReportEntry[i++] = "GoNL_Alt";
		mappingReportEntry[i++] = "GoNL_AF";
		mappingReportEntry[i++] = "GoNL_AlleleCount";
		mappingReportEntry[i++] = "GoNL_Snp";
		mappingReportEntry[i++] = "GoNL_MultiVarAtPos";
		mappingReportEntry[i++] = "GoNL_InDelOverLap";
		mappingReportEntry[i++] = "GoNL_VarInProbeEnd";
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

			if (nm + clipping.getClippingTotal() > maxTotalEditDistance) {
				continue;
			}

			int leftProbePos = record.getAlignmentStart() - clipping.getClippingLeft();
			int probeLength = record.getReadLength();

			int snpPos = reverseStrand ? leftProbePos - 1 : leftProbePos + probeLength;

			String chr = record.getReferenceName();

			//Check for mismatch in n bases next to SNP

			int mismatchesProbeEnd;
			if (referenceGenome.loadedChr(chr)) {
				mismatchesProbeEnd = 0;

				String probeSeq = record.getReadString();

				int checkStartProbe = reverseStrand ? 0 : probeLength - probeBasesToCheckForMismatch;
				int checkStartGenome = reverseStrand ? leftProbePos : leftProbePos + probeLength - probeBasesToCheckForMismatch;

//				if (snpName.equals("imm_6_20856706") || snpName.equals("1kg_6_111885392") || snpName.equals("1kg_1_100995654") || snpName.equals("1kg_1_101019989")) {
//					
//					System.out.println("-----");
//					System.out.println("snp: " + snpName);
//					System.out.println("leftProbePos: " + leftProbePos);
//					System.out.println("reverse: " + reverseStrand);
//					System.out.println("cigar: " + cigar);
//					System.out.println("edit: " + nm);
//					System.out.println("snp clipping: " + clippingAtSnp);
//					System.out.println("probe: " + probeSeq);
//					
//					StringBuilder refSeq = new StringBuilder();
//					for (int x = leftProbePos; x < leftProbePos + probeLength; ++x) {
//						refSeq.append(referenceGenome.getNucleotide(chr, x));
//					}
//					System.out.println("ref:   " + refSeq);
//					
//				}

				for (int j = 0; j < probeBasesToCheckForMismatch; ++j) {

					char probeN = probeSeq.charAt(checkStartProbe + j);
					char refN = referenceGenome.getNucleotide(chr, checkStartGenome + j);

//					if (snpName.equals("imm_6_20856706") || snpName.equals("1kg_6_111885392") || snpName.equals("1kg_1_100995654") || snpName.equals("1kg_1_101019989")) {
//						System.out.println(probeN + "-" + refN);
//					}

					if (probeN != refN) {
						++mismatchesProbeEnd;
					}

				}

			} else {
				mismatchesProbeEnd = -1;
			}


			String refAllele = referenceGenome.loadedChr(chr) ? String.valueOf(referenceGenome.getNucleotide(chr, snpPos)) : "";

			GeneticVariant probe1000gVariant = null;
			boolean probe1000gOverlap = false;
			boolean multipleVarAtPos = false;

			for (GeneticVariant g1000Var : g1000.getVariantsByPos(chr, snpPos)) {

				if (g1000Var.getStartPos() == snpPos) {
					if (probe1000gVariant == null) {
						probe1000gVariant = g1000Var;
					} else {
						multipleVarAtPos = true;
					}
				} else {
					for (float af : (ArrayList<Float>) g1000Var.getAnnotationValues().get("EUR_AF")) {
						if (af >= mafCutoffOverlappingVariants) {
							probe1000gOverlap = true;
						}
					}
				}

			}

			String refAllele1000g = "";
			String altAllele1000g = "";
			float af1000g = Float.NaN;
			int alleleCount = 0;
			boolean snp = false;

			if (probe1000gVariant != null && !multipleVarAtPos) {

				if (probe1000gVariant.isBiallelic()) {
					refAllele1000g = probe1000gVariant.getVariantAlleles().get(0).toString();
					altAllele1000g = probe1000gVariant.getVariantAlleles().get(1).toString();
					af1000g = ((ArrayList<Float>) probe1000gVariant.getAnnotationValues().get("EUR_AF")).get(0);
				}
				alleleCount = probe1000gVariant.getVariantAlleles().getAlleleCount();
				snp = probe1000gVariant.isSnp();

			}



			GeneticVariant probeGonlVariant = null;
			boolean probeGonlOverlap = false;
			boolean multipleVarAtPosGonl = false;

			for (GeneticVariant gonlVar : gonl.getVariantsByPos(chr, snpPos)) {

				if (gonlVar.getStartPos() == snpPos) {
					if (probeGonlVariant == null) {
						probeGonlVariant = gonlVar;
					} else {
						multipleVarAtPosGonl = true;
					}
				} else {
					for (float af : (ArrayList<Float>) gonlVar.getAnnotationValues().get("AF")) {
						if (af >= mafCutoffOverlappingVariants) {
							probeGonlOverlap = true;
						}
					}
				}

			}

			String refAlleleGonl = "";
			String altAlleleGonl = "";
			float afGonl = Float.NaN;
			int alleleCountGonl = 0;
			boolean snpGonl = false;

			if (probeGonlVariant != null && !multipleVarAtPosGonl) {

				if (probeGonlVariant.isBiallelic()) {
					refAlleleGonl = probeGonlVariant.getVariantAlleles().get(0).toString();
					altAlleleGonl = probeGonlVariant.getVariantAlleles().get(1).toString();
					afGonl = ((ArrayList<Float>) probeGonlVariant.getAnnotationValues().get("AF")).get(0);
				}
				alleleCountGonl = probeGonlVariant.getVariantAlleles().getAlleleCount();
				snpGonl = probeGonlVariant.isSnp();

			}



			//Testing voor variants in probe end
			int beginQuery;
			int endQuery;

			if (reverseStrand) {

				beginQuery = snpPos + 1;
				endQuery = snpPos + 1 + probeBasesToCheckForMismatch;

			} else {

				beginQuery = snpPos - probeBasesToCheckForMismatch - 1;
				endQuery = snpPos;

			}

			boolean g1000SnpInProbeEnd = false;
			for (GeneticVariant g1000Var : g1000.getVariantsByRange(chr, beginQuery, endQuery)) {
				for (float af : (ArrayList<Float>) g1000Var.getAnnotationValues().get("EUR_AF")) {
					if (af >= mafCutoffOverlappingVariants) {
						g1000SnpInProbeEnd = true;
					}
				}
			}

			boolean gonlSnpInProbeEnd = false;
			for (GeneticVariant gonlVar : gonl.getVariantsByRange(chr, beginQuery, endQuery)) {
				for (float af : (ArrayList<Float>) gonlVar.getAnnotationValues().get("AF")) {
					if (af >= mafCutoffOverlappingVariants) {
						gonlSnpInProbeEnd = true;
					}
				}
			}

			IlluminaProbeInfo illuminaInfo = illuminaProbeInfo.get(snpName);


			i = 0;
			mappingReportEntry[i++] = snpName;
			mappingReportEntry[i++] = chr;
			mappingReportEntry[i++] = Integer.toString(snpPos);
			mappingReportEntry[i++] = refAllele;
			mappingReportEntry[i++] = Integer.toString(leftProbePos);
			mappingReportEntry[i++] = reverseStrand ? "-" : "+";
			mappingReportEntry[i++] = Integer.toString(probeLength);
			mappingReportEntry[i++] = record.getCigarString();
			mappingReportEntry[i++] = Boolean.toString(unexpectedCigarOperator);
			mappingReportEntry[i++] = Integer.toString(clipping.getClippingLeft());
			mappingReportEntry[i++] = Integer.toString(clipping.getClippingRight());
			mappingReportEntry[i++] = Integer.toString(nm);
			mappingReportEntry[i++] = Integer.toString(nm + clipping.getClippingTotal());
			mappingReportEntry[i++] = mismatchesProbeEnd == -1 ? "NA" : Integer.toString(mismatchesProbeEnd);
			mappingReportEntry[i++] = Boolean.toString(clippingAtSnp);
			mappingReportEntry[i++] = Integer.toString(probeMatchedCounter.get(snpName));
			mappingReportEntry[i++] = Integer.toString(probeWeakMatchedCounter.get(snpName));
			mappingReportEntry[i++] = String.valueOf(illuminaInfo.getA1());
			mappingReportEntry[i++] = String.valueOf(illuminaInfo.getA2());
			mappingReportEntry[i++] = String.valueOf(illuminaInfo.isSnp());
			mappingReportEntry[i++] = String.valueOf(illuminaInfo.getA1_b());
			mappingReportEntry[i++] = String.valueOf(illuminaInfo.getA2_b());
			mappingReportEntry[i++] = refAllele1000g;
			mappingReportEntry[i++] = altAllele1000g;
			mappingReportEntry[i++] = String.valueOf(af1000g);
			mappingReportEntry[i++] = String.valueOf(alleleCount);
			mappingReportEntry[i++] = String.valueOf(snp);
			mappingReportEntry[i++] = String.valueOf(multipleVarAtPos);
			mappingReportEntry[i++] = String.valueOf(probe1000gOverlap);
			mappingReportEntry[i++] = String.valueOf(g1000SnpInProbeEnd);
			mappingReportEntry[i++] = refAlleleGonl;
			mappingReportEntry[i++] = altAlleleGonl;
			mappingReportEntry[i++] = String.valueOf(afGonl);
			mappingReportEntry[i++] = String.valueOf(alleleCountGonl);
			mappingReportEntry[i++] = String.valueOf(snpGonl);
			mappingReportEntry[i++] = String.valueOf(multipleVarAtPosGonl);
			mappingReportEntry[i++] = String.valueOf(probeGonlOverlap);
			mappingReportEntry[i++] = String.valueOf(gonlSnpInProbeEnd);

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

	private static char complement(char n) throws Exception {

		switch (n) {

			case 'A':
				return 'T';
			case 'T':
				return 'A';
			case 'C':
				return 'G';
			case 'G':
				return 'C';

			case 'a':
				return 't';
			case 't':
				return 'a';
			case 'c':
				return 'g';
			case 'g':
				return 'c';

			default:
				throw new Exception("Unexpected allele: " + n);

		}

	}

	private static class IlluminaProbeInfo {

		private final char a1;
		private final char a2;
		private final String a1_b;
		private final String a2_b;
		private static final Pattern ALLELE_PATTERN = Pattern.compile("\\[(.*)/(.*)\\]");

		public IlluminaProbeInfo(String SNP, String SourceSeq) {

			a1 = SNP.charAt(1);
			a2 = SNP.charAt(3);

			Matcher alleleMatcher = ALLELE_PATTERN.matcher(SourceSeq);
			alleleMatcher.find();
			a1_b = alleleMatcher.group(1);
			a2_b = alleleMatcher.group(2);

		}

		public char getA1() {
			return a1;
		}

		public char getA2() {
			return a2;
		}

		public String getA1_b() {
			return a1_b;
		}

		public String getA2_b() {
			return a2_b;
		}

		public boolean isSnp() {
			return isSnpAllele(a1) && isSnpAllele(a2);
		}

		private boolean isSnpAllele(char a) {
			return a == 'A' || a == 'T' || a == 'G' || a == 'C';
		}
	}
}
