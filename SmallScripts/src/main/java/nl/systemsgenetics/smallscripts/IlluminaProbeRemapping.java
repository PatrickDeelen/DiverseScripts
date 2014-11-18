package nl.systemsgenetics.smallscripts;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.apache.mahout.math.Arrays;
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

//		File samFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\ImmunoProbes_gapped.sam");
//		File illuminaMappingFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\Immuno_BeadChip_11419691_B.csv");
//		String outputPrefix = "D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\ImmunoProbes_gapped";


//		File samFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\CytoSNP-12v2.1Probes_gapped.sam");
//		File illuminaMappingFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\HumanCytoSNP-12v2-1_L.csv");
//		String outputPrefix = "D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\CytoSNP-12v2.1Probes_gapped";

		File samFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\CytoSNP-12v2Probes_gapped.sam");
		File illuminaMappingFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\HumanCytoSNP-12v2_H.csv");
		String outputPrefix = "D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\CytoSNP-12v2Probes_gapped";

		File referenceGenomeFile = new File("D:\\UMCG\\Genetica\\Projects\\LifeLinesDeep\\genotypingRelease3\\remappingProbes\\human_g1k_v37.fasta");

		int maxTotalEditDistance = 5; //Including clipping	
		int probeBasesToCheckForMismatch = 10; // Bases nearest to SNP to check for mismatch
		double mafCutoffOverlappingVariants = 0.01;

		File mappingReportFile = new File(outputPrefix + "MappingReport.txt");
		File passedProbesFile = new File(outputPrefix + "PassedProbes.txt");
		File failedProbesFile = new File(outputPrefix + "FailedProbes.txt");

		String vcf1000gPath = "D:\\UMCG\\Genetica\\Data\\ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites";
		String vcfFolderGoNLPath = "D:\\UMCG\\Genetica\\Data\\GoNL\\release5.2\\01_public_SNVs_and_InDels";

		RandomAccessGenotypeData g1000 = RandomAccessGenotypeDataReaderFormats.VCF.createGenotypeData(vcf1000gPath);
		RandomAccessGenotypeData gonl = RandomAccessGenotypeDataReaderFormats.VCF_FOLDER.createGenotypeData(vcfFolderGoNLPath);

		HashMap<String, IlluminaProbeInfo> illuminaProbeInfo = new HashMap<String, IlluminaProbeInfo>();
		HashSet<String> probesInManifestNotMapped = new HashSet<String>();

		TObjectIntHashMap<FailReason> failCounter = new TObjectIntHashMap<FailReason>();

		CSVReader reader = new CSVReader(new FileReader(illuminaMappingFile), ',', '\0');
		boolean header = true;
		String[] nextLine;
		while ((nextLine = reader.readNext()) != null) {
			if (header) {
				if (nextLine[0].equals("[Assay]")) {
					header = false;
				}
			} else {
				if (nextLine[0].equals("[Controls]")) {
					break;
				}
				if (nextLine.length >= 17 && !nextLine[0].equals("IlmnID")) {
					try {
						String snp = new String(nextLine[1]);
						probesInManifestNotMapped.add(snp);
						illuminaProbeInfo.put(snp, new IlluminaProbeInfo(nextLine[3], nextLine[16]));
					} catch (Exception e) {
						System.err.println(Arrays.toString(nextLine));
						throw e;
					}
				}
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

		final String[] mappingReportEntry = new String[48];
		int i = 0;
		mappingReportEntry[i++] = "SNP";
		mappingReportEntry[i++] = "Chr";
		mappingReportEntry[i++] = "SnpPos";
		mappingReportEntry[i++] = "RefAllele";
		mappingReportEntry[i++] = "LeftProbePos";
		mappingReportEntry[i++] = "Strand";
		mappingReportEntry[i++] = "ProbeLength";
		mappingReportEntry[i++] = "Cigar";
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
		mappingReportEntry[i++] = "1000G_StartPos";
		mappingReportEntry[i++] = "1000G_AF";
		mappingReportEntry[i++] = "1000G_AlleleCount";
		mappingReportEntry[i++] = "1000G_Snp";
		mappingReportEntry[i++] = "1000G_MultiVarAtPos";
		mappingReportEntry[i++] = "1000G_InDelOverLap";
		mappingReportEntry[i++] = "1000G_VarInProbeEnd";
		mappingReportEntry[i++] = "GoNL_Ref";
		mappingReportEntry[i++] = "GoNL_Alt";
		mappingReportEntry[i++] = "GoNL_StartPos";
		mappingReportEntry[i++] = "GoNL_AF";
		mappingReportEntry[i++] = "GoNL_AlleleCount";
		mappingReportEntry[i++] = "GoNL_Snp";
		mappingReportEntry[i++] = "GoNL_MultiVarAtPos";
		mappingReportEntry[i++] = "GoNL_InDelOverLap";
		mappingReportEntry[i++] = "GoNL_VarInProbeEnd";
		mappingReportEntry[i++] = "VariantRef";
		mappingReportEntry[i++] = "VariantAlt";
		mappingReportEntry[i++] = "A_Allele";
		mappingReportEntry[i++] = "B_Allele";
		mappingReportEntry[i++] = "BiallelicVariant";
		mappingReportEntry[i++] = "SnpVariant";
		mappingReportEntry[i++] = "AlleleMatch";
		mappingReportEntry[i++] = "AlleleMatchedReversed";
		mappingReportEntry[i++] = "ProbePassedQc";
		mappingReportWriter.writeNext(mappingReportEntry);


		CSVWriter passedProbesWriter = new CSVWriter(new FileWriter(passedProbesFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] passedProbeEntry = new String[7];
		i = 0;
		passedProbeEntry[i++] = "Id";
		passedProbeEntry[i++] = "Chr";
		passedProbeEntry[i++] = "Pos";
		passedProbeEntry[i++] = "Ref";
		passedProbeEntry[i++] = "Alt";
		passedProbeEntry[i++] = "A";
		passedProbeEntry[i++] = "B";
		passedProbesWriter.writeNext(passedProbeEntry);

		HashSet<String> failReported = new HashSet<String>();
		CSVWriter failedProbesWriter = new CSVWriter(new FileWriter(failedProbesFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);

		final String[] failedProbeEntry = new String[2];
		i = 0;
		failedProbeEntry[i++] = "Id";
		failedProbeEntry[i++] = "Reason";
		failedProbesWriter.writeNext(failedProbeEntry);




		inputSam = new SAMFileReader(samFile);
		for (SAMRecord record : inputSam) {

			final String snpName = record.getReadName();
			final boolean reverseStrand = record.getReadNegativeStrandFlag();
			final Integer nm = record.getIntegerAttribute("NM");
			final boolean clippingAtSnp;
			final int leftProbePos;
			final int probeLength;
			int varPosition;
			String chr;
			int mismatchesProbeEnd;

			if (!illuminaProbeInfo.containsKey(snpName)) {
				throw new Exception(snpName + " not found in manifest");
			}


			Cigar cigar = record.getCigar();

			Clipping clipping = getClipping(cigar);

//			boolean unexpectedCigarOperator = false;
//			for (CigarElement cigarElement : cigar.getCigarElements()) {
//				if (cigarElement.getOperator() != CigarOperator.M && cigarElement.getOperator() != CigarOperator.S) {
//					unexpectedCigarOperator = true;
//				}
//			}

			clippingAtSnp = (clipping.getClippingLeft() > 0 && reverseStrand) || (clipping.getClippingRight() > 0 && !reverseStrand);


			if (nm + clipping.getClippingTotal() > maxTotalEditDistance) {
				continue;
			}

			probesInManifestNotMapped.remove(snpName);

			leftProbePos = record.getAlignmentStart() - clipping.getClippingLeft();

			probeLength = record.getReadLength();

			varPosition = reverseStrand ? leftProbePos - 1 : leftProbePos + probeLength;

			chr = record.getReferenceName();

			IlluminaProbeInfo illuminaInfo = illuminaProbeInfo.get(snpName);

			if (illuminaInfo.hasNoAlleles()) {

				if (!failReported.contains(snpName)) {
					failReported.add(snpName);

					failCounter.adjustOrPutValue(FailReason.ILLUMINA_EXCLUDED, 1, 1);
					i = 0;
					failedProbeEntry[i++] = snpName;
					failedProbeEntry[i++] = FailReason.ILLUMINA_EXCLUDED.getReason();
					failedProbesWriter.writeNext(failedProbeEntry);

				}

				continue;
			}

			//Check for mismatch in n bases next to SNP


			if (probeMatchedCounter.get(snpName) == 1 && referenceGenome.loadedChr(chr)) {
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


			String refAllele = referenceGenome.loadedChr(chr) ? String.valueOf(referenceGenome.getNucleotide(chr, varPosition)) : "";

			GeneticVariant probe1000gVariant = null;
			boolean probe1000gOverlap = false;
			boolean multipleVarAtPos = false;


			String refAllele1000g = "";
			String altAllele1000g = "";
			int pos1000gVariant = -1;
			float af1000g = Float.NaN;
			int alleleCount1000g = 0;
			boolean snp = false;




			GeneticVariant probeGonlVariant = null;
			boolean probeGonlOverlap = false;
			boolean multipleVarAtPosGonl = false;



			String refAlleleGonl = "";
			String altAlleleGonl = "";
			float afGonl = Float.NaN;
			int posGonlVariant = -1;
			int alleleCountGonl = 0;
			boolean snpGonl = false;






			boolean g1000SnpInProbeEnd = false;
			boolean gonlSnpInProbeEnd = false;

			if (probeMatchedCounter.get(snpName) == 1) {


				for (GeneticVariant g1000Var : g1000.getVariantsByPos(chr, varPosition)) {

					if (g1000Var.getStartPos() == varPosition || (!reverseStrand && !illuminaInfo.isSnp() && g1000Var.getStartPos() == varPosition - 1)) {
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

				if (probe1000gVariant != null && !multipleVarAtPos) {

					if (probe1000gVariant.isBiallelic()) {
						refAllele1000g = probe1000gVariant.getVariantAlleles().get(0).toString();
						altAllele1000g = probe1000gVariant.getVariantAlleles().get(1).toString();
						pos1000gVariant = probe1000gVariant.getStartPos();
						af1000g = ((ArrayList<Float>) probe1000gVariant.getAnnotationValues().get("EUR_AF")).get(0);
					}
					alleleCount1000g = probe1000gVariant.getVariantAlleles().getAlleleCount();
					snp = probe1000gVariant.isSnp();

				}


				for (GeneticVariant gonlVar : gonl.getVariantsByPos(chr, varPosition)) {

					if (gonlVar.getStartPos() == varPosition || (!reverseStrand && !illuminaInfo.isSnp() && gonlVar.getStartPos() == varPosition - 1)) {
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

				if (probeGonlVariant != null && !multipleVarAtPosGonl) {

					if (probeGonlVariant.isBiallelic()) {
						refAlleleGonl = probeGonlVariant.getVariantAlleles().get(0).toString();
						altAlleleGonl = probeGonlVariant.getVariantAlleles().get(1).toString();
						posGonlVariant = probeGonlVariant.getStartPos();
						afGonl = ((ArrayList<Float>) probeGonlVariant.getAnnotationValues().get("AF")).get(0);
					}
					alleleCountGonl = probeGonlVariant.getVariantAlleles().getAlleleCount();
					snpGonl = probeGonlVariant.isSnp();

				}




				//Testing voor variants in probe end
				int beginQuery;
				int endQuery;

				if (reverseStrand) {

					beginQuery = varPosition + 1;
					endQuery = varPosition + probeBasesToCheckForMismatch;

				} else {

					beginQuery = varPosition - probeBasesToCheckForMismatch - 1;
					endQuery = varPosition - 1;

				}

				for (GeneticVariant g1000Var : g1000.getVariantsByRange(chr, beginQuery, endQuery)) {
					if (g1000Var.getStartPos() == pos1000gVariant) {
						continue;//This can happen for indels
					}
					for (float af : (ArrayList<Float>) g1000Var.getAnnotationValues().get("EUR_AF")) {
						if (af >= mafCutoffOverlappingVariants) {
							g1000SnpInProbeEnd = true;
						}
					}
				}



				for (GeneticVariant gonlVar : gonl.getVariantsByRange(chr, beginQuery, endQuery)) {
					if (gonlVar.getStartPos() == posGonlVariant) {
						continue;//This can happen for indels
					}
					for (float af : (ArrayList<Float>) gonlVar.getAnnotationValues().get("AF")) {
						if (af >= mafCutoffOverlappingVariants) {
							gonlSnpInProbeEnd = true;
						}
					}
				}
			}



			final String variantRefAllele;
			final String variantAltAllele;

			final boolean bialleleicVariant;

			if (alleleCount1000g > 2 || alleleCountGonl > 2) {
				//Yes first check this for next loop to work
				variantRefAllele = "";
				variantAltAllele = "";
				bialleleicVariant = false;
			} else {
				if (alleleCount1000g == 2 || alleleCountGonl == 2) {

					if (alleleCount1000g == 2 && alleleCountGonl == 2) {

						if (refAllele1000g.equals(refAlleleGonl) && altAllele1000g.equals(altAlleleGonl)) {
							variantRefAllele = refAllele1000g;
							variantAltAllele = altAllele1000g;
							bialleleicVariant = true;
						} else {
							variantRefAllele = "";
							variantAltAllele = "";
							bialleleicVariant = false;
						}

					} else if (alleleCount1000g == 2) {
						variantRefAllele = refAllele1000g;
						variantAltAllele = altAllele1000g;
						bialleleicVariant = true;
					} else { //alleleCountGonl == 2
						variantRefAllele = refAlleleGonl;
						variantAltAllele = altAlleleGonl;
						bialleleicVariant = true;
					}


				} else {
					variantRefAllele = "";
					variantAltAllele = "";
					bialleleicVariant = false;
				}
			}


			final boolean isSnpVariant = variantAltAllele.length() == 1 && variantRefAllele.length() == 1;

			final boolean alleleMatch;
			final boolean alleleMatchedAfterComplement;

			if (illuminaInfo.isSnp() && isSnpVariant && bialleleicVariant) {

				char aRef = variantRefAllele.charAt(0);
				char aAlt = variantAltAllele.charAt(0);

				char aIllum1 = illuminaInfo.getA1();
				char aIllum2 = illuminaInfo.getA2();

				if ((aRef == aIllum1 && aAlt == aIllum2) || (aAlt == aIllum1 && aRef == aIllum2)) {
					alleleMatch = true;
					alleleMatchedAfterComplement = false;
				} else if ((aRef == complement(aIllum1) && aAlt == complement(aIllum2)) || (aAlt == complement(aIllum1) && aRef == complement(aIllum2))) {
					alleleMatch = true;
					alleleMatchedAfterComplement = true;
				} else {
					alleleMatch = false;
					alleleMatchedAfterComplement = false;
				}

			} else if (!illuminaInfo.isSnp() && !isSnpVariant && bialleleicVariant) {

				if (!illuminaInfo.getA1_b().equals("-")) {
					throw new Exception("");
				}

				if (variantAltAllele.length() == 1) {
					if (illuminaInfo.getA2_b().equals(variantRefAllele.substring(1))) {
						alleleMatch = true;
						alleleMatchedAfterComplement = false;
					} else {
						alleleMatch = false;
						alleleMatchedAfterComplement = false;
					}
				} else {
					if (illuminaInfo.getA2_b().equals(variantAltAllele.substring(1))) {
						alleleMatch = true;
						alleleMatchedAfterComplement = false;
					} else {
						alleleMatch = false;
						alleleMatchedAfterComplement = false;
					}
				}


			} else {
				alleleMatch = false;
				alleleMatchedAfterComplement = false;
			}

			final boolean probePassedQc;

			if (probeMatchedCounter.get(snpName) == 1
					&& !g1000SnpInProbeEnd
					&& !gonlSnpInProbeEnd
					&& mismatchesProbeEnd == 0
					&& alleleMatch
					&& !probe1000gOverlap
					&& !probeGonlOverlap) {
				probePassedQc = true;
			} else {
				probePassedQc = false;
			}


			final String AAllele;
			final String BAllele;
			if (bialleleicVariant) {
				if (isSnpVariant) {

					if ((variantAltAllele.equals("G") && variantRefAllele.equals("C"))
							|| (variantAltAllele.equals("C") && variantRefAllele.equals("G"))
							|| (variantAltAllele.equals("A") && variantRefAllele.equals("T"))
							|| (variantAltAllele.equals("T") && variantRefAllele.equals("A"))) {

						final int measurePos;
						if (reverseStrand) {
							measurePos = varPosition + 1;
						} else {
							measurePos = varPosition - 1;
						}
						final char refMeasureAllele = referenceGenome.getNucleotide(chr, measurePos);
						if (refMeasureAllele == 'A' || refMeasureAllele == 'T') {
							AAllele = variantRefAllele;
							BAllele = variantAltAllele;
						} else {
							AAllele = variantAltAllele;
							BAllele = variantRefAllele;
						}

					} else {
						if (variantRefAllele.equals("T") || variantRefAllele.equals("A")) {
							AAllele = variantRefAllele;
							BAllele = variantAltAllele;
						} else {
							AAllele = variantAltAllele;
							BAllele = variantRefAllele;
						}
					}


				} else {
					//indels
					final char refMeasureAllele = referenceGenome.getNucleotide(chr, varPosition);
					if (refMeasureAllele == 'A' || refMeasureAllele == 'T') {
						if (variantRefAllele.length() > 1) {
							AAllele = "I";
							BAllele = "D";
						} else {
							AAllele = "D";
							BAllele = "I";
						}
					} else {
						if (variantRefAllele.length() > 1) {
							AAllele = "D";
							BAllele = "I";
						} else {
							AAllele = "I";
							BAllele = "D";
						}
					}
				}
			} else {
				AAllele = "";
				BAllele = "";
			}


			if (chr.equals("X") && (ParB37.PAR1.isInChrXRange(varPosition) || ParB37.PAR2.isInChrXRange(varPosition))) {
				chr = "XY";
			}

			if (!isSnpVariant) {
				if (pos1000gVariant > 0) {
					varPosition = pos1000gVariant;
				} else if (posGonlVariant > 0) {
					varPosition = posGonlVariant;
				}
			}

			i = 0;
			mappingReportEntry[i++] = snpName;
			mappingReportEntry[i++] = chr;
			mappingReportEntry[i++] = Integer.toString(varPosition);
			mappingReportEntry[i++] = refAllele;
			mappingReportEntry[i++] = Integer.toString(leftProbePos);
			mappingReportEntry[i++] = reverseStrand ? "-" : "+";
			mappingReportEntry[i++] = Integer.toString(probeLength);
			mappingReportEntry[i++] = record.getCigarString();
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
			mappingReportEntry[i++] = String.valueOf(pos1000gVariant);
			mappingReportEntry[i++] = String.valueOf(af1000g);
			mappingReportEntry[i++] = String.valueOf(alleleCount1000g);
			mappingReportEntry[i++] = String.valueOf(snp);
			mappingReportEntry[i++] = String.valueOf(multipleVarAtPos);
			mappingReportEntry[i++] = String.valueOf(probe1000gOverlap);
			mappingReportEntry[i++] = String.valueOf(g1000SnpInProbeEnd);
			mappingReportEntry[i++] = refAlleleGonl;
			mappingReportEntry[i++] = altAlleleGonl;
			mappingReportEntry[i++] = String.valueOf(posGonlVariant);
			mappingReportEntry[i++] = String.valueOf(afGonl);
			mappingReportEntry[i++] = String.valueOf(alleleCountGonl);
			mappingReportEntry[i++] = String.valueOf(snpGonl);
			mappingReportEntry[i++] = String.valueOf(multipleVarAtPosGonl);
			mappingReportEntry[i++] = String.valueOf(probeGonlOverlap);
			mappingReportEntry[i++] = String.valueOf(gonlSnpInProbeEnd);
			mappingReportEntry[i++] = variantRefAllele;
			mappingReportEntry[i++] = variantAltAllele;
			mappingReportEntry[i++] = AAllele;
			mappingReportEntry[i++] = BAllele;
			mappingReportEntry[i++] = String.valueOf(bialleleicVariant);
			mappingReportEntry[i++] = String.valueOf(isSnpVariant);
			mappingReportEntry[i++] = String.valueOf(alleleMatch);
			mappingReportEntry[i++] = String.valueOf(alleleMatchedAfterComplement);
			mappingReportEntry[i++] = String.valueOf(probePassedQc);

			mappingReportWriter.writeNext(mappingReportEntry);

			if (probePassedQc) {
				i = 0;
				passedProbeEntry[i++] = snpName;
				passedProbeEntry[i++] = chr;
				passedProbeEntry[i++] = Integer.toString(varPosition);
				passedProbeEntry[i++] = variantRefAllele;
				passedProbeEntry[i++] = variantAltAllele;
				passedProbeEntry[i++] = AAllele;
				passedProbeEntry[i++] = BAllele;
				passedProbesWriter.writeNext(passedProbeEntry);
			} else if (!failReported.contains(snpName)) {



				failReported.add(snpName);

				FailReason reason;


				/*
				 * 
				 * probeMatchedCounter.get(snpName) == 1
				 && !g1000SnpInProbeEnd
				 && !gonlSnpInProbeEnd
				 && mismatchesProbeEnd == 0
				 && alleleMatch
				 && !probe1000gOverlap
				 && !probeGonlOverlap) {
				 * 
				 */
				if (probeMatchedCounter.get(snpName) > 1) {
					reason = FailReason.MULTIPLE_MAPPING;
				} else if (alleleCount1000g == 0 && alleleCountGonl == 0) {
					reason = FailReason.NO_VAR_REFERENCE;
				} else if (!bialleleicVariant) {
					reason = FailReason.NOT_BIALLELIC_REFERENCE;
				} else if (!alleleMatch) {
					reason = FailReason.INCONSISTENT_ALLELES;
				} else if (mismatchesProbeEnd >= 1) {
					reason = FailReason.MISMATCH_TO_REFERENCE_IN_END;
				} else if (g1000SnpInProbeEnd || gonlSnpInProbeEnd) {
					reason = FailReason.VARIATION_PROBE_END;
				} else if (probe1000gOverlap || probeGonlOverlap) {
					reason = FailReason.OVERLAPPING_INDEL;
				} else {
					mappingReportWriter.close();
					throw new Exception("Probe " + snpName + " failed without defined reason. (This could be a bug in the software)");
				}

				failCounter.adjustOrPutValue(reason, 1, 1);
				i = 0;
				failedProbeEntry[i++] = snpName;
				failedProbeEntry[i++] = reason.getReason();
				failedProbesWriter.writeNext(failedProbeEntry);

			}
		}

		mappingReportWriter.close();
		passedProbesWriter.close();

		for (String probe : probesInManifestNotMapped) {

			failCounter.adjustOrPutValue(FailReason.NO_MAPPING, 1, 1);
			i = 0;
			failedProbeEntry[i++] = probe;
			failedProbeEntry[i++] = FailReason.NO_MAPPING.getReason();
			failedProbesWriter.writeNext(failedProbeEntry);

		}

		failedProbesWriter.close();

		System.out.println(
				"Failed probe summary:");
		for (FailReason failedReason : FailReason.values()) {
			System.out.println(" - " + failedReason.getReason() + ": " + failCounter.get(failedReason));
		}

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

			if (SNP.equals("[N/A]")) {
				a1 = '\0';
				a2 = '\0';
				a1_b = "";
				a2_b = "";
			} else if (SNP.equals("")) {
				a1 = '\0';
				a2 = '\0';
				a1_b = "";
				a2_b = "";
			} else {
				a1 = SNP.charAt(1);
				a2 = SNP.charAt(3);

				Matcher alleleMatcher = ALLELE_PATTERN.matcher(SourceSeq);
				alleleMatcher.find();
				a1_b = alleleMatcher.group(1);
				a2_b = alleleMatcher.group(2);
			}


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

		private boolean hasNoAlleles() {
			return a1 == '\0' || a2 == '\0';
		}
	}

	private static enum ParB37 {

		/*
		 * http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/
		 * 
		 * HG37
		 * 
		 * Name 	Chr 	Start 	Stop
		 * PAR#1 	X 	60,001 	2,699,520
		 * PAR#2 	X 	154,931,044 	155,260,560
		 * PAR#1 	Y 	10,001 	2,649,520
		 * PAR#2 	Y 	59,034,050 	59,363,566
		 */
		PAR1(60001, 2699520, 10001, 2649520),
		PAR2(154931044, 155260560, 59034050, 59363566);
		private final int startX;
		private final int stopX;
		private final int startY;
		private final int stopY;

		ParB37(int startX, int stopX, int startY, int stopY) {
			this.startX = startX;
			this.stopX = stopX;
			this.startY = startY;
			this.stopY = stopY;
		}

		public boolean isInChrXRange(int pos) {
			return pos >= this.startX && pos < this.stopX;
		}

		public boolean isInChrYRange(int pos) {
			return pos >= this.startY && pos < this.stopY;
		}

		public int posInY(int posX) {
			if (!isInChrXRange(posX)) {
				throw new IllegalArgumentException();
			}

			return startY + posX - startX;

		}

		public int posInX(int posY) {
			if (!isInChrYRange(posY)) {
				throw new IllegalArgumentException();
			}

			return startX + posY - startY;

		}
	}

	private static enum FailReason {

		NO_MAPPING("Not mapped"),
		MULTIPLE_MAPPING("Mapped to multiple locations"),
		MISMATCH_TO_REFERENCE_IN_END("Missmatch in probe end"),
		OVERLAPPING_INDEL("Indel / rearrangement overlapping variant MAF >= 1%"),
		NO_VAR_REFERENCE("Variant not found in reference"),
		NOT_BIALLELIC_REFERENCE("Variant is not biallelic"),
		VARIATION_PROBE_END("Reference variants in probe end MAF >= 1%"),
		ILLUMINA_EXCLUDED("Not variant in Illumina manifest"),
		INCONSISTENT_ALLELES("Unexpected alleles in reference");
		private final String reason;

		FailReason(String reason) {
			this.reason = reason;
		}

		public String getReason() {
			return reason;
		}
	}
}
