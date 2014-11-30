package nl.systemsgenetics.smallscripts;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Iterator;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;

/**
 *
 * @author Patrick Deelen
 */
public class FindEqtlMatchedRandom {

	public static void main(String[] args) throws Exception {

		System.out.println("eqtl file: " + args[0]);
		System.out.println("vcf folder genotypes: " + args[1]);
		System.out.println("result: " + args[2]);
		System.out.println("window " + args[3]);
		System.out.println("maf diff " + args[4]);
		System.out.println("min r2 diff " + args[5]);

		int window = Integer.parseInt(args[3]);
		double mafDiffPercentage = Double.parseDouble(args[4]);
		double r2Diff = Double.parseDouble(args[5]);


		QTLTextFile eQTLsTextFile = new QTLTextFile(args[0], false);
		RandomAccessGenotypeData genotypes = RandomAccessGenotypeDataReaderFormats.VCF_FOLDER.createFilteredGenotypeData(args[1], 1000, null, null);

		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(args[2])));
		writer.append("eQTL_SNP\teQTL_SNP_Chr_\teQTL_SNP_Pos\tControl_SNP\tControl_SNP_Chr\tControl_SNP_Pos\tControl_SNP_A1\tControl_SNP_A2\n");

		int eQtlsWithoutMatch = 0;

		for (Iterator<EQTL> eQtlIt = eQTLsTextFile.getEQtlIterator(); eQtlIt.hasNext();) {

			EQTL eQtl = eQtlIt.next();
			String chr = eQtl.getRsChr() + "";
			int pos = eQtl.getRsChrPos();

			GeneticVariant eQtlVariant = genotypes.getSnpVariantByPos(chr, pos);

			GeneticVariant controlVariant = null;

			double minMafOther = eQtlVariant.getMinorAlleleFrequency() - mafDiffPercentage;
			double maxMafOther = eQtlVariant.getMinorAlleleFrequency() + mafDiffPercentage;

			for (GeneticVariant otherVariant : genotypes.getVariantsByRange(chr, pos - window < 0 ? 0 : pos - window, pos + window)) {
				if (!otherVariant.isBiallelic() || !otherVariant.isSnp()){
					continue;
				}
				if (otherVariant.getMinorAlleleFrequency() < minMafOther || otherVariant.getMinorAlleleFrequency() > maxMafOther) {
					continue;
				}
				if (eQtlVariant.calculateLd(otherVariant).getR2() >= r2Diff) {
					continue;
				}
				if (controlVariant == null) {
					controlVariant = otherVariant;
				} else if (Math.abs(otherVariant.getStartPos() - pos) < Math.abs(controlVariant.getStartPos() - pos)) {
					controlVariant = otherVariant;
				}
			}

			writer.append(eQtlVariant.getPrimaryVariantId());
			writer.append('\t');
			writer.append(eQtlVariant.getSequenceName());
			writer.append('\t');
			writer.append(String.valueOf(eQtlVariant.getStartPos()));

			if (controlVariant == null) {
				++eQtlsWithoutMatch;
				writer.append("\tNA\t0\t0\t0\t0\n");
			} else {
				//GeneticVariant randomControl = potentialRandomControlVariants.get(random.nextInt());
				writer.append('\t');
				if (controlVariant.getPrimaryVariantId() == null) {
					writer.append('.');
				} else {
					writer.append(controlVariant.getPrimaryVariantId());

				}
				writer.append('\t');
				writer.append(controlVariant.getSequenceName());
				writer.append('\t');
				writer.append(String.valueOf(controlVariant.getStartPos()));
				writer.append('\t');
				writer.append(controlVariant.getVariantAlleles().get(0).toString());
				writer.append('\t');
				writer.append(controlVariant.getVariantAlleles().get(1).toString());
				writer.append('\n');
			}


		}
		writer.close();
		System.out.println("eQTLs without matching control: " + eQtlsWithoutMatch);
	}
}
