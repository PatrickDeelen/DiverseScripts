package nl.systemsgenetics.smallscripts;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilterWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.eQTLTextFile;

/**
 *
 * @author Patrick Deelen
 */
public class FindEqtlMatchedRandom {

	public static void main(String[] args) throws Exception {
		
		System.out.println("eqtl file: " + args[0]);
		System.out.println("trityper genotypes: " + args[1]);
		System.out.println("result: " + args[2]);
		System.out.println("window" + args[3]);
		System.out.println("maf diff %" + args[4]);
		System.out.println("random seed" + args[5]);
		
		int window = Integer.parseInt(args[3]);
		double mafDiffPercentage = Double.parseDouble(args[4]);
		long randomSeed = Long.parseLong(args[5]);
		
		Random random = new Random(randomSeed);
		
		eQTLTextFile eQTLsTextFile = new eQTLTextFile(args[0], false);
		RandomAccessGenotypeData genotypes = RandomAccessGenotypeDataReaderFormats.TRITYPER.createFilteredGenotypeData(args[1], 1000, null, null);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(args[2])));
		writer.append("eQTL_SNP\tControl_SNP\tControl_SNP_Chr\tControl_SNP_Pos\n");
		
		int eQtlsWithoutMatch = 0;
		
		for (Iterator<EQTL> eQtlIt = eQTLsTextFile.getEQtlIterator(); eQtlIt.hasNext();) {
			
			EQTL eQtl = eQtlIt.next();
			String chr = eQtl.getRsChr() + "";
			int pos = eQtl.getRsChrPos();
			
			GeneticVariant eQtlVariant = genotypes.getSnpVariantByPos(chr, pos);
			
			ArrayList<GeneticVariant> potentialRandomControlVariants = new ArrayList<GeneticVariant>();
			
			double minMafOther = eQtlVariant.getMinorAlleleFrequency() - mafDiffPercentage;
			double maxMafOther = eQtlVariant.getMinorAlleleFrequency() + mafDiffPercentage;
			
			for(GeneticVariant otherVariant : genotypes.getVariantsByRange(chr, pos - window, pos + window)){
				if(otherVariant.getMinorAlleleFrequency() < minMafOther || otherVariant.getMinorAlleleFrequency() > maxMafOther){
					continue;
				}
				if(eQtlVariant.calculateLd(otherVariant).getR2() >= 0.1){
					continue;
				}
				potentialRandomControlVariants.add(otherVariant);
			}
			
			writer.append(eQtlVariant.getPrimaryVariantId());
			if(potentialRandomControlVariants.isEmpty()){
				++eQtlsWithoutMatch;
				writer.append("\tNA\t0\t0\n");
			} else {
				GeneticVariant randomControl = potentialRandomControlVariants.get(random.nextInt());
				writer.append('\t');
				writer.append(randomControl.getPrimaryVariantId());
				writer.append('\t');
				writer.append(randomControl.getSequenceName());
				writer.append('\t');
				writer.append(String.valueOf(randomControl.getStartPos()));
				writer.append('\n');
			}
			
			
		}
		
		System.out.println("eQTLs without matching control: " + eQtlsWithoutMatch);
	}
	
	
	
}
