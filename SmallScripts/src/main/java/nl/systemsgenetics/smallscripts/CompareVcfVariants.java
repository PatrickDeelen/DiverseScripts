package nl.systemsgenetics.smallscripts;

import java.text.NumberFormat;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilterSeq;

/**
 * Hello world!
 *
 */
public class CompareVcfVariants 
{
	
	public static final NumberFormat DEFAULT_NUMBER_FORMATTER = NumberFormat.getInstance();
	
	
    public static void main( String[] args ) throws Exception
    {
		
		VariantFilterSeq variantFilterSeq = new VariantFilterSeq("2");
		
		final String vcf1Path = "D:\\UMCG\\Genetica\\Projects\\RnaSeqEqtl\\ALL.phase1_release_v3.20101123.snps.passed.maf0.01.clinvarAdded.sorted.filteredRepRNAedJunc.vcf.gz";
		final String vcf2Path = "D:\\UMCG\\Genetica\\Projects\\RnaSeqEqtl\\GoNL+1000G+ClinVar+UMCG_allSNPs.filtered.vcf.gz";
		
		RandomAccessGenotypeData vcf1 = RandomAccessGenotypeDataReaderFormats.VCF.createFilteredGenotypeData(vcf1Path, 100, variantFilterSeq, null);
		RandomAccessGenotypeData vcf2 = RandomAccessGenotypeDataReaderFormats.VCF.createGenotypeData(vcf2Path);
		
		int vcf1VariantCount = 0;
		int vcf1VariantFoundInVcf2 = 0;
		
		for(GeneticVariant vcf1Variant : vcf1){
			
			++vcf1VariantCount;
			
			if(vcf1VariantCount % 100000 == 0){
				System.out.println("Parsed: " + DEFAULT_NUMBER_FORMATTER.format(vcf1VariantCount) + " variants from VCF1");
			}
			
			if(vcf2.getSnpVariantByPos(vcf1Variant.getSequenceName(), vcf1Variant.getStartPos()) != null){
				++vcf1VariantFoundInVcf2;
			}
		}
		
		System.out.println("Total variants in VCF1 after filter: " + DEFAULT_NUMBER_FORMATTER.format(vcf1VariantCount));
		System.out.println("Variants from VCF1 found in VCF2: " + DEFAULT_NUMBER_FORMATTER.format(vcf1VariantFoundInVcf2));
		
    }
}
