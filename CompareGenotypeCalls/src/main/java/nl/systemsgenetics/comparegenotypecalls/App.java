package nl.systemsgenetics.comparegenotypecalls;

import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Pattern;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.sampleFilter.SampleFilter;
import org.molgenis.genotype.sampleFilter.SampleIdIncludeFilter;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilterSeqPos;
import org.molgenis.genotype.variantFilter.VariantFilterableGenotypeDataDecorator;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import org.molgenis.genotype.variantFilter.VariantQcChecker;

/**
 * Hello world!
 *
 */
public class App {

	static Pattern TAB_PATTERN = Pattern.compile("\\t");

	public static void main(String[] args) throws Exception {

		
		String data1Type = args[0];
		String data1Path = args[1];
		String data2Type = args[2];
		String data2Path = args[3];
		String idMapPath = args[4];
		String outputFilePath = args[5];
		String mafFilterData2 = args.length >= 7 ? args[6] : "0";
		String snpIdFilterData1 = args.length >= 8 && !args[7].equals("null") ? args[7] : null;
		String data1ProbCall = args.length == 9 ? args[8] : "0.4";
		
		System.out.println("Data1 " + data1Type + " at " + data1Path);
		System.out.println("Data2 " + data2Type + " at " + data2Path);
		System.out.println("SampleIdMap (data1Id\tdata2Id): " + idMapPath);
		System.out.println("output file: " + outputFilePath);
		System.out.println("Maf filter data 2: " + mafFilterData2);
		if(snpIdFilterData1 != null){
			System.out.println("SNP ID filter data 1: " + snpIdFilterData1);
		}
		System.out.println("Data 1 prob call: " + data1ProbCall);
		

		HashMap<String, String> sampleIdMap = new HashMap<String, String>();
		BufferedReader sampleMapReader = new BufferedReader(new FileReader(idMapPath));
		String line;
		String[] elements;
		while ((line = sampleMapReader.readLine()) != null) {
			elements = TAB_PATTERN.split(line);
			sampleIdMap.put(elements[0], elements[1]);
		}

		VariantIdIncludeFilter snpIdFilter = null;
		if(snpIdFilterData1 != null){
			HashSet<String> snps = new HashSet<String>();
			BufferedReader snpIdFilterReader = new BufferedReader(new FileReader(snpIdFilterData1));
			while ((line = snpIdFilterReader.readLine()) != null) {
				snps.add(line);
			}
			snpIdFilter = new VariantIdIncludeFilter(snps);
		}
		

		SampleFilter data1SampleFilter = new SampleIdIncludeFilter(sampleIdMap.keySet());
		SampleFilter data2SampleFilter = new SampleIdIncludeFilter(sampleIdMap.values());

		
		
		RandomAccessGenotypeData data1 = RandomAccessGenotypeDataReaderFormats.valueOf(data1Type.toUpperCase()).createFilteredGenotypeData(data1Path, 1024, snpIdFilter, data1SampleFilter, null, Double.parseDouble(data1ProbCall));

		VariantFilterSeqPos seqPosFilter = new VariantFilterSeqPos();
		for (GeneticVariant data1Var : data1) {
			seqPosFilter.addSeqPos(data1Var);
		}

		RandomAccessGenotypeData data2 = RandomAccessGenotypeDataReaderFormats.valueOf(data2Type.toUpperCase()).createFilteredGenotypeData(data2Path, 1024, seqPosFilter, null);

		//Do here to optimize trityper 
		data2 = new VariantFilterableGenotypeDataDecorator(data2, new VariantQcChecker(Float.valueOf(mafFilterData2), 0, 0));
		
		
		System.out.println("Reading data 2 samples");

		String[] data2SamplesNames = data2.getSampleNames();
		TObjectIntHashMap<String> data2SampleIndex = new TObjectIntHashMap<String>();
		for (int i = 0; i < data2SamplesNames.length; ++i) {
			data2SampleIndex.put(data2SamplesNames[i], i);
			//System.out.println(data2SamplesNames[i] + " " + i);
		}

		System.out.println("Reading data 1 samples");

		String[] data1SamplesNames = data1.getSampleNames();
		TObjectIntHashMap<String> data1SampleIndex = new TObjectIntHashMap<String>();
		TObjectIntHashMap sharedSamplesIndex = new TObjectIntHashMap<String>();
		ArrayList<String> sharedSamples = new ArrayList<String>();
		for (int i = 0; i < data1SamplesNames.length; ++i) {
			data1SampleIndex.put(data1SamplesNames[i], i);
			//System.out.println(data1SamplesNames[i] + " " + i);
			if (sampleIdMap.containsKey(data1SamplesNames[i]) && data2SampleIndex.containsKey(sampleIdMap.get(data1SamplesNames[i]))) {
				sharedSamplesIndex.put(data1SamplesNames[i], sharedSamples.size());
				sharedSamples.add(data1SamplesNames[i]);
			}
		}

		if (sharedSamples.size() != sampleIdMap.size()) {
			System.out.println("Expected " + sampleIdMap.size() + " found in mapping but only " + sharedSamples.size() + " are found te be shared");
		}

		int[] sampleIdenticalGenotypeCount = new int[sharedSamples.size()];
		int[] sampleTestedGenotypeCount = new int[sharedSamples.size()];

		List<Alleles> data1VarAlleles;
		List<Alleles> data2VarAlleles;

		int i = 1;
		int skippedVar = 0;

		for (GeneticVariant data1Var : data1) {

			if (i % 1000 == 0) {
				System.out.println("Variant: " + i);
			}
			++i;

			GeneticVariant data2Var;
			if ((data2Var = data2.getSnpVariantByPos(data1Var.getSequenceName(), data1Var.getStartPos())) != null) {

				if(!data1Var.getVariantAlleles().sameAlleles(data2Var.getVariantAlleles())){
					System.err.println("Different alleles for " + data1Var.getPrimaryVariantId() + " " + data1Var.getVariantAlleles() + " vs " + data2Var.getVariantAlleles() + ". " + data1Var.getSequenceName() + ":" + data1Var.getStartPos() + " vs " + data2Var.getSequenceName() + ":" + data2Var.getStartPos() );
					
					++skippedVar;
					
					continue;
				}
				
				data1VarAlleles = data1Var.getSampleVariants();
				data2VarAlleles = data2Var.getSampleVariants();

				for (String sample : sharedSamples) {

					Alleles data1SampleAlleles = data1VarAlleles.get(data1SampleIndex.get(sample));
					Alleles data2SampleAlleles = data2VarAlleles.get(data2SampleIndex.get(sampleIdMap.get(sample)));
					
					//System.out.println(sample + " " + sampleIdMap.get(sample));
					
					//System.out.println("   data1 index: " + data1SampleIndex.get(sample));
					//System.out.println("   data2 index: " + data2SampleIndex.get(sampleIdMap.get(sample)));

					if (data1SampleAlleles.contains(Allele.ZERO) || data2SampleAlleles.contains(Allele.ZERO)) {
						continue;
					}

					int sampleIndex = sharedSamplesIndex.get(sample);

					sampleTestedGenotypeCount[sampleIndex]++;

					if (data1SampleAlleles.sameAlleles(data2SampleAlleles)) {
						
						//System.out.println(data1SampleAlleles + " vs " + data2SampleAlleles);
						sampleIdenticalGenotypeCount[sampleIndex]++;
					} 

				}

			}
			
			//break;

		}
		
		System.out.println("Skipped vars due to incosistant alleles: " + skippedVar);

		BufferedWriter out = new BufferedWriter(new FileWriter(outputFilePath));

		out.append("sample");
		out.append('\t');
		out.append("sampleTestedCount");
		out.append('\t');
		out.append("sampleIdenticalCount");
		out.append('\t');
		out.append("identical");
		out.append('\n');

		for (String sample : sharedSamples) {

			int sampleIndex = sharedSamplesIndex.get(sample);
			int sampleTestedCount = sampleTestedGenotypeCount[sampleIndex];
			int sampleIdenticalCount = sampleIdenticalGenotypeCount[sampleIndex];
			double identical = (sampleIdenticalCount * 100d) / sampleTestedCount;

			out.append(sample);
			out.append('\t');
			out.append(Integer.toString(sampleTestedCount));
			out.append('\t');
			out.append(Integer.toString(sampleIdenticalCount));
			out.append('\t');
			out.append(Double.toString(identical));
			out.append('\n');

		}

		out.close();

	}
}
