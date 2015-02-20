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
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.math3.stat.regression.SimpleRegression;
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
public class App_1 {

	static final Pattern TAB_PATTERN = Pattern.compile("\\t");
	static final String DEFAULT_CALL_P = "0.7";
	
	private static final Options OPTIONS;

	private static final String HEADER =
			"  /---------------------------------------\\\n"
			+ "  |         Compare genotype calls        |\n"
			+ "  |                                       |\n"
			+ "  |             Patrick Deelen            |\n"
			+ "  |        patrickdeelen@gmail.com        |\n"
			+ "  |                                       |\n"
			+ "  |     Genomics Coordication Center      |\n"
			+ "  |        Department of Genetics         |\n"
			+ "  |  University Medical Center Groningen  |\n"
			+ "  \\---------------------------------------/";
	
	static {

		OPTIONS = new Options();
		
		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Gentoype data 1");
		OptionBuilder.withLongOpt("data1");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("d1"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The input data 1 type.\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* GENFOLDER - matches all Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("data1Type");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("D1"));
		
		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Gentoype data 2");
		OptionBuilder.withLongOpt("data2");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("d2"));

		OptionBuilder.withArgName("type");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("The input data 1 type.\n"
				+ "* PED_MAP - plink PED MAP files.\n"
				+ "* PLINK_BED - plink BED BIM FAM files.\n"
				+ "* VCF - bgziped vcf with tabix index file\n"
				+ "* VCFFOLDER - matches all bgziped vcf files + tabix index in a folder\n"
				+ "* SHAPEIT2 - shapeit2 phased haplotypes .haps & .sample\n"
				+ "* GEN - Oxford .gen & .sample\n"
				+ "* GENFOLDER - matches all Oxford .gen & .sample\n"
				+ "* TRITYPER - TriTyper format folder");
		OptionBuilder.withLongOpt("data2Type");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("D2"));
		
		OptionBuilder.withArgName("basePath");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Output base path");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("o"));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Sample ID mappings (data1Id\tdata2Id)");
		OptionBuilder.withLongOpt("sampleMap");
		OptionBuilder.isRequired();
		OPTIONS.addOption(OptionBuilder.create("s"));
		
		OptionBuilder.withArgName("double");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Min posterior probability to make a call in data 1. Default: " + DEFAULT_CALL_P);
		OptionBuilder.withLongOpt("prob1");
		OPTIONS.addOption(OptionBuilder.create("p1"));
		
		OptionBuilder.withArgName("double");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Min posterior probability to make a call in data 2. Default: " + DEFAULT_CALL_P);
		OptionBuilder.withLongOpt("prob2");
		OPTIONS.addOption(OptionBuilder.create("p2"));
		
		OptionBuilder.withArgName("path");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Variant include filter file data 1");
		OptionBuilder.withLongOpt("variantInclude");
		OPTIONS.addOption(OptionBuilder.create("v"));
		
		OptionBuilder.withArgName("string");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Force chr");
		OptionBuilder.withLongOpt("chr");
		OPTIONS.addOption(OptionBuilder.create("c"));
		
		OptionBuilder.withArgName("double");
		OptionBuilder.hasArgs();
		OptionBuilder.withDescription("Maf filter data 2");
		OptionBuilder.withLongOpt("maf2");
		OPTIONS.addOption(OptionBuilder.create("m2"));
		
	}


	public static void main(String[] args) throws Exception {

		
		System.out.println(HEADER);
		System.out.println();
		System.out.flush(); //flush to make sure header is before errors
		try {
			Thread.sleep(25); //Allows flush to complete
		} catch (InterruptedException ex) {
		}
		

		final String data1Type;
		final String data1Path;
		final String data2Type;
		final String data2Path;
		final String idMapPath;
		final String outputFilePath;
		final String mafFilterData2;
		final String snpIdFilterData1;
		final String data1ProbCall;
		final String data2ProbCall;
		final String forceChr;
		
		try {
			final CommandLine commandLine = new PosixParser().parse(OPTIONS, args, false);

			data1Type = commandLine.getOptionValue("D1");
			data1Path = commandLine.getOptionValue("d1");
			data2Type = commandLine.getOptionValue("D2");
			data2Path = commandLine.getOptionValue("d2");
			idMapPath = commandLine.getOptionValue("s");
			outputFilePath = commandLine.getOptionValue("o");
			mafFilterData2 = commandLine.getOptionValue("m2", "0");
			snpIdFilterData1 = commandLine.getOptionValue("v", null);
			data1ProbCall = commandLine.getOptionValue("p1", DEFAULT_CALL_P);
			data2ProbCall = commandLine.getOptionValue("p2", DEFAULT_CALL_P);
			forceChr = commandLine.getOptionValue("c", null);

		} catch (ParseException ex) {
			System.err.println("Invalid command line arguments: ");
			System.err.println(ex.getMessage());
			System.err.println();
			new HelpFormatter().printHelp(" ", OPTIONS);
			System.exit(1);
			return;
		}

		System.out.println("Data1 " + data1Type + " at " + data1Path);
		System.out.println("Data2 " + data2Type + " at " + data2Path);
		System.out.println("SampleIdMap (data1Id\tdata2Id): " + idMapPath);
		System.out.println("output file: " + outputFilePath);
		System.out.println("Maf filter data 2: " + mafFilterData2);
		if (snpIdFilterData1 != null) {
			System.out.println("SNP ID filter data 1: " + snpIdFilterData1);
		}
		System.out.println("Data 1 prob call: " + data1ProbCall);
		System.out.println("Data 2 prob call: " + data2ProbCall);
		if (forceChr != null) {
			System.out.println("Force Chr: " + forceChr);
		}

		HashMap<String, String> sampleIdMap = new HashMap<String, String>();
		BufferedReader sampleMapReader = new BufferedReader(new FileReader(idMapPath));
		String line;
		String[] elements;
		while ((line = sampleMapReader.readLine()) != null) {
			elements = TAB_PATTERN.split(line);
			sampleIdMap.put(elements[0], elements[1]);
		}
		
		SampleFilter data1SampleFilter = new SampleIdIncludeFilter(sampleIdMap.keySet());
		SampleFilter data2SampleFilter = new SampleIdIncludeFilter(sampleIdMap.values());

		VariantIdIncludeFilter snpIdFilter = null;
		if (snpIdFilterData1 != null) {
			HashSet<String> snps = new HashSet<String>();
			BufferedReader snpIdFilterReader = new BufferedReader(new FileReader(snpIdFilterData1));
			while ((line = snpIdFilterReader.readLine()) != null) {
				snps.add(line);
			}
			snpIdFilter = new VariantIdIncludeFilter(snps);
		}


		
		



		RandomAccessGenotypeData data1 = RandomAccessGenotypeDataReaderFormats.valueOfSmart(data1Type.toUpperCase()).createFilteredGenotypeData(data1Path, 1024, snpIdFilter, data1SampleFilter, forceChr, Double.parseDouble(data1ProbCall));

		VariantFilterSeqPos seqPosFilter = new VariantFilterSeqPos();
		for (GeneticVariant data1Var : data1) {
			seqPosFilter.addSeqPos(data1Var);
		}

		
		
		RandomAccessGenotypeData data2 = RandomAccessGenotypeDataReaderFormats.valueOfSmart(data2Type.toUpperCase()).createFilteredGenotypeData(data2Path, 1024, seqPosFilter, data2SampleFilter, forceChr, Double.parseDouble(data2ProbCall));

		//Do here to optimize trityper 
		data2 = new VariantFilterableGenotypeDataDecorator(data2, new VariantQcChecker(Float.valueOf(mafFilterData2), 0, 0));


		System.out.println("Reading data 2 samples, total: " + data2.getSampleNames().length);

		String[] data2SamplesNames = data2.getSampleNames();
		TObjectIntHashMap<String> data2SampleIndex = new TObjectIntHashMap<String>();
		for (int i = 0; i < data2SamplesNames.length; ++i) {
			data2SampleIndex.put(data2SamplesNames[i], i);
			//System.out.println(data2SamplesNames[i] + " " + i);
		}

		System.out.println("Reading data 1 samples, total: " + data1.getSampleNames().length);

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
			System.out.println("");
		}

		float[] data1VarDosages;
		float[] data2VarDosages;

		List<Alleles> data1VarAlleles;
		List<Alleles> data2VarAlleles;

		int i = 0;
		int skippedVar = 0;
		int comparedVar = 0;


		BufferedWriter outSnp = new BufferedWriter(new FileWriter(outputFilePath + ".snp"));
		outSnp.append("snp\tchr\tpos\talleles\tr2\tidenticalCall\tmafData2\tsampleCount\tcallData1\tcallData2\n");

		HashMap<String, SimpleRegression> sampleCor = new HashMap<String, SimpleRegression>(sharedSamples.size());
		HashMap<String, AtomicInteger> sampleVarCount = new HashMap<String, AtomicInteger>(sharedSamples.size());
		HashMap<String, AtomicInteger> sampleData1CallCount = new HashMap<String, AtomicInteger>(sharedSamples.size());
		HashMap<String, AtomicInteger> sampleData2CallCount = new HashMap<String, AtomicInteger>(sharedSamples.size());
		HashMap<String, AtomicInteger> sampleIndenticalCallCount = new HashMap<String, AtomicInteger>(sharedSamples.size());

		for (String sample : sharedSamples) {
			sampleCor.put(sample, new SimpleRegression());
			sampleVarCount.put(sample, new AtomicInteger());
			sampleData1CallCount.put(sample, new AtomicInteger());
			sampleData2CallCount.put(sample, new AtomicInteger());
			sampleIndenticalCallCount.put(sample, new AtomicInteger());
		}


		for (GeneticVariant data1Var : data1) {

			++i;
			if (i % 1000 == 0) {
				System.out.println("Variant: " + i);
			}

			GeneticVariant data2Var;
			if ((data2Var = data2.getSnpVariantByPos(data1Var.getSequenceName(), data1Var.getStartPos())) != null) {


				if (!data1Var.getVariantAlleles().sameAlleles(data2Var.getVariantAlleles())) {
					System.err.println("Different alleles for " + data1Var.getPrimaryVariantId() + " " + data1Var.getVariantAlleles() + " vs " + data2Var.getVariantAlleles() + " " + data1Var.getSequenceName() + ":" + data1Var.getStartPos() + " vs " + data2Var.getSequenceName() + ":" + data2Var.getStartPos());

					++skippedVar;

					continue;
				}

				++comparedVar;

				data1VarDosages = data1Var.getSampleDosages();
				data2VarDosages = data2Var.getSampleDosages();

				data1VarAlleles = data1Var.getSampleVariants();
				data2VarAlleles = data2Var.getSampleVariants();

				if (data1Var.getVariantAlleles().get(0) != data2Var.getVariantAlleles().get(0)) {
					for (int j = 0; j < data2VarDosages.length; ++j) {
						data2VarDosages[j] = (data2VarDosages[j] - 2) * -1;
					}
				}

				SimpleRegression varCor = new SimpleRegression();

				int snpData1CallCount = 0;
				int snpData2CallCount = 0;
				int snpIndenticalCall = 0;

				for (String sample : sharedSamples) {

					double data1SampleDosage = data1VarDosages[data1SampleIndex.get(sample)];
					double data2SampleDosage = data2VarDosages[data2SampleIndex.get(sampleIdMap.get(sample))];

					Alleles data1SampleAllele = data1VarAlleles.get(data1SampleIndex.get(sample));
					Alleles data2SampleAllele = data2VarAlleles.get(data2SampleIndex.get(sampleIdMap.get(sample)));

					boolean data1Missing = data1SampleAllele.contains(Allele.ZERO) || data1SampleAllele.getAlleleCount() == 0 || data1SampleDosage < 0;
					boolean data2Missing = data2SampleAllele.contains(Allele.ZERO) || data2SampleAllele.getAlleleCount() == 0 || data2SampleDosage < 0;

					if (!data1Missing) {
						++snpData1CallCount;
						sampleData1CallCount.get(sample).getAndIncrement();
					}

					if (!data2Missing) {
						++snpData2CallCount;
						sampleData2CallCount.get(sample).getAndIncrement();
					}

					if (data1Missing || data2Missing) {
						continue;
					}

					if (data1SampleAllele.sameAlleles(data2SampleAllele)) {
						++snpIndenticalCall;
						sampleIndenticalCallCount.get(sample).getAndIncrement();
					}

					sampleVarCount.get(sample).getAndIncrement();
					varCor.addData(data1SampleDosage, data2SampleDosage);
					sampleCor.get(sample).addData(data1SampleDosage, data2SampleDosage);

				}


				outSnp.append(data1Var.getPrimaryVariantId());
				outSnp.append('\t');
				outSnp.append(data1Var.getSequenceName());
				outSnp.append('\t');
				outSnp.append(String.valueOf(data1Var.getStartPos()));
				outSnp.append('\t');
				outSnp.append(data1Var.getVariantAlleles().toString());
				outSnp.append('\t');
				outSnp.append(String.valueOf(varCor.getRSquare()));
				outSnp.append('\t');
				outSnp.append(String.valueOf(snpIndenticalCall / (double) varCor.getN()));
				outSnp.append('\t');
				outSnp.append(String.valueOf(data2Var.getMinorAlleleFrequency()));
				outSnp.append('\t');
				outSnp.append(String.valueOf(varCor.getN()));
				outSnp.append('\t');
				outSnp.append(String.valueOf(snpData1CallCount / (double) sharedSamples.size()));
				outSnp.append('\t');
				outSnp.append(String.valueOf(snpData2CallCount / (double) sharedSamples.size()));
				outSnp.append('\n');



			}

			//break;

		}

		outSnp.close();

		System.out.println("Number of SNPs compared: " + comparedVar);
		System.out.println("Skipped vars due to incosistant alleles: " + skippedVar);

		BufferedWriter out = new BufferedWriter(new FileWriter(outputFilePath + ".sample"));

		out.append("sample");
		out.append('\t');
		out.append("R2");
		out.append('\t');
		out.append("IdenticalCall");
		out.append('\t');
		out.append("VarCount");
		out.append('\t');
		out.append("CallData1");
		out.append('\t');
		out.append("CallData2");
		out.append('\n');

		for (String sample : sharedSamples) {

			out.append(sample);
			out.append('\t');
			out.append(Double.toString(sampleCor.get(sample).getRSquare()));
			out.append('\t');
			out.append(Double.toString(sampleIndenticalCallCount.get(sample).get() / (double) comparedVar));
			out.append('\t');
			out.append(sampleVarCount.get(sample).toString());
			out.append('\t');
			out.append(String.valueOf(sampleData1CallCount.get(sample).get() / (double) comparedVar));
			out.append('\t');
			out.append(String.valueOf(sampleData2CallCount.get(sample).get() / (double) comparedVar));
			out.append('\n');

		}

		out.close();

	}
}
