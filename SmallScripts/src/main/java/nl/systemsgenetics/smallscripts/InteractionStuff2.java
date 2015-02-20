package nl.systemsgenetics.smallscripts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author Patrick Deelen
 */
public class InteractionStuff2 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
        
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream("D:\\UMCG\\Genetica\\Projects\\rp3\\interactions\\interactionMetaZ8.txt"), "UTF-8"));
		
		HashMap<VariantGene, Interaction> interactions = new HashMap<VariantGene, Interaction>();
		HashMap<String, Interaction> covariateInteractions = new HashMap<String, Interaction>();

		reader.readLine();
		String line;
		while ((line = reader.readLine()) != null) {
			String[] elements = StringUtils.split(line, '\t');

			final String variant = elements[0];
			final String gene = elements[1];
			final String covarivate = elements[2];
			final double zscore = Double.parseDouble(elements[3]);
			
			VariantGene variantGene = new VariantGene(variant, gene);
			
			Interaction covariateInteraction = covariateInteractions.get(covarivate);
			Interaction interaction = interactions.get(variantGene);
			
			if(interaction == null){
				interaction = new Interaction();
				interactions.put(variantGene, interaction);
			}
			if(covariateInteraction == null){
				covariateInteraction = new Interaction();
				covariateInteractions.put(covarivate, covariateInteraction);
			}
			if(zscore < 0){
				interaction.addNegative(covarivate);
				covariateInteraction.addNegative(gene);
			} else {
				interaction.addPositive(covarivate);
				covariateInteraction.addPositive(gene);
			}
			
			

		}
		
		BufferedWriter writer = new BufferedWriter(new FileWriter("D:\\UMCG\\Genetica\\Projects\\rp3\\interactions\\qtlInteractionsZ8.txt"));
		
		writer.write("Variant\tGene\tPositiveCount\tNegativeCount\tPositive\tNegative\n");
		
		for(Map.Entry<VariantGene, Interaction> interactionEntry : interactions.entrySet()){
			
			VariantGene variantGene = interactionEntry.getKey();
			Interaction interaction = interactionEntry.getValue();
			
			writer.write(variantGene.getVariant());
			writer.write('\t');
			writer.write(variantGene.getGene());
			writer.write('\t');
			writer.write(String.valueOf(interaction.getPositive().size()));
			writer.write('\t');
			writer.write(String.valueOf(interaction.getNegative().size()));
			writer.write('\t');
			
			
			boolean first = true;
			for(String covariate : interaction.getPositive()){
				if(first){
					first = false;
				} else {
					writer.write(';');
				}
				writer.write(covariate);
			}
			writer.write('\t');
			
			first = true;
			for(String covariate : interaction.getNegative()){
				if(first){
					first = false;
				} else {
					writer.write(';');
				}
				writer.write(covariate);
			}
			writer.write('\n');
			
		}
		writer.close();
		
		
		BufferedWriter writer2 = new BufferedWriter(new FileWriter("D:\\UMCG\\Genetica\\Projects\\rp3\\interactions\\covariateInteractionsZ8.txt"));
		
		writer2.write("Covariate\tPositiveCount\tNegativeCount\tPositive\tNegative\n");
		
		for(Map.Entry<String, Interaction> interactionEntry : covariateInteractions.entrySet()){
			
			String covariate = interactionEntry.getKey();
			Interaction interaction = interactionEntry.getValue();
			
			writer2.write(covariate);
			writer2.write('\t');
			writer2.write(String.valueOf(interaction.getPositive().size()));
			writer2.write('\t');
			writer2.write(String.valueOf(interaction.getNegative().size()));
			writer2.write('\t');
			
			
			boolean first = true;
			for(String gene : interaction.getPositive()){
				if(first){
					first = false;
				} else {
					writer2.write(';');
				}
				writer2.write(gene);
			}
			writer2.write('\t');
			
			first = true;
			for(String gene : interaction.getNegative()){
				if(first){
					first = false;
				} else {
					writer2.write(';');
				}
				writer2.write(gene);
			}
			writer2.write('\n');
			
		}
		writer2.close();
		
    }
	
	private static class Interaction{
		
		private final ArrayList<String> positive;
		private final ArrayList<String> negative;

		public Interaction() {
			this.positive = new ArrayList<String>();
			this.negative = new ArrayList<String>();
		}

		public ArrayList<String> getPositive() {
			return positive;
		}

		public ArrayList<String> getNegative() {
			return negative;
		}
		
		public void addPositive(String covariate){
			positive.add(covariate);
		}
		
		public void addNegative(String covariate){
			negative.add(covariate);
		}
		
		
	}
	
	private static class VariantGene{
		
		private final String variant;
		private final String gene;

		public VariantGene(String variant, String gene) {
			this.variant = variant;
			this.gene = gene;
		}

		public String getVariant() {
			return variant;
		}

		public String getGene() {
			return gene;
		}

		@Override
		public int hashCode() {
			return this.variant.hashCode();
		}

		@Override
		public boolean equals(Object obj) {
			if (obj == null) {
				return false;
			}
			if (getClass() != obj.getClass()) {
				return false;
			}
			final VariantGene other = (VariantGene) obj;
			if ((this.variant == null) ? (other.variant != null) : !this.variant.equals(other.variant)) {
				return false;
			}
			if ((this.gene == null) ? (other.gene != null) : !this.gene.equals(other.gene)) {
				return false;
			}
			return true;
		}
		
		
			
	}

}
