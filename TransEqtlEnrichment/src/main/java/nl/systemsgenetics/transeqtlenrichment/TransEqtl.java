/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.transeqtlenrichment;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author patri
 */
public class TransEqtl {
	
	final String snp;
	final String transGene;
	final HashMap<String, ArrayList<TransEqtlExplanation>> transEqtlExplanations = new HashMap<>();

	public TransEqtl(String snp, String transGene) {
		this.snp = snp;
		this.transGene = transGene;
	}

	public String getSnp() {
		return snp;
	}

	public String getTransGene() {
		return transGene;
	}

	public void addExplanation(TransEqtlExplanation transEqtlExplanation){
		
		ArrayList<TransEqtlExplanation> transEqtlExplanationsList = transEqtlExplanations.get(transEqtlExplanation.getExplanation());
		if(transEqtlExplanationsList == null){
			transEqtlExplanationsList = new ArrayList<>();
			transEqtlExplanations.put(transEqtlExplanation.getExplanation(), transEqtlExplanationsList);
		}
		transEqtlExplanationsList.add(transEqtlExplanation);
	}

	public HashMap<String, ArrayList<TransEqtlExplanation>> getTransEqtlExplanations() {
		return transEqtlExplanations;
	}
	
}
