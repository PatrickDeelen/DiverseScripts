/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.transeqtlenrichment;

/**
 *
 * @author patri
 */
public class TransEqtlExplanation {
	
	final String cisGene;
	final String explanation;

	public TransEqtlExplanation(String cisGene, String explanation) {
		this.cisGene = cisGene;
		this.explanation = explanation;
	}

	public String getCisGene() {
		return cisGene;
	}

	public String getExplanation() {
		return explanation;
	}
	
}
