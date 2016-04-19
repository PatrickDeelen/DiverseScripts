/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.umcg.collapsegenotypedsamples;

import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author patri
 */
public class Collapse {
	
	private final HashSet<String> samples = new HashSet<>();

	public Collapse(String sample) {
		samples.add(sample);
	}
	
	public Set<String> getSamples() {
		return samples;
	}
	
	public boolean addSample(String sample){
		return samples.add(sample);
	}
	
	public void addSamples(Set<String> samplesToAdd){
		this.samples.addAll(samplesToAdd);
	}
	
}
