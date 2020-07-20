/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.transeqtlenrichment;

import java.io.File;

/**
 *
 * @author patri
 */
public class TransEqtlExplanationSource {
	
	final String name;
	final File path;

	public TransEqtlExplanationSource(String name, File path) {
		this.name = name;
		this.path = path;
	}

	public String getName() {
		return name;
	}

	public File getPath() {
		return path;
	}
	
}
