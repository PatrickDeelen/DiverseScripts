/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.smallscripts;

import java.io.IOException;
import org.ujmp.core.DenseMatrix;
import org.ujmp.core.Matrix;
import org.ujmp.core.doublematrix.impl.DenseFileMatrix;



/**
 *
 * @author patri
 */
public class TestUjmp {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {
		
		Matrix dense = DenseMatrix.Factory.linkToBinaryFile("C:\\UMCG\\Genetica\\test.bin", 11,11);
		
		//dense.setAsDouble(10, 10,10);
		
		
		System.out.println(dense.getAsDouble(10, 10));

		
		
	}
	
}
