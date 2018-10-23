/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.smallscripts;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.TTest;

/**
 *
 * @author patri
 */
public class Test {

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) {
		MannWhitneyUTest x = new MannWhitneyUTest();
		
		
		double[] a = {10,5,3};
		double[] b = {2,4,6};
		
		System.out.println(x.mannWhitneyU(a, b));
		System.out.println(x.mannWhitneyUTest(a, b));
		
		System.out.println(x.mannWhitneyU(b, a));
		System.out.println(x.mannWhitneyUTest(b, a));
		
		
		TTest y = new TTest();
		//y.tTest(a, a)
		
		
	}
	
}
