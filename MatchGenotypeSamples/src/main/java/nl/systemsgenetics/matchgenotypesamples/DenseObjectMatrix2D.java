package nl.systemsgenetics.matchgenotypesamples;

/**
 *
 * @author Patrick Deelen
 */
public class DenseObjectMatrix2D<E> {

	private final cern.colt.matrix.tobject.impl.DenseObjectMatrix2D matrix;

	public DenseObjectMatrix2D(int rows, int columns) {
		matrix = new cern.colt.matrix.tobject.impl.DenseObjectMatrix2D(rows, columns);
	}
	
	public E getQuick(int row, int column){
		return (E) matrix.getQuick(row, column);
	}
	
	public void setQuick(int row, int column, E value){
		matrix.setQuick(row, column, value);
	}
	
}
