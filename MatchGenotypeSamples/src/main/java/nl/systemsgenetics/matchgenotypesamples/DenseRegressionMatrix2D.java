package nl.systemsgenetics.matchgenotypesamples;

import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 *
 * @author Patrick Deelen
 */
    /**
     * Constructs a view with the given parameters.
     * 
     * rows
     *            the number of rows the matrix shall have.
     * columns
     *            the number of columns the matrix shall have.
     * elements
     *            the cells.
     * rowZero
     *            the position of the first element.
     * columnZero
     *            the position of the first element.
     * rowStride
     *            the number of elements between two rows, i.e.
     *            <tt>index(i+1,j)-index(i,j)</tt>.
     * columnStride
     *            the number of elements between two columns, i.e.
     *            <tt>index(i,j+1)-index(i,j)</tt>.
     * 
     * */

public class DenseRegressionMatrix2D {

    private final SimpleRegression[] elements;
    private final int cols;
    
    public DenseRegressionMatrix2D(int rows, int columns) {
        
        cols = columns;
        
        elements = new SimpleRegression[rows * columns];
        
    }
    
    public void setQuick(int row, int column, SimpleRegression value) {
        
        
        elements[row * cols + column] = value;
    }

    public SimpleRegression getQuick(int row, int column) {
        
        return elements[row * cols + column];
        
    }

}
