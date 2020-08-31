// --- BEGIN LICENSE BLOCK ---
/* 
 * Copyright (c) 2009, Mikio L. Braun
 * Copyright (c) 2008, Johannes Schaback
 * Copyright (c) 2009, Jan Saputra Müller
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of the Technische Universität Berlin nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
// --- END LICENSE BLOCK ---
package org.jblas;

import units.qual.Dimensionless;
import org.jblas.exceptions.SizeException;
import org.jblas.ranges.Range;
import org.jblas.util.Random;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.StringWriter;
import java.util.AbstractList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * A general matrix class for <tt>double</tt> typed values.
 * 
 * Don't be intimidated by the large number of methods this function defines. Most
 * are overloads provided for ease of use. For example, for each arithmetic operation,
 * up to six overloaded versions exist to handle in-place computations, and
 * scalar arguments (like adding a number to all elements of a matrix).
 * 
 * <h3>Construction</h3>
 * 
 * <p>To construct a two-dimensional matrices, you can use the following constructors
 * and static methods.</p>
 * 
 * <table class="my">
 * <tr><th>Method<th>Description
 * <tr><td>DoubleMatrix(m,n, [value1, value2, value3...])<td>Values are filled in column by column.
 * <tr><td>DoubleMatrix(new double[][] {{value1, value2, ...}, ...}<td>Inner arrays are rows.
 * <tr><td>DoubleMatrix.zeros(m,n) <td>Initial values set to 0.0.
 * <tr><td>DoubleMatrix.ones(m,n) <td>Initial values set to 1.0.
 * <tr><td>DoubleMatrix.rand(m,n) <td>Values drawn at random between 0.0 and 1.0.
 * <tr><td>DoubleMatrix.randn(m,n) <td>Values drawn from normal distribution.
 * <tr><td>DoubleMatrix.eye(n) <td>Unit matrix (values 0.0 except for 1.0 on the diagonal).
 * <tr><td>DoubleMatrix.diag(array) <td>Diagonal matrix with given diagonal elements.
 * </table>
 * 
 * <p>Alternatively, you can construct (column) vectors, if you just supply the length
 * using the following constructors and static methods.</p>
 * 
 * <table class="my">
 * <tr><th>Method</th>                        <th>Description</th></tr>
 * <tr><td>DoubleMatrix(m)</td>               <td>Constructs a column vector.</td></tr>
 * <tr><td>DoubleMatrix(new double[] {value1, value2, ...})</td><td>Constructs a column vector.</td></tr>
 * <tr><td>DoubleMatrix.zeros(m)</td>         <td>Initial values set to 0.0.</td></tr>
 * <tr><td>DoubleMatrix.ones(m)</td>          <td>Initial values set to 1.0.</td></tr>
 * <tr><td>DoubleMatrix.rand(m)</td>          <td>Values drawn at random between 0.0 and 1.0.</td></tr>
 * <tr><td>DoubleMatrix.randn(m)</td>         <td>Values drawn from normal distribution.</td></tr>
 * <tr><td>DoubleMatrix.linspace(a, b, n)</td><td>n linearly spaced values from a to b.</td></tr>
 * <tr><td>DoubleMatrix.logspace(a, b, n)</td><td>n logarithmically spaced values form 10^a to 10^b.</td></tr>
 * </table>
 * 
 * <p>You can also construct new matrices by concatenating matrices either horziontally
 * or vertically:</p>
 * 
 * <table class="my">
 * <tr><th>Method<th>Description
 * <tr><td>x.concatHorizontally(y)<td>New matrix will be x next to y.
 * <tr><td>x.concatVertically(y)<td>New matrix will be x atop y.
 * </table>
 * 
 * <h3>Element Access, Copying and Duplication</h3>
 * 
 * <p>To access individual elements, or whole rows and columns, use the following
 * methods:<p>
 * 
 * <table class="my">
 * <tr><th>x.Method<th>Description
 * <tr><td>x.get(i,j)<td>Get element in row i and column j.
 * <tr><td>x.put(i, j, v)<td>Set element in row i and column j to value v
 * <tr><td>x.get(i)<td>Get the ith element of the matrix (traversing rows first).
 * <tr><td>x.put(i, v)<td>Set the ith element of the matrix (traversing rows first).
 * <tr><td>x.getColumn(i)<td>Get a copy of column i.
 * <tr><td>x.putColumn(i, c)<td>Put matrix c into column i.
 * <tr><td>x.getRow(i)<td>Get a copy of row i.
 * <tr><td>x.putRow(i, c)<td>Put matrix c into row i.
 * <tr><td>x.swapColumns(i, j)<td>Swap the contents of columns i and j.
 * <tr><td>x.swapRows(i, j)<td>Swap the contents of rows i and j.
 * </table>
 * 
 * <p>For <tt>get</tt> and <tt>put</tt>, you can also pass integer arrays,
 * DoubleMatrix objects, or Range objects, which then specify the indices used 
 * as follows:
 * 
 * <ul>
 * <li><em>integer array:</em> the elements will be used as indices.
 * <li><em>DoubleMatrix object:</em> non-zero entries specify the indices.
 * <li><em>Range object:</em> see below.
 * </ul>
 * 
 * <p>When using <tt>put</tt> with multiple indices, the assigned object must
 * have the correct size or be a scalar.</p>
 *
 * <p>There exist the following Range objects. The Class <tt>RangeUtils</tt> also
 * contains the a number of handy helper methods for constructing these ranges.</p>
 * <table class="my">
 * <tr><th>Class <th>RangeUtils method <th>Indices
 * <tr><td>AllRange <td>all() <td>All legal indices.
 * <tr><td>PointRange <td>point(i) <td> A single point.
 * <tr><td>IntervalRange <td>interval(a, b)<td> All indices from a to b (inclusive)
 * <tr><td rowspan=3>IndicesRange <td>indices(int[])<td> The specified indices.
 * <tr><td>indices(DoubleMatrix)<td>The specified indices.
 * <tr><td>find(DoubleMatrix)<td>The non-zero entries of the matrix.
 * </table>
 * 
 * <p>The following methods can be used for duplicating and copying matrices.</p>
 * 
 * <table class="my">
 * <tr><th>Method<th>Description
 * <tr><td>x.dup()<td>Get a copy of x.
 * <tr><td>x.copy(y)<td>Copy the contents of y to x (possible resizing x).
 * </table>
 *    
 * <h3>Size and Shape</h3>
 * 
 * <p>The following methods permit to access the size of a matrix and change its size or shape.</p>
 * 
 * <table class="my">
 * <tr><th>x.Method<th>Description
 * <tr><td>x.rows<td>Number of rows.
 * <tr><td>x.columns<td>Number of columns.
 * <tr><td>x.length<td>Total number of elements.
 * <tr><td>x.isEmpty()<td>Checks whether rows == 0 and columns == 0.
 * <tr><td>x.isRowVector()<td>Checks whether rows == 1.
 * <tr><td>x.isColumnVector()<td>Checks whether columns == 1.
 * <tr><td>x.isVector()<td>Checks whether rows == 1 or columns == 1.
 * <tr><td>x.isSquare()<td>Checks whether rows == columns.
 * <tr><td>x.isScalar()<td>Checks whether length == 1.
 * <tr><td>x.resize(r, c)<td>Resize the matrix to r rows and c columns, discarding the content.
 * <tr><td>x.reshape(r, c)<td>Resize the matrix to r rows and c columns.<br> Number of elements must not change.
 * </table>
 * 
 * <p>The size is stored in the <tt>rows</tt> and <tt>columns</tt> member variables.
 * The total number of elements is stored in <tt>length</tt>. Do not change these
 * values unless you know what you're doing!</p>
 * 
 * <h3>Arithmetics</h3>
 * 
 * <p>The usual arithmetic operations are implemented. Each operation exists in a
 * in-place version, recognizable by the suffix <tt>"i"</tt>, to which you can supply
 * the result matrix (or <tt>this</tt> is used, if missing). Using in-place operations
 * can also lead to a smaller memory footprint, as the number of temporary objects is
 * reduced (although the JVM garbage collector is usually pretty good at reusing these
 * temporary object immediately with little overhead.)</p>
 * 
 * <p>Whenever you specify a result vector, the result vector must already have the
 * correct dimensions.</p>
 * 
 * <p>For example, you can add two matrices using the <tt>add</tt> method. If you want
 * to store the result in of <tt>x + y</tt> in <tt>z</tt>, type
 * <span class=code>
 * x.addi(y, z)   // computes x = y + z.
 * </span>
 * Even in-place methods return the result, such that you can easily chain in-place methods,
 * for example:
 * <span class=code>
 * x.addi(y).addi(z) // computes x += y; x += z
 * </span></p> 
 *
 * <p>Methods which operate element-wise only make sure that the length of the matrices
 * is correct. Therefore, you can add a 3 * 3 matrix to a 1 * 9 matrix, for example.</p>
 * 
 * <p>Finally, there exist versions which take doubles instead of DoubleMatrix Objects
 * as arguments. These then compute the operation with the same value as the
 * right-hand-side. The same effect can be achieved by passing a DoubleMatrix with
 * exactly one element.</p>
 * 
 * <table class="my">
 * <tr><th>Operation <th>Method <th>Comment
 * <tr><td>x + y <td>x.add(y)			<td>
 * <tr><td>x - y <td>x.sub(y), y.rsub(x) <td>rsub subtracts left from right hand side
 * <tr><td rowspan=3>x * y 	<td>x.mul(y) <td>element-wise multiplication 
 * <tr>                     <td>x.mmul(y)<td>matrix-matrix multiplication
 * <tr>                     <td>x.dot(y) <td>scalar-product
 * <tr><td>x / y <td>x.div(y), y.rdiv(x) <td>rdiv divides right hand side by left hand side.
 * <tr><td>- x	 <td>x.neg()				<td>
 * </table>
 * 
 * <p>There also exist operations which work on whole columns or rows.</p>
 * 
 * <table class="my">
 * <tr><th>Method</th>           <th>Description</th></tr>
 * <tr><td>x.addRowVector</td>   <td>adds a vector to each row (addiRowVector works in-place)</td></tr>
 * <tr><td>x.addColumnVector</td><td>adds a vector to each column</td></tr>
 * <tr><td>x.subRowVector</td>   <td>subtracts a vector from each row</td></tr>
 * <tr><td>x.subColumnVector</td><td>subtracts a vector from each column</td></tr>
 * <tr><td>x.mulRowVector</td>   <td>Multiplies each row by a vector (elementwise)</td></tr>
 * <tr><td>x.mulColumnVector</td><td>Multiplies each column by a vector (elementwise)</td></tr>
 * <tr><td>x.divRowVector</td>   <td>Divide each row by a vector (elementwise)</td></tr>
 * <tr><td>x.divColumnVector</td><td>Divide each column by a vector (elementwise)</td></tr>
 * <tr><td>x.mulRow</td>         <td>Multiplies a row by a scalar</td></tr>
 * <tr><td>x.mulColumn</td>      <td>Multiplies a column by a scalar</td></tr>
 * </table>
 * 
 * <p>In principle, you could achieve the same result by first calling getColumn(), 
 * adding, and then calling putColumn, but these methods are much faster.</p>
 * 
 * <p>The following comparison operations are available</p>
 *  
 * <table class="my">
 * <tr><th>Operation <th>Method
 * <tr><td>x &lt; y		<td>x.lt(y)
 * <tr><td>x &lt;= y	<td>x.le(y)
 * <tr><td>x &gt; y		<td>x.gt(y)
 * <tr><td>x &gt;= y	<td>x.ge(y)
 * <tr><td>x == y		<td>x.eq(y)
 * <tr><td>x != y		<td>x.ne(y)
 * </table>
 *
 * <p> Logical operations are also supported. For these operations, a value different from
 * zero is treated as "true" and zero is treated as "false". All operations are carried
 * out elementwise.</p>
 * 
 * <table class="my">
 * <tr><th>Operation <th>Method
 * <tr><td>x & y 	<td>x.and(y)
 * <tr><td>x | y 	<td>x.or(y)
 * <tr><td>x ^ y	<td>x.xor(y)
 * <tr><td>! x		<td>x.not()
 * </table>
 * 
 * <p>Finally, there are a few more methods to compute various things:</p>
 * 
 * <table class="my">
 * <tr><th>Method <th>Description
 * <tr><td>x.max() <td>Return maximal element
 * <tr><td>x.argmax() <td>Return index of largest element
 * <tr><td>x.min() <td>Return minimal element
 * <tr><td>x.argmin() <td>Return index of largest element
 * <tr><td>x.columnMins() <td>Return column-wise minima
 * <tr><td>x.columnArgmins() <td>Return column-wise index of minima
 * <tr><td>x.columnMaxs() <td>Return column-wise maxima
 * <tr><td>x.columnArgmaxs() <td>Return column-wise index of maxima
 * </table>
 * 
 * @author Mikio Braun, Johannes Schaback
 */
public class DoubleMatrix implements Serializable {

    /** Number of rows. */
    public int rows;
    /** Number of columns. */
    public int columns;
    /** Total number of elements (for convenience). */
    public @Dimensionless int length;
    /** The actual data stored by rows (that is, row 0, row 1...). */
    public @Dimensionless double @Dimensionless [] data = null; // rows are contiguous
    public static final @Dimensionless DoubleMatrix EMPTY = new @Dimensionless DoubleMatrix();
    static final @Dimensionless long serialVersionUID = ((@Dimensionless long) (-1249281332731183060L));

    // Precompile regex patterns
    private static final @Dimensionless Pattern SEMICOLON = Pattern.compile(";");
    private static final @Dimensionless Pattern WHITESPACES = Pattern.compile("\\s+");
    private static final @Dimensionless Pattern COMMA = Pattern.compile(",");

    /**************************************************************************
     *
     * Constructors and factory functions
     *
     **************************************************************************/
    /** Create a new matrix with <i>newRows</i> rows, <i>newColumns</i> columns
     * using <i>newData></i> as the data. Note that any change to the DoubleMatrix
     * will change the input array, too.
     *
     * @param newRows the number of rows of the new matrix
     * @param newColumns the number of columns of the new matrix
     * @param newData the data array to be used. Data must be stored by column (column-major)
     */
    public DoubleMatrix(@Dimensionless int newRows, @Dimensionless int newColumns, @Dimensionless double @Dimensionless ... newData) {
        rows = newRows;
        columns = newColumns;
        length = rows * columns;

        if (newData != null && newData.length != newRows * newColumns) {
            throw new @Dimensionless IllegalArgumentException(
                    "Passed data must match matrix dimensions.");
        }

        data = newData;
        //System.err.printf("%d * %d matrix created\n", rows, columns);
    }

    /**
     * Creates a new <i>n</i> times <i>m</i> <tt>DoubleMatrix</tt>.
     *
     * @param newRows the number of rows (<i>n</i>) of the new matrix.
     * @param newColumns the number of columns (<i>m</i>) of the new matrix.
     */
    public DoubleMatrix(@Dimensionless int newRows, @Dimensionless int newColumns) {
        this(newRows, newColumns, new @Dimensionless double @Dimensionless [newRows * newColumns]);
    }

    /**
     * Creates a new <tt>DoubleMatrix</tt> of size 0 times 0.
     */
    public DoubleMatrix() {
        this(((@Dimensionless int) (0)), ((@Dimensionless int) (0)), (@Dimensionless double @Dimensionless []) null);
    }

    /**
     * Create a Matrix of length <tt>len</tt>. This creates a column vector.
     *
     * @param len
     */
    public DoubleMatrix(int len) {
        this(len, ((@Dimensionless int) (1)), new @Dimensionless double @Dimensionless [len]);
    }
    
    /**
     * Create a a column vector using <i>newData</i> as the data array.
     * Note that any change to the created DoubleMatrix will change in input array.
     */
    public DoubleMatrix(@Dimensionless double @Dimensionless [] newData) {
        this(newData.length, ((@Dimensionless int) (1)), newData);
    }

    /**
     * Creates a new matrix by reading it from a file.
     *
     * @param filename the path and name of the file to read the matrix from
     * @throws IOException
     */
    public DoubleMatrix(@Dimensionless String filename) throws IOException {
        load(filename);
    }

    /**
     * Creates a new <i>n</i> times <i>m</i> <tt>DoubleMatrix</tt> from
     * the given <i>n</i> times <i>m</i> 2D data array. Note that the input array
     * is copied and any change to the DoubleMatrix will not change the input array.
     * The first dimension of the array makes the
     * rows (<i>n</i>) and the second dimension the columns (<i>m</i>). For example, the
     * given code <br/><br/>
     * <code>new DoubleMatrix(new double[][]{{1d, 2d, 3d}, {4d, 5d, 6d}, {7d, 8d, 9d}}).print();</code><br/><br/>
     * will constructs the following matrix:
     * <pre>
     * 1.0	2.0	3.0
     * 4.0	5.0	6.0
     * 7.0	8.0	9.0
     * </pre>.
     * @param data <i>n</i> times <i>m</i> data array
     */
    public DoubleMatrix(@Dimensionless double @Dimensionless [] @Dimensionless [] data) {
        this(data.length, data[((@Dimensionless int) (0))].length);

        for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
            assert (data[r].length == columns);
        }

        for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
                put(r, c, data[r][c]);
            }
        }
    }

    /**
     * Creates a DoubleMatrix column vector from the given List&lt;Double&rt;.
     *
     * @param data data from which the entries are taken.
     */
    public DoubleMatrix(@Dimensionless List<@Dimensionless Double> data) {
        this(data.size());

        @Dimensionless
        int c = ((@Dimensionless int) (0));
        for (@Dimensionless double d : data) {
            put(c++, d);
        }
    }

    /**
     * Construct DoubleMatrix from ASCII representation.
     *
     * This is not very fast, but can be quiet useful when
     * you want to "just" construct a matrix, for example
     * when testing.
     *
     * The format is semicolon separated rows of space separated values,
     * for example "1 2 3; 4 5 6; 7 8 9".
     */
    public static @Dimensionless DoubleMatrix valueOf(@Dimensionless String text) {
        @Dimensionless
        String @Dimensionless [] rowValues = SEMICOLON.split(text);

        @Dimensionless
        DoubleMatrix result = null;

        // process rest
        for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rowValues.length; r++) {
          @Dimensionless
          String @Dimensionless [] columnValues = WHITESPACES.split(rowValues[r].trim());

            if (r == ((@Dimensionless int) (0))) {
                result = new @Dimensionless DoubleMatrix(rowValues.length, columnValues.length);
            }

            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columnValues.length; c++) {
                result.put(r, c, Double.valueOf(columnValues[c]));
            }
        }

        return result;
    }

    /**
     * Serialization
     */
    private void writeObject(@Dimensionless DoubleMatrix this, @Dimensionless ObjectOutputStream s) throws IOException {
        s.defaultWriteObject();
    }

    private void readObject(@Dimensionless DoubleMatrix this, @Dimensionless ObjectInputStream s) throws IOException, ClassNotFoundException {
        s.defaultReadObject();
    }

    /** Create matrix with random values uniformly in 0..1. */
    public static @Dimensionless DoubleMatrix rand(@Dimensionless int rows, @Dimensionless int columns) {
        @Dimensionless
        DoubleMatrix m = new @Dimensionless DoubleMatrix(rows, columns);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows * columns; i++) {
            m.data[i] = (@Dimensionless double) Random.nextDouble();
        }

        return m;
    }

    /** Creates a column vector with random values uniformly in 0..1. */
    public static @Dimensionless DoubleMatrix rand(@Dimensionless int len) {
        return rand(len, ((@Dimensionless int) (1)));
    }

    /** Create matrix with normally distributed random values. */
    public static @Dimensionless DoubleMatrix randn(@Dimensionless int rows, @Dimensionless int columns) {
        @Dimensionless
        DoubleMatrix m = new @Dimensionless DoubleMatrix(rows, columns);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows * columns; i++) {
            m.data[i] = (@Dimensionless double) Random.nextGaussian();
        }

        return m;
    }

    /** Create column vector with normally distributed random values. */
    public static @Dimensionless DoubleMatrix randn(@Dimensionless int len) {
        return randn(len, ((@Dimensionless int) (1)));
    }

    /** Creates a new matrix in which all values are equal 0. */
    public static @Dimensionless DoubleMatrix zeros(@Dimensionless int rows, @Dimensionless int columns) {
        return new @Dimensionless DoubleMatrix(rows, columns);
    }

    /** Creates a column vector of given length. */
    public static @Dimensionless DoubleMatrix zeros(@Dimensionless int length) {
        return zeros(length, ((@Dimensionless int) (1)));
    }

    /** Creates a new matrix in which all values are equal 1. */
    public static @Dimensionless DoubleMatrix ones(@Dimensionless int rows, @Dimensionless int columns) {
        @Dimensionless
        DoubleMatrix m = new @Dimensionless DoubleMatrix(rows, columns);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows * columns; i++) {
            m.put(i, ((@Dimensionless double) (1.0)));
        }

        return m;
    }

    /** Creates a column vector with all elements equal to 1. */
    public static @Dimensionless DoubleMatrix ones(@Dimensionless int length) {
        return ones(length, ((@Dimensionless int) (1)));
    }

    /** Construct a new n-by-n identity matrix. */
    public static @Dimensionless DoubleMatrix eye(@Dimensionless int n) {
        @Dimensionless
        DoubleMatrix m = new @Dimensionless DoubleMatrix(n, n);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
            m.put(i, i, ((@Dimensionless double) (1.0)));
        }

        return m;
    }

    /**
     * Creates a new matrix where the values of the given vector are the diagonal values of
     * the matrix.
     */
    public static @Dimensionless DoubleMatrix diag(@Dimensionless DoubleMatrix x) {
        @Dimensionless
        DoubleMatrix m = new @Dimensionless DoubleMatrix(x.length, x.length);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++) {
            m.put(i, i, x.get(i));
        }

        return m;
    }

  /**
   * Construct a matrix of arbitrary shape and set the diagonal according
   * to a passed vector.
   *
   * length of needs to be smaller than rows or columns.
   *
   * @param x vector to fill the diagonal with
   * @param rows number of rows of the resulting matrix
   * @param columns number of columns of the resulting matrix
   * @return a matrix with dimensions rows * columns whose diagonal elements are filled by x
   */
    public static @Dimensionless DoubleMatrix diag(@Dimensionless DoubleMatrix x, @Dimensionless int rows, @Dimensionless int columns) {
      @Dimensionless
      DoubleMatrix m = new @Dimensionless DoubleMatrix(rows, columns);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++) {
          m.put(i, i, x.get(i));
      }

      return m;
    }

    /**
     * Create a 1-by-1 matrix. For many operations, this matrix functions like a
     * normal double.
     */
    public static @Dimensionless DoubleMatrix scalar(@Dimensionless double s) {
        @Dimensionless
        DoubleMatrix m = new @Dimensionless DoubleMatrix(((@Dimensionless int) (1)), ((@Dimensionless int) (1)));
        m.put(((@Dimensionless int) (0)), ((@Dimensionless int) (0)), s);
        return m;
    }

    /** Test whether a matrix is scalar. */
    public @Dimensionless boolean isScalar(@Dimensionless DoubleMatrix this) {
        return length == ((@Dimensionless int) (1));
    }

    /** Return the first element of the matrix. */
    public @Dimensionless double scalar(@Dimensionless DoubleMatrix this) {
        return get(((@Dimensionless int) (0)));
    }

    /**
     * Construct a column vector whose entries are logarithmically spaced points from
     * 10^lower to 10^upper using the specified number of steps
     *
     * @param lower starting exponent
     * @param upper ending exponent
     * @param size number of steps
     * @return a column vector with (10^lower, ... 10^upper) with size many entries.
     */
    public static @Dimensionless DoubleMatrix logspace(@Dimensionless double lower, @Dimensionless double upper, @Dimensionless int size) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(size);
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < size; i++) {
            @Dimensionless
            double t = (@Dimensionless double) i / (size - ((@Dimensionless int) (1)));
            @Dimensionless
            double e = lower * (((@Dimensionless int) (1)) - t) + t * upper;
            result.put(i, (@Dimensionless double) Math.pow(((@Dimensionless double) (10.0)), e));
        }
        return result;
    }

    /**
     * Construct a column vector whose entries are linearly spaced points from lower to upper with size
     * many steps.
     *
     * @param lower starting value
     * @param upper end value
     * @param size number of steps
     * @return a column vector of size (lower, ..., upper) with size many entries.
     */
    public static @Dimensionless DoubleMatrix linspace(@Dimensionless int lower, @Dimensionless int upper, @Dimensionless int size) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(size);
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < size; i++) {
            @Dimensionless
            double t = (@Dimensionless double) i / (size - ((@Dimensionless int) (1)));
            result.put(i, lower * (((@Dimensionless int) (1)) - t) + t * upper);
        }
        return result;
    }

    /**
     * Concatenates two matrices horizontally. Matrices must have identical
     * numbers of rows.
     */
    public static @Dimensionless DoubleMatrix concatHorizontally(@Dimensionless DoubleMatrix A, @Dimensionless DoubleMatrix B) {
        if (A.rows != B.rows) {
            throw new @Dimensionless SizeException("Matrices don't have same number of rows.");
        }

        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(A.rows, A.columns + B.columns);
        SimpleBlas.copy(A, result);
        JavaBlas.rcopy(B.length, B.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), result.data, A.length, ((@Dimensionless int) (1)));
        return result;
    }

    /**
     * Concatenates two matrices vertically. Matrices must have identical
     * numbers of columns.
     */
    public static @Dimensionless DoubleMatrix concatVertically(@Dimensionless DoubleMatrix A, @Dimensionless DoubleMatrix B) {
        if (A.columns != B.columns) {
            throw new @Dimensionless SizeException("Matrices don't have same number of columns (" + A.columns + " != " + B.columns + ".");
        }

        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(A.rows + B.rows, A.columns);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < A.columns; i++) {
            JavaBlas.rcopy(A.rows, A.data, A.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)), result.data, result.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)));
            JavaBlas.rcopy(B.rows, B.data, B.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)), result.data, result.index(A.rows, i), ((@Dimensionless int) (1)));
        }

        return result;
    }

    /**************************************************************************
     * Working with slices (Man! 30+ methods just to make this a bit flexible...)
     */
    /** Get all elements specified by the linear indices. */
    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] indices) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(indices.length);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++) {
            result.put(i, get(indices[i]));
        }

        return result;
    }

    /** Get all elements for a given row and the specified columns. */
    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(((@Dimensionless int) (1)), indices.length);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++) {
            result.put(i, get(r, indices[i]));
        }

        return result;
    }

    /** Get all elements for a given column and the specified rows. */
    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(indices.length, ((@Dimensionless int) (1)));

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++) {
            result.put(i, get(indices[i], c));
        }

        return result;
    }

    /** Get all elements from the specified rows and columns. */
    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(rindices.length, cindices.length);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++) {
            for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++) {
                result.put(i, j, get(rindices[i], cindices[j]));
            }
        }

        return result;
    }

    /** Get elements from specified rows and columns. */
    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless Range rs, @Dimensionless Range cs) {
        rs.init(((@Dimensionless int) (0)), rows);
        cs.init(((@Dimensionless int) (0)), columns);
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(rs.length(), cs.length());

        for (; rs.hasMore(); rs.next()) {
            cs.init(((@Dimensionless int) (0)), columns);
            for (; cs.hasMore(); cs.next()) {
                result.put(rs.index(), cs.index(), get(rs.value(), cs.value()));
            }
        }

        return result;
    }

    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless Range rs, @Dimensionless int c) {
        rs.init(((@Dimensionless int) (0)), rows);
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(rs.length(), ((@Dimensionless int) (1)));

        for (; rs.hasMore(); rs.next()) {
            result.put(rs.index(), ((@Dimensionless int) (0)), get(rs.value(), c));
        }

        return result;
    }

    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless int r, @Dimensionless Range cs) {
        cs.init(((@Dimensionless int) (0)), columns);
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(((@Dimensionless int) (1)), cs.length());

        for (; cs.hasMore(); cs.next()) {
            result.put(((@Dimensionless int) (0)), cs.index(), get(r, cs.value()));
        }

        return result;

    }

    /** Get elements specified by the non-zero entries of the passed matrix. */
    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix indices) {
        return get(indices.findIndices());
    }

    /**
     * Get elements from a row and columns as specified by the non-zero entries of
     * a matrix.
     */
    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless int r, @Dimensionless DoubleMatrix indices) {
        return get(r, indices.findIndices());
    }

    /**
     * Get elements from a column and rows as specified by the non-zero entries of
     * a matrix.
     */
    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix indices, @Dimensionless int c) {
        return get(indices.findIndices(), c);
    }

    /**
     * Get elements from columns and rows as specified by the non-zero entries of
     * the passed matrices.
     */
    public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix rindices, @Dimensionless DoubleMatrix cindices) {
        return get(rindices.findIndices(), cindices.findIndices());
    }

    /** Return all elements with linear index a, a + 1, ..., b - 1.*/
    public @Dimensionless DoubleMatrix getRange(@Dimensionless DoubleMatrix this, @Dimensionless int a, @Dimensionless int b) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(b - a);

        for (@Dimensionless int k = ((@Dimensionless int) (0)); k < b - a; k++) {
            result.put(k, get(a + k));
        }

        return result;
    }

    /** Get elements from a row and columns <tt>a</tt> to <tt>b</tt>. */
    public @Dimensionless DoubleMatrix getColumnRange(@Dimensionless DoubleMatrix this, @Dimensionless int r, @Dimensionless int a, @Dimensionless int b) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(((@Dimensionless int) (1)), b - a);

        for (@Dimensionless int k = ((@Dimensionless int) (0)); k < b - a; k++) {
            result.put(k, get(r, a + k));
        }

        return result;
    }

    /** Get elements from a column and rows <tt>a/tt> to <tt>b</tt>. */
    public @Dimensionless DoubleMatrix getRowRange(@Dimensionless DoubleMatrix this, @Dimensionless int a, @Dimensionless int b, @Dimensionless int c) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(b - a);

        for (@Dimensionless int k = ((@Dimensionless int) (0)); k < b - a; k++) {
            result.put(k, get(a + k, c));
        }

        return result;
    }

    /**
     * Get elements from rows <tt>ra</tt> to <tt>rb</tt> and
     * columns <tt>ca</tt> to <tt>cb</tt>.
     */
    public @Dimensionless DoubleMatrix getRange(@Dimensionless DoubleMatrix this, @Dimensionless int ra, @Dimensionless int rb, @Dimensionless int ca, @Dimensionless int cb) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(rb - ra, cb - ca);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rb - ra; i++) {
            for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cb - ca; j++) {
                result.put(i, j, get(ra + i, ca + j));
            }
        }

        return result;
    }

    /** Get whole rows from the passed indices. */
    public @Dimensionless DoubleMatrix getRows(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] rindices) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(rindices.length, columns);
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++) {
            JavaBlas.rcopy(columns, data, index(rindices[i], ((@Dimensionless int) (0))), rows, result.data, result.index(i, ((@Dimensionless int) (0))), result.rows);
        }
        return result;
    }

    /** Get whole rows as specified by the non-zero entries of a matrix. */
    public @Dimensionless DoubleMatrix getRows(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix rindices) {
        return getRows(rindices.findIndices());
    }

    public @Dimensionless DoubleMatrix getRows(@Dimensionless DoubleMatrix this, @Dimensionless Range indices, @Dimensionless DoubleMatrix result) {
        indices.init(((@Dimensionless int) (0)), rows);
        if (result.rows < indices.length()) {
            throw new @Dimensionless SizeException("Result matrix does not have enough rows (" + result.rows + " < " + indices.length() + ")");
        }
        result.checkColumns(columns);

      indices.init(((@Dimensionless int) (0)), rows);
      for (@Dimensionless int r = ((@Dimensionless int) (0)); indices.hasMore(); indices.next(), r++) {
            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
                result.put(r, c, get(indices.value(), c));
            }
        }
        return result;
    }

    public @Dimensionless DoubleMatrix getRows(@Dimensionless DoubleMatrix this, @Dimensionless Range indices) {
        indices.init(((@Dimensionless int) (0)), rows);
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(indices.length(), columns);
        return getRows(indices, result);
    }

    /** Get whole columns from the passed indices. */
    public @Dimensionless DoubleMatrix getColumns(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] cindices) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(rows, cindices.length);
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < cindices.length; i++) {
            JavaBlas.rcopy(rows, data, index(((@Dimensionless int) (0)), cindices[i]), ((@Dimensionless int) (1)), result.data, result.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)));
        }
        return result;
    }

    /** Get whole columns as specified by the non-zero entries of a matrix. */
    public @Dimensionless DoubleMatrix getColumns(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix cindices) {
        return getColumns(cindices.findIndices());
    }


    /** Get whole columns as specified by Range. */
    public @Dimensionless DoubleMatrix getColumns(@Dimensionless DoubleMatrix this, @Dimensionless Range indices, @Dimensionless DoubleMatrix result) {
        indices.init(((@Dimensionless int) (0)), columns);
        if (result.columns < indices.length()) {
            throw new @Dimensionless SizeException("Result matrix does not have enough columns (" + result.columns + " < " + indices.length() + ")");
        }
        result.checkRows(rows);

        indices.init(((@Dimensionless int) (0)), columns);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); indices.hasMore(); indices.next(), c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                result.put(r, c, get(r, indices.value()));
            }
        }
        return result;
    }

    public @Dimensionless DoubleMatrix getColumns(@Dimensionless DoubleMatrix this, @Dimensionless Range indices) {
        indices.init(((@Dimensionless int) (0)), columns);
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(rows, indices.length());
        return getColumns(indices, result);
    }

    /**
     * Assert that the matrix has a certain length.
     * @throws SizeException
     */
    public void checkLength(@Dimensionless DoubleMatrix this, @Dimensionless int l) {
        if (length != l) {
            throw new @Dimensionless SizeException("Matrix does not have the necessary length (" + length + " != " + l + ").");
        }
    }

    /**
     * Asserts that the matrix has a certain number of rows.
     * @throws SizeException
     */
    public void checkRows(@Dimensionless DoubleMatrix this, @Dimensionless int r) {
        if (rows != r) {
            throw new @Dimensionless SizeException("Matrix does not have the necessary number of rows (" + rows + " != " + r + ").");
        }
    }

    /**
     * Asserts that the amtrix has a certain number of columns.
     * @throws SizeException
     */
    public void checkColumns(@Dimensionless DoubleMatrix this, @Dimensionless int c) {
        if (columns != c) {
            throw new @Dimensionless SizeException("Matrix does not have the necessary number of columns (" + columns + " != " + c + ").");
        }
    }


    /**
     * Set elements in linear ordering in the specified indices.
     *
     * For example, <code>a.put(new int[]{ 1, 2, 0 }, new DoubleMatrix(3, 1, 2.0, 4.0, 8.0)</code>
     * does <code>a.put(1, 2.0), a.put(2, 4.0), a.put(0, 8.0)</code>.
     */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless DoubleMatrix x) {
        if (x.isScalar()) {
            return put(indices, x.scalar());
        }
        x.checkLength(indices.length);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++) {
            put(indices[i], x.get(i));
        }

        return this;
    }

    /** Set multiple elements in a row. */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless DoubleMatrix x) {
        if (x.isScalar()) {
            return put(r, indices, x.scalar());
        }
        x.checkColumns(indices.length);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++) {
            put(r, indices[i], x.get(i));
        }

        return this;
    }

    /** Set multiple elements in a row. */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless DoubleMatrix x) {
        if (x.isScalar()) {
            return put(indices, c, x.scalar());
        }
        x.checkRows(indices.length);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++) {
            put(indices[i], c, x.get(i));
        }

        return this;
    }

    /** Put a sub-matrix as specified by the indices. */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless DoubleMatrix x) {
        if (x.isScalar()) {
            return put(rindices, cindices, x.scalar());
        }
        x.checkRows(rindices.length);
        x.checkColumns(cindices.length);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++) {
            for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++) {
                put(rindices[i], cindices[j], x.get(i, j));
            }
        }

        return this;
    }

    /** Put a matrix into specified indices. */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless Range rs, @Dimensionless Range cs, @Dimensionless DoubleMatrix x) {
        rs.init(((@Dimensionless int) (0)), rows);
        cs.init(((@Dimensionless int) (0)), columns);

        x.checkRows(rs.length());
        x.checkColumns(cs.length());

        for (; rs.hasMore(); rs.next()) {
            cs.init(((@Dimensionless int) (0)), columns);
            for (; cs.hasMore(); cs.next()) {
                put(rs.value(), cs.value(), x.get(rs.index(), cs.index()));
            }
        }

        return this;
    }

    /** Put a single value into the specified indices (linear adressing). */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless double v) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++) {
            put(indices[i], v);
        }

        return this;
    }

    /** Put a single value into a row and the specified columns. */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless double v) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++) {
            put(r, indices[i], v);
        }

        return this;
    }

    /** Put a single value into the specified rows of a column. */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless double v) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++) {
            put(indices[i], c, v);
        }

        return this;
    }

    /** Put a single value into the specified rows and columns. */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless double v) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++) {
            for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++) {
                put(rindices[i], cindices[j], v);
            }
        }

        return this;
    }

    /**
     * Put a sub-matrix into the indices specified by the non-zero entries
     * of <tt>indices</tt> (linear adressing).
     */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix indices, @Dimensionless DoubleMatrix v) {
        return put(indices.findIndices(), v);
    }

    /** Put a sub-vector into the specified columns (non-zero entries of <tt>indices</tt>) of a row. */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int r, @Dimensionless DoubleMatrix indices, @Dimensionless DoubleMatrix v) {
        return put(r, indices.findIndices(), v);
    }

    /** Put a sub-vector into the specified rows (non-zero entries of <tt>indices</tt>) of a column. */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix indices, @Dimensionless int c, @Dimensionless DoubleMatrix v) {
        return put(indices.findIndices(), c, v);
    }

    /**
     * Put a sub-matrix into the specified rows and columns (non-zero entries of
     * <tt>rindices</tt> and <tt>cindices</tt>.
     */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix rindices, @Dimensionless DoubleMatrix cindices, @Dimensionless DoubleMatrix v) {
        return put(rindices.findIndices(), cindices.findIndices(), v);
    }

    /**
     * Put a single value into the elements specified by the non-zero
     * entries of <tt>indices</tt> (linear adressing).
     */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix indices, @Dimensionless double v) {
        return put(indices.findIndices(), v);
    }

    /**
     * Put a single value into the specified columns (non-zero entries of
     * <tt>indices</tt>) of a row.
     */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int r, @Dimensionless DoubleMatrix indices, @Dimensionless double v) {
        return put(r, indices.findIndices(), v);
    }

    /**
     * Put a single value into the specified rows (non-zero entries of
     * <tt>indices</tt>) of a column.
     */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix indices, @Dimensionless int c, @Dimensionless double v) {
        return put(indices.findIndices(), c, v);
    }

    /**
     * Put a single value in the specified rows and columns (non-zero entries
     * of <tt>rindices</tt> and <tt>cindices</tt>.
     */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix rindices, @Dimensionless DoubleMatrix cindices, @Dimensionless double v) {
        return put(rindices.findIndices(), cindices.findIndices(), v);
    }

    /** Find the linear indices of all non-zero elements. */
    public @Dimensionless int @Dimensionless [] findIndices(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int len = ((@Dimensionless int) (0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            if (get(i) != ((@Dimensionless double) (0.0))) {
                len++;
            }
        }

        @Dimensionless
        int @Dimensionless [] indices = new @Dimensionless int @Dimensionless [len];
        @Dimensionless
        int c = ((@Dimensionless int) (0));

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            if (get(i) != ((@Dimensionless double) (0.0))) {
                indices[c++] = i;
            }
        }

        return indices;
    }

    /**************************************************************************
     * Basic operations (copying, resizing, element access)
     */
    /** Return transposed copy of this matrix. */
    public DoubleMatrix transpose(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(columns, rows);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows; i++) {
            for (@Dimensionless int j = ((@Dimensionless int) (0)); j < columns; j++) {
                result.put(j, i, get(i, j));
            }
        }

        return result;
    }

    /**
     * Compare two matrices. Returns true if and only if other is also a
     * DoubleMatrix which has the same size and the maximal absolute
     * difference in matrix elements is smaller than the specified tolerance
     */
    public @Dimensionless boolean compare(@Dimensionless DoubleMatrix this, @Dimensionless Object o, @Dimensionless double tolerance) {
        if (!(o instanceof DoubleMatrix)) {
            return false;
        }

        @Dimensionless
        DoubleMatrix other = (@Dimensionless DoubleMatrix) o;

        if (!sameSize(other)) {
            return false;
        }

        @Dimensionless
        DoubleMatrix diff = MatrixFunctions.absi(sub(other));

        return diff.max() / (rows * columns) < tolerance;
    }

    @Override
    public @Dimensionless boolean equals(@Dimensionless DoubleMatrix this, @Dimensionless Object o) {
        if (!(o instanceof DoubleMatrix)) {
            return false;
        }

        @Dimensionless
        DoubleMatrix other = (@Dimensionless DoubleMatrix) o;

        if (!sameSize(other)) {
            return false;
        } else {
            return Arrays.equals(data, other.data);
        }
    }

    @Override
    public @Dimensionless int hashCode(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int hash = ((@Dimensionless int) (7));
        hash = ((@Dimensionless int) (83)) * hash + this.rows;
        hash = ((@Dimensionless int) (83)) * hash + this.columns;
        hash = ((@Dimensionless int) (83)) * hash + Arrays.hashCode(this.data);
        return hash;
    }
    
    /** Resize the matrix. All elements will be set to zero. */
    public void resize(@Dimensionless int newRows, @Dimensionless int newColumns) {
        rows = newRows;
        columns = newColumns;
        length = newRows * newColumns;
        data = new @Dimensionless double @Dimensionless [rows * columns];
    }

    /** Reshape the matrix. Number of elements must not change. */
    public @Dimensionless DoubleMatrix reshape(@Dimensionless DoubleMatrix this, @Dimensionless int newRows, @Dimensionless int newColumns) {
        if (length != newRows * newColumns) {
            throw new @Dimensionless IllegalArgumentException(
                    "Number of elements must not change.");
        }

        rows = newRows;
        columns = newColumns;

        return this;
    }

    /** Generate a new matrix which has the given number of replications of this. */
    public @Dimensionless DoubleMatrix repmat(@Dimensionless DoubleMatrix this, @Dimensionless int rowMult, @Dimensionless int columnMult) {
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(rows * rowMult, columns * columnMult);

        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columnMult; c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rowMult; r++) {
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows; i++) {
                    for (@Dimensionless int j = ((@Dimensionless int) (0)); j < columns; j++) {
                        result.put(r * rows + i, c * columns + j, get(i, j));
                    }
                }
            }
        }
        return result;
    }

    /** Checks whether two matrices have the same size. */
    public @Dimensionless boolean sameSize(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix a) {
        return rows == a.rows && columns == a.columns;
    }

    /** Throws SizeException unless two matrices have the same size. */
    public void assertSameSize(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix a) {
        if (!sameSize(a)) {
            throw new @Dimensionless SizeException("Matrices must have the same size.");
        }
    }

    /** Checks whether two matrices can be multiplied (that is, number of columns of
     * this must equal number of rows of a. */
    public @Dimensionless boolean multipliesWith(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix a) {
        return columns == a.rows;
    }

    /** Throws SizeException unless matrices can be multiplied with one another. */
    public void assertMultipliesWith(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix a) {
        if (!multipliesWith(a)) {
            throw new @Dimensionless SizeException("Number of columns of left matrix must be equal to number of rows of right matrix.");
        }
    }

    /** Checks whether two matrices have the same length. */
    public @Dimensionless boolean sameLength(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix a) {
        return length == a.length;
    }

    /** Throws SizeException unless matrices have the same length. */
    public void assertSameLength(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix a) {
        if (!sameLength(a)) {
            throw new @Dimensionless SizeException("Matrices must have same length (is: " + length + " and " + a.length + ")");
        }
    }

    /** Copy DoubleMatrix a to this. this a is resized if necessary. */
    public @Dimensionless DoubleMatrix copy(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix a) {
        if (!sameSize(a)) {
            resize(a.rows, a.columns);
        }

        System.arraycopy(a.data, ((@Dimensionless int) (0)), data, ((@Dimensionless int) (0)), length);
        return a;
    }

    /**
     * Returns a duplicate of this matrix. Geometry is the same (including offsets, transpose, etc.),
     * but the buffer is not shared.
     */
    public @Dimensionless DoubleMatrix dup(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        DoubleMatrix out = new @Dimensionless DoubleMatrix(rows, columns);

        JavaBlas.rcopy(length, data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), out.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));

        return out;
    }

    /** Swap two columns of a matrix. */
    public @Dimensionless DoubleMatrix swapColumns(@Dimensionless DoubleMatrix this, @Dimensionless int i, @Dimensionless int j) {
        NativeBlas.dswap(rows, data, index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), j), ((@Dimensionless int) (1)));
        return this;
    }

    /** Swap two rows of a matrix. */
    public @Dimensionless DoubleMatrix swapRows(@Dimensionless DoubleMatrix this, @Dimensionless int i, @Dimensionless int j) {
        NativeBlas.dswap(columns, data, index(i, ((@Dimensionless int) (0))), rows, data, index(j, ((@Dimensionless int) (0))), rows);
        return this;
    }

    /** Set matrix element */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless double value) {
        data[index(rowIndex, columnIndex)] = value;
        return this;
    }

    /** Retrieve matrix element */
    public @Dimensionless double get(@Dimensionless DoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex) {
        return data[index(rowIndex, columnIndex)];
    }

    /** Get index of an element */
    public @Dimensionless int index(@Dimensionless DoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex) {
        return rowIndex + rows * columnIndex;
    }

    /** Compute the row index of a linear index. */
    public @Dimensionless int indexRows(@Dimensionless DoubleMatrix this, @Dimensionless int i) {
		  return i - indexColumns(i) * rows;
    }

    /** Compute the column index of a linear index. */
    public @Dimensionless int indexColumns(@Dimensionless DoubleMatrix this, @Dimensionless int i) {
		  return i / rows;
    }

    /** Get a matrix element (linear indexing). */
    public @Dimensionless double get(@Dimensionless DoubleMatrix this, @Dimensionless int i) {
        return data[i];
    }

    /** Set a matrix element (linear indexing). */
    public @Dimensionless DoubleMatrix put(@Dimensionless DoubleMatrix this, @Dimensionless int i, @Dimensionless double v) {
        data[i] = v;
        return this;
    }

    /** Set all elements to a value. */
    public @Dimensionless DoubleMatrix fill(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            put(i, value);
        }
        return this;
    }

    /** Get number of rows. */
    public @Dimensionless int getRows(@Dimensionless DoubleMatrix this) {
        return rows;
    }

    /** Get number of columns. */
    public @Dimensionless int getColumns(@Dimensionless DoubleMatrix this) {
        return columns;
    }

    /** Get total number of elements. */
    public @Dimensionless int getLength(@Dimensionless DoubleMatrix this) {
        return length;
    }

    /** Checks whether the matrix is empty. */
    public @Dimensionless boolean isEmpty(@Dimensionless DoubleMatrix this) {
        return columns == ((@Dimensionless int) (0)) || rows == ((@Dimensionless int) (0));
    }

    /** Checks whether the matrix is square. */
    public @Dimensionless boolean isSquare(@Dimensionless DoubleMatrix this) {
        return columns == rows;
    }

    /** Throw SizeException unless matrix is square. */
    public void assertSquare(@Dimensionless DoubleMatrix this) {
        if (!isSquare()) {
            throw new @Dimensionless SizeException("Matrix must be square!");
        }
    }

    /** Checks whether the matrix is a vector. */
    public @Dimensionless boolean isVector(@Dimensionless DoubleMatrix this) {
        return columns == ((@Dimensionless int) (1)) || rows == ((@Dimensionless int) (1));
    }

    /** Checks whether the matrix is a row vector. */
    public @Dimensionless boolean isRowVector(@Dimensionless DoubleMatrix this) {
        return rows == ((@Dimensionless int) (1));
    }

    /** Checks whether the matrix is a column vector. */
    public @Dimensionless boolean isColumnVector(@Dimensionless DoubleMatrix this) {
        return columns == ((@Dimensionless int) (1));
    }

    /** Returns the diagonal of the matrix. */
    public @Dimensionless DoubleMatrix diag(@Dimensionless DoubleMatrix this) {
        assertSquare();
        @Dimensionless
        DoubleMatrix d = new @Dimensionless DoubleMatrix(rows);
        JavaBlas.rcopy(rows, data, ((@Dimensionless int) (0)), rows + ((@Dimensionless int) (1)), d.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
        return d;
    }

    /** Pretty-print this matrix to <tt>System.out</tt>. */
    public void print(@Dimensionless DoubleMatrix this) {
        System.out.println(toString());
    }

    /** Generate string representation of the matrix. */
    @Override
    public @Dimensionless String toString(@Dimensionless DoubleMatrix this) {
        return toString("%f");
    }

    /**
     * Generate string representation of the matrix, with specified
     * format for the entries. For example, <code>x.toString("%.1f")</code>
     * generates a string representations having only one position after the
     * decimal point.
     */
    public @Dimensionless String toString(@Dimensionless DoubleMatrix this, @Dimensionless String fmt) {
        return toString(fmt, "[", "]", ", ", "; ");
    }

  /**
   * Generate string representation of the matrix, with specified
   * format for the entries, and delimiters.
   *
   * @param fmt entry format (passed to String.format())
   * @param open opening parenthesis
   * @param close closing parenthesis
   * @param colSep separator between columns
   * @param rowSep separator between rows
   */
    public @Dimensionless String toString(@Dimensionless DoubleMatrix this, @Dimensionless String fmt, @Dimensionless String open, @Dimensionless String close, @Dimensionless String colSep, @Dimensionless String rowSep) {
        @Dimensionless
        StringWriter s = new @Dimensionless StringWriter();
        @Dimensionless
        PrintWriter p = new @Dimensionless PrintWriter(s);

        p.print(open);

        for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
                p.printf(fmt, get(r, c));
                if (c < columns - ((@Dimensionless int) (1))) {
                    p.print(colSep);
                }
            }
            if (r < rows - ((@Dimensionless int) (1))) {
                p.print(rowSep);
            }
        }

        p.print(close);

        return s.toString();
    }

    /** Converts the matrix to a one-dimensional array of doubles. */
    public @Dimensionless double @Dimensionless [] toArray(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        double @Dimensionless [] array = new @Dimensionless double @Dimensionless [length];

        System.arraycopy(data, ((@Dimensionless int) (0)), array, ((@Dimensionless int) (0)), length);

        return array;
    }

    /** Converts the matrix to a two-dimensional array of doubles. */
    public @Dimensionless double @Dimensionless [] @Dimensionless [] toArray2(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        double @Dimensionless [] @Dimensionless [] array = new @Dimensionless double @Dimensionless [rows] @Dimensionless [columns];

        for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
                array[r][c] = get(r, c);
            }
        }

        return array;
    }

    /** Converts the matrix to a one-dimensional array of integers. */
    public @Dimensionless int @Dimensionless [] toIntArray(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int @Dimensionless [] array = new @Dimensionless int @Dimensionless [length];

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            array[i] = (@Dimensionless int) Math.rint(get(i));
        }

        return array;
    }

    /** Convert the matrix to a two-dimensional array of integers. */
    public @Dimensionless int @Dimensionless [] @Dimensionless [] toIntArray2(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int @Dimensionless [] @Dimensionless [] array = new @Dimensionless int @Dimensionless [rows] @Dimensionless [columns];

        for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
                array[r][c] = (@Dimensionless int) Math.rint(get(r, c));
            }
        }

        return array;
    }

    /** Convert the matrix to a one-dimensional array of boolean values. */
    public @Dimensionless boolean @Dimensionless [] toBooleanArray(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        boolean @Dimensionless [] array = new @Dimensionless boolean @Dimensionless [length];

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            array[i] = get(i) != ((@Dimensionless double) (0.0));
        }

        return array;
    }

    /** Convert the matrix to a two-dimensional array of boolean values. */
    public @Dimensionless boolean @Dimensionless [] @Dimensionless [] toBooleanArray2(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        boolean @Dimensionless [] @Dimensionless [] array = new @Dimensionless boolean @Dimensionless [rows] @Dimensionless [columns];

        for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
                array[r][c] = get(r, c) != ((@Dimensionless double) (0.0));
            }
        }

        return array;
    }

    public @Dimensionless FloatMatrix toFloat(@Dimensionless DoubleMatrix this) {
         @Dimensionless
         FloatMatrix result = new @Dimensionless FloatMatrix(rows, columns);
         for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result.put(i, (@Dimensionless float) get(i));
         }
         return result;
    }

    /**
     * A wrapper which allows to view a matrix as a List of Doubles (read-only!).
     * Also implements the {@link ConvertsToDoubleMatrix} interface.
     */
    @Dimensionless
    public class ElementsAsListView extends @Dimensionless AbstractList<@Dimensionless Double> implements ConvertsToDoubleMatrix {

        private final @Dimensionless DoubleMatrix me;

        public ElementsAsListView(@Dimensionless DoubleMatrix me) {
            this.me = me;
        }

        @Override
        public @Dimensionless Double get(@Dimensionless DoubleMatrix.ElementsAsListView this, @Dimensionless int index) {
            return me.get(index);
        }

        @Override
        public @Dimensionless int size(@Dimensionless DoubleMatrix.ElementsAsListView this) {
            return me.length;
        }

        public @Dimensionless DoubleMatrix convertToDoubleMatrix(@Dimensionless DoubleMatrix.ElementsAsListView this) {
            return me;
        }
    }

    @Dimensionless
    public class RowsAsListView extends @Dimensionless AbstractList<@Dimensionless DoubleMatrix> implements ConvertsToDoubleMatrix {

        private final @Dimensionless DoubleMatrix me;

        public RowsAsListView(@Dimensionless DoubleMatrix me) {
            this.me = me;
        }

        @Override
        public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix.RowsAsListView this, @Dimensionless int index) {
            return getRow(index);
        }

        @Override
        public @Dimensionless int size(@Dimensionless DoubleMatrix.RowsAsListView this) {
            return rows;
        }

        public @Dimensionless DoubleMatrix convertToDoubleMatrix(@Dimensionless DoubleMatrix.RowsAsListView this) {
            return me;
        }
    }

    @Dimensionless
    public class ColumnsAsListView extends @Dimensionless AbstractList<@Dimensionless DoubleMatrix> implements ConvertsToDoubleMatrix {

        private final @Dimensionless DoubleMatrix me;

        public ColumnsAsListView(@Dimensionless DoubleMatrix me) {
            this.me = me;
        }

        @Override
        public @Dimensionless DoubleMatrix get(@Dimensionless DoubleMatrix.ColumnsAsListView this, @Dimensionless int index) {
            return getColumn(index);
        }

        @Override
        public @Dimensionless int size(@Dimensionless DoubleMatrix.ColumnsAsListView this) {
            return columns;
        }

        public @Dimensionless DoubleMatrix convertToDoubleMatrix(@Dimensionless DoubleMatrix.ColumnsAsListView this) {
            return me;
        }
    }

    public @Dimensionless List<@Dimensionless Double> elementsAsList(@Dimensionless DoubleMatrix this) {
        return new @Dimensionless ElementsAsListView(this);
    }

    public @Dimensionless List<@Dimensionless DoubleMatrix> rowsAsList(@Dimensionless DoubleMatrix this) {
        return new @Dimensionless RowsAsListView(this);
    }

    public @Dimensionless List<@Dimensionless DoubleMatrix> columnsAsList(@Dimensionless DoubleMatrix this) {
        return new @Dimensionless ColumnsAsListView(this);
    }

    /**************************************************************************
     * Arithmetic Operations
     */
    /**
     * Ensures that the result vector has the same length as this. If not,
     * resizing result is tried, which fails if result == this or result == other.
     */
    private void ensureResultLength(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        if (!sameLength(result)) {
            if (result == this || result == other) {
                throw new @Dimensionless SizeException("Cannot resize result matrix because it is used in-place.");
            }
            result.resize(rows, columns);
        }
    }

    /** Add two matrices (in-place). */
    public @Dimensionless DoubleMatrix addi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        if (other.isScalar()) {
            return addi(other.scalar(), result);
        }
        if (isScalar()) {
            return other.addi(scalar(), result);
        }

        assertSameLength(other);
        ensureResultLength(other, result);

        if (result == this) {
            SimpleBlas.axpy(((@Dimensionless double) (1.0)), other, result);
        } else if (result == other) {
            SimpleBlas.axpy(((@Dimensionless double) (1.0)), this, result);
        } else {
            /*SimpleBlas.copy(this, result);
            SimpleBlas.axpy(1.0, other, result);*/
            JavaBlas.rzgxpy(length, result.data, data, other.data);
        }

        return result;
    }

    /** Add a scalar to a matrix (in-place). */
    public @Dimensionless DoubleMatrix addi(@Dimensionless DoubleMatrix this, @Dimensionless double v, @Dimensionless DoubleMatrix result) {
        ensureResultLength(null, result);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result.put(i, get(i) + v);
        }
        return result;
    }

    /** Subtract two matrices (in-place). */
    public @Dimensionless DoubleMatrix subi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        if (other.isScalar()) {
            return subi(other.scalar(), result);
        }
        if (isScalar()) {
            return other.rsubi(scalar(), result);
        }

        assertSameLength(other);
        ensureResultLength(other, result);

        if (result == this) {
            SimpleBlas.axpy(- ((@Dimensionless double) (1.0)), other, result);
        } else if (result == other) {
            SimpleBlas.scal(- ((@Dimensionless double) (1.0)), result);
            SimpleBlas.axpy(((@Dimensionless double) (1.0)), this, result);
        } else {
            SimpleBlas.copy(this, result);
            SimpleBlas.axpy(- ((@Dimensionless double) (1.0)), other, result);
        }
        return result;
    }

    /** Subtract a scalar from a matrix (in-place). */
    public @Dimensionless DoubleMatrix subi(@Dimensionless DoubleMatrix this, @Dimensionless double v, @Dimensionless DoubleMatrix result) {
        ensureResultLength(null, result);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result.put(i, get(i) - v);
        }
        return result;
    }

    /**
     * Subtract two matrices, but subtract first from second matrix, that is,
     * compute <em>result = other - this</em> (in-place).
     * */
    public @Dimensionless DoubleMatrix rsubi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        return other.subi(this, result);
    }

    /** Subtract a matrix from a scalar (in-place). */
    public @Dimensionless DoubleMatrix rsubi(@Dimensionless DoubleMatrix this, @Dimensionless double a, @Dimensionless DoubleMatrix result) {
        ensureResultLength(null, result);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result.put(i, a - get(i));
        }
        return result;
    }

    /** Elementwise multiplication (in-place). */
    public @Dimensionless DoubleMatrix muli(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        if (other.isScalar()) {
            return muli(other.scalar(), result);
        }
        if (isScalar()) {
            return other.muli(scalar(), result);
        }

        assertSameLength(other);
        ensureResultLength(other, result);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result.put(i, get(i) * other.get(i));
        }
        return result;
    }

    /** Elementwise multiplication with a scalar (in-place). */
    public @Dimensionless DoubleMatrix muli(@Dimensionless DoubleMatrix this, @Dimensionless double v, @Dimensionless DoubleMatrix result) {
        ensureResultLength(null, result);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result.put(i, get(i) * v);
        }
        return result;
    }

    /** Matrix-matrix multiplication (in-place). */
    public @Dimensionless DoubleMatrix mmuli(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        if (other.isScalar()) {
            return muli(other.scalar(), result);
        }
        if (isScalar()) {
            return other.muli(scalar(), result);
        }

        /* check sizes and resize if necessary */
        assertMultipliesWith(other);
        if (result.rows != rows || result.columns != other.columns) {
            if (result != this && result != other) {
                result.resize(rows, other.columns);
            } else {
                throw new @Dimensionless SizeException("Cannot resize result matrix because it is used in-place.");
            }
        }

        if (result == this || result == other) {
            /* actually, blas cannot do multiplications in-place. Therefore, we will fake by
             * allocating a temporary object on the side and copy the result later.
             */
            @Dimensionless
            DoubleMatrix temp = new @Dimensionless DoubleMatrix(result.rows, result.columns);
            if (other.columns == ((@Dimensionless int) (1))) {
                SimpleBlas.gemv(((@Dimensionless double) (1.0)), this, other, ((@Dimensionless double) (0.0)), temp);
            } else {
                SimpleBlas.gemm(((@Dimensionless double) (1.0)), this, other, ((@Dimensionless double) (0.0)), temp);
            }
            SimpleBlas.copy(temp, result);
        } else {
            if (other.columns == ((@Dimensionless int) (1))) {
                SimpleBlas.gemv(((@Dimensionless double) (1.0)), this, other, ((@Dimensionless double) (0.0)), result);
            } else {
                SimpleBlas.gemm(((@Dimensionless double) (1.0)), this, other, ((@Dimensionless double) (0.0)), result);
            }
        }
        return result;
    }

    /** Matrix-matrix multiplication with a scalar (for symmetry, does the
     * same as <code>muli(scalar)</code> (in-place).
     */
    public @Dimensionless DoubleMatrix mmuli(@Dimensionless DoubleMatrix this, @Dimensionless double v, @Dimensionless DoubleMatrix result) {
        return muli(v, result);
    }

    /** Elementwise division (in-place). */
    public @Dimensionless DoubleMatrix divi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        if (other.isScalar()) {
            return divi(other.scalar(), result);
        }
        if (isScalar()) {
            return other.rdivi(scalar(), result);
        }

        assertSameLength(other);
        ensureResultLength(other, result);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result.put(i, get(i) / other.get(i));
        }
        return result;
    }

    /** Elementwise division with a scalar (in-place). */
    public @Dimensionless DoubleMatrix divi(@Dimensionless DoubleMatrix this, @Dimensionless double a, @Dimensionless DoubleMatrix result) {
        ensureResultLength(null, result);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result.put(i, get(i) / a);
        }
        return result;
    }

    /**
     * Elementwise division, with operands switched. Computes
     * <code>result = other / this</code> (in-place). */
    public @Dimensionless DoubleMatrix rdivi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        return other.divi(this, result);
    }

    /** (Elementwise) division with a scalar, with operands switched. Computes
     * <code>result = a / this</code> (in-place). */
    public @Dimensionless DoubleMatrix rdivi(@Dimensionless DoubleMatrix this, @Dimensionless double a, @Dimensionless DoubleMatrix result) {
        ensureResultLength(null, result);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result.put(i, a / get(i));
        }
        return result;
    }

    /** Negate each element (in-place). */
    public @Dimensionless DoubleMatrix negi(@Dimensionless DoubleMatrix this) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            put(i, -get(i));
        }
        return this;
    }

    /** Negate each element. */
    public @Dimensionless DoubleMatrix neg(@Dimensionless DoubleMatrix this) {
        return dup().negi();
    }

    /** Maps zero to 1.0 and all non-zero values to 0.0 (in-place). */
    public @Dimensionless DoubleMatrix noti(@Dimensionless DoubleMatrix this) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            put(i, get(i) == ((@Dimensionless double) (0.0)) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
        }
        return this;
    }

    /** Maps zero to 1.0 and all non-zero values to 0.0. */
    public @Dimensionless DoubleMatrix not(@Dimensionless DoubleMatrix this) {
        return dup().noti();
    }

    /** Maps zero to 0.0 and all non-zero values to 1.0 (in-place). */
    public @Dimensionless DoubleMatrix truthi(@Dimensionless DoubleMatrix this) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            put(i, get(i) == ((@Dimensionless double) (0.0)) ? ((@Dimensionless double) (0.0)) : ((@Dimensionless double) (1.0)));
        }
        return this;
    }

    /** Maps zero to 0.0 and all non-zero values to 1.0. */
    public @Dimensionless DoubleMatrix truth(@Dimensionless DoubleMatrix this) {
        return dup().truthi();
    }

    public @Dimensionless DoubleMatrix isNaNi(@Dimensionless DoubleMatrix this) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            put(i, Double.isNaN(get(i)) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
        }
        return this;
    }

    public @Dimensionless DoubleMatrix isNaN(@Dimensionless DoubleMatrix this) {
        return dup().isNaNi();
    }

    public @Dimensionless DoubleMatrix isInfinitei(@Dimensionless DoubleMatrix this) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            put(i, Double.isInfinite(get(i)) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
        }
        return this;
    }

    public @Dimensionless DoubleMatrix isInfinite(@Dimensionless DoubleMatrix this) {
        return dup().isInfinitei();
    }

    /** Checks whether all entries (i, j) with i >= j are zero. */
    public @Dimensionless boolean isLowerTriangular(@Dimensionless DoubleMatrix this) {
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows; i++)
        for (@Dimensionless int j = i+ ((@Dimensionless int) (1)); j < columns; j++) {
          if (get(i, j) != ((@Dimensionless double) (0.0)))
            return false;
        }

      return true;
    }

  /**
   * Checks whether all entries (i, j) with i <= j are zero.
   */
    public @Dimensionless boolean isUpperTriangular(@Dimensionless DoubleMatrix this) {
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows; i++)
        for (@Dimensionless int j = ((@Dimensionless int) (0)); j < i && j < columns; j++) {
          if (get(i, j) != ((@Dimensionless double) (0.0)))
            return false;
        }

      return true;
    }

    public @Dimensionless DoubleMatrix selecti(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix where) {
        checkLength(where.length);
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            if (where.get(i) == ((@Dimensionless double) (0.0))) {
                put(i, ((@Dimensionless double) (0.0)));
            }
        }
        return this;
    }

    public @Dimensionless DoubleMatrix select(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix where) {
        return dup().selecti(where);
    }

    /****************************************************************
     * Rank one-updates
     */
    /** Computes a rank-1-update A = A + alpha * x * y'. */
    public @Dimensionless DoubleMatrix rankOneUpdate(@Dimensionless DoubleMatrix this, @Dimensionless double alpha, @Dimensionless DoubleMatrix x, @Dimensionless DoubleMatrix y) {
        if (rows != x.length) {
            throw new @Dimensionless SizeException("Vector x has wrong length (" + x.length + " != " + rows + ").");
        }
        if (columns != y.length) {
            throw new @Dimensionless SizeException("Vector y has wrong length (" + x.length + " != " + columns + ").");
        }

        SimpleBlas.ger(alpha, x, y, this);
        return this;
    }

    /** Computes a rank-1-update A = A + alpha * x * x'. */
    public @Dimensionless DoubleMatrix rankOneUpdate(@Dimensionless DoubleMatrix this, @Dimensionless double alpha, @Dimensionless DoubleMatrix x) {
        return rankOneUpdate(alpha, x, x);
    }

    /** Computes a rank-1-update A = A + x * x'. */
    public @Dimensionless DoubleMatrix rankOneUpdate(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        return rankOneUpdate(((@Dimensionless double) (1.0)), x, x);
    }

    /** Computes a rank-1-update A = A + x * y'. */
    public @Dimensionless DoubleMatrix rankOneUpdate(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x, @Dimensionless DoubleMatrix y) {
        return rankOneUpdate(((@Dimensionless double) (1.0)), x, y);
    }

    /****************************************************************
     * Logical operations
     */
    /** Returns the minimal element of the matrix. */
    public @Dimensionless double min(@Dimensionless DoubleMatrix this) {
        if (isEmpty()) {
            return ((@Dimensionless double) (Double.POSITIVE_INFINITY));
        }
        @Dimensionless
        double v = ((@Dimensionless double) (Double.POSITIVE_INFINITY));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            if (!Double.isNaN(get(i)) && get(i) < v) {
                v = get(i);
            }
        }

        return v;
    }

    /**
     * Returns the linear index of the minimal element. If there are
     * more than one elements with this value, the first one is returned.
     */
    public @Dimensionless int argmin(@Dimensionless DoubleMatrix this) {
        if (isEmpty()) {
            return ((@Dimensionless int) (-1));
        }
        @Dimensionless
        double v = ((@Dimensionless double) (Double.POSITIVE_INFINITY));
        @Dimensionless
        int a = ((@Dimensionless int) (-1));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            if (!Double.isNaN(get(i)) && get(i) < v) {
                v = get(i);
                a = i;
            }
        }

        return a;
    }

    /**
     * Computes the minimum between two matrices. Returns the smaller of the
     * corresponding elements in the matrix (in-place).
     */
    public @Dimensionless DoubleMatrix mini(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        if (result == this) {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                if (get(i) > other.get(i)) {
                    put(i, other.get(i));
                }
            }
        } else {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                if (get(i) > other.get(i)) {
                    result.put(i, other.get(i));
                } else {
                    result.put(i, get(i));
                }
            }
        }
        return result;
    }

    /**
     * Computes the minimum between two matrices. Returns the smaller of the
     * corresponding elements in the matrix (in-place on this).
     */
    public @Dimensionless DoubleMatrix mini(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
        return mini(other, this);
    }

    /**
     * Computes the minimum between two matrices. Returns the smaller of the
     * corresponding elements in the matrix (in-place on this).
     */
    public @Dimensionless DoubleMatrix min(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
        return mini(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    public @Dimensionless DoubleMatrix mini(@Dimensionless DoubleMatrix this, @Dimensionless double v, @Dimensionless DoubleMatrix result) {
        if (result == this) {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                if (get(i) > v) {
                    result.put(i, v);
                }
            }
        } else {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                if (get(i) > v) {
                    result.put(i, v);
                } else {
                    result.put(i, get(i));
                }
            }

        }
        return result;
    }

    public @Dimensionless DoubleMatrix mini(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
        return mini(v, this);
    }

    public @Dimensionless DoubleMatrix min(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
        return mini(v, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Returns the maximal element of the matrix. */
    public @Dimensionless double max(@Dimensionless DoubleMatrix this) {
        if (isEmpty()) {
            return ((@Dimensionless double) (Double.NEGATIVE_INFINITY));
        }
        @Dimensionless
        double v = ((@Dimensionless double) (Double.NEGATIVE_INFINITY));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            if (!Double.isNaN(get(i)) && get(i) > v) {
                v = get(i);
            }
        }
        return v;
    }

    /**
     * Returns the linear index of the maximal element of the matrix. If
     * there are more than one elements with this value, the first one
     * is returned.
     */
    public @Dimensionless int argmax(@Dimensionless DoubleMatrix this) {
        if (isEmpty()) {
            return ((@Dimensionless int) (-1));
        }
        @Dimensionless
        double v = ((@Dimensionless double) (Double.NEGATIVE_INFINITY));
        @Dimensionless
        int a = ((@Dimensionless int) (-1));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            if (!Double.isNaN(get(i)) && get(i) > v) {
                v = get(i);
                a = i;
            }
        }

        return a;
    }

    /**
     * Computes the maximum between two matrices. Returns the larger of the
     * corresponding elements in the matrix (in-place).
     */
    public @Dimensionless DoubleMatrix maxi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
        if (result == this) {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                if (get(i) < other.get(i)) {
                    put(i, other.get(i));
                }
            }
        } else {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                if (get(i) < other.get(i)) {
                    result.put(i, other.get(i));
                } else {
                    result.put(i, get(i));
                }
            }
        }
        return result;
    }

    /**
     * Computes the maximum between two matrices. Returns the smaller of the
     * corresponding elements in the matrix (in-place on this).
     */
    public @Dimensionless DoubleMatrix maxi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
        return maxi(other, this);
    }

    /**
     * Computes the maximum between two matrices. Returns the larger of the
     * corresponding elements in the matrix (in-place on this).
     */
    public @Dimensionless DoubleMatrix max(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
        return maxi(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    public @Dimensionless DoubleMatrix maxi(@Dimensionless DoubleMatrix this, @Dimensionless double v, @Dimensionless DoubleMatrix result) {
        if (result == this) {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                if (get(i) < v) {
                    result.put(i, v);
                }
            }
        } else {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                if (get(i) < v) {
                    result.put(i, v);
                } else {
                    result.put(i, get(i));
                }
            }

        }
        return result;
    }

    public @Dimensionless DoubleMatrix maxi(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
        return maxi(v, this);
    }

    public @Dimensionless DoubleMatrix max(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
        return maxi(v, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Computes the sum of all elements of the matrix. */
    public @Dimensionless double sum(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        double s = ((@Dimensionless double) (0.0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            s += get(i);
        }
        return s;
    }

    /** Computes the product of all elements of the matrix */
    public @Dimensionless double prod(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        double p = ((@Dimensionless double) (1.0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            p *= get(i);
        }
        return p;
    }

    /**
     * Computes the mean value of all elements in the matrix,
     * that is, <code>x.sum() / x.length</code>.
     */
    public double mean(@Dimensionless DoubleMatrix this) {
        return sum() / length;
    }

    /**
     * Computes the cumulative sum, that is, the sum of all elements
     * of the matrix up to a given index in linear addressing (in-place).
     */
    public @Dimensionless DoubleMatrix cumulativeSumi(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        double s = ((@Dimensionless double) (0.0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            s += get(i);
            put(i, s);
        }
        return this;
    }

    /**
     * Computes the cumulative sum, that is, the sum of all elements
     * of the matrix up to a given index in linear addressing.
     */
    public @Dimensionless DoubleMatrix cumulativeSum(@Dimensionless DoubleMatrix this) {
        return dup().cumulativeSumi();
    }

    /** The scalar product of this with other. */
    public @Dimensionless double dot(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
        return SimpleBlas.dot(this, other);
    }

    /** 
     * Computes the projection coefficient of other on this.
     *
     * The returned scalar times <tt>this</tt> is the orthogonal projection
     * of <tt>other</tt> on <tt>this</tt>.
     */
    public @Dimensionless double project(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
        other.checkLength(length);
        double norm = ((@Dimensionless int) (0)), dot = ((@Dimensionless int) (0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < this.length; i++) {
            @Dimensionless
            double x = get(i);
            norm += x * x;
            dot += x * other.get(i);
        }
        return dot / norm;
    }

    /**
     * The Euclidean norm of the matrix as vector, also the Frobenius
     * norm of the matrix.
     */
    public double norm2(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        double norm = ((@Dimensionless double) (0.0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            norm += get(i) * get(i);
        }
        return (@Dimensionless double) Math.sqrt(norm);
    }

    /**
     * The maximum norm of the matrix (maximal absolute value of the elements).
     */
    public @Dimensionless double normmax(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        double max = ((@Dimensionless double) (0.0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            @Dimensionless
            double a = Math.abs(get(i));
            if (a > max) {
                max = a;
            }
        }
        return max;
    }

    /**
     * The 1-norm of the matrix as vector (sum of absolute values of elements).
     */
    public @Dimensionless double norm1(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        double norm = ((@Dimensionless double) (0.0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            norm += Math.abs(get(i));
        }
        return norm;
    }

    /**
     * Returns the squared (Euclidean) distance.
     */
    public @Dimensionless double squaredDistance(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
        other.checkLength(length);
        @Dimensionless
        double sd = ((@Dimensionless double) (0.0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            @Dimensionless
            double d = get(i) - other.get(i);
            sd += d * d;
        }
        return sd;
    }

    /**
     * Returns the (euclidean) distance.
     */
    public @Dimensionless double distance2(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
        return (@Dimensionless double) Math.sqrt(squaredDistance(other));
    }

    /**
     * Returns the (1-norm) distance.
     */
    public @Dimensionless double distance1(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
        other.checkLength(length);
        @Dimensionless
        double d = ((@Dimensionless double) (0.0));
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            d += Math.abs(get(i) - other.get(i));
        }
        return d;
    }

    /**
     * Return a new matrix with all elements sorted.
     */
    public @Dimensionless DoubleMatrix sort(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        double array @Dimensionless [] = toArray();
        java.util.Arrays.sort(array);
        return new @Dimensionless DoubleMatrix(rows, columns, array);
    }

    /**
     * Sort elements in-place.
     */
    public @Dimensionless DoubleMatrix sorti(@Dimensionless DoubleMatrix this) {
        Arrays.sort(data);
        return this;
    }

    /**
     * Get the sorting permutation.
     *
     * @return an int[] array such that which indexes the elements in sorted
     * order.
     */
    public @Dimensionless int @Dimensionless [] sortingPermutation(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        Integer @Dimensionless [] indices = new @Dimensionless Integer @Dimensionless [length];

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            indices[i] = i;
        }

        final @Dimensionless double @Dimensionless [] array = data;

        Arrays.sort(indices, new @Dimensionless Comparator() @Dimensionless {

            public @Dimensionless int compare(@Dimensionless Object o1, @Dimensionless Object o2) {
                @Dimensionless
                int i = (@Dimensionless Integer) o1;
                @Dimensionless
                int j = (@Dimensionless Integer) o2;
                if (array[i] < array[j]) {
                    return ((@Dimensionless int) (-1));
                } else if (array[i] == array[j]) {
                    return ((@Dimensionless int) (0));
                } else {
                    return ((@Dimensionless int) (1));
                }
            }
        });

        @Dimensionless
        int @Dimensionless [] result = new @Dimensionless int @Dimensionless [length];

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
            result[i] = indices[i];
        }

        return result;
    }

    /**
     * Sort columns (in-place).
     */
    public @Dimensionless DoubleMatrix sortColumnsi(@Dimensionless DoubleMatrix this) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i += rows) {
            Arrays.sort(data, i, i + rows);
        }
        return this;
    }

    /** Sort columns. */
    public @Dimensionless DoubleMatrix sortColumns(@Dimensionless DoubleMatrix this) {
        return dup().sortColumnsi();
    }

    /** Return matrix of indices which sort all columns. */
    public @Dimensionless int @Dimensionless [] @Dimensionless [] columnSortingPermutations(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int @Dimensionless [] @Dimensionless [] result = new @Dimensionless int @Dimensionless [columns] @Dimensionless [];

        @Dimensionless
        DoubleMatrix temp = new @Dimensionless DoubleMatrix(rows);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            result[c] = getColumn(c, temp).sortingPermutation();
        }

        return result;
    }

    /** Sort rows (in-place). */
    public @Dimensionless DoubleMatrix sortRowsi(@Dimensionless DoubleMatrix this) {
        // actually, this is much harder because the data is not consecutive
        // in memory...
        @Dimensionless
        DoubleMatrix temp = new @Dimensionless DoubleMatrix(columns);
        for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
            putRow(r, getRow(r, temp).sorti());
        }
        return this;
    }

    /** Sort rows. */
    public @Dimensionless DoubleMatrix sortRows(@Dimensionless DoubleMatrix this) {
        return dup().sortRowsi();
    }

    /** Return matrix of indices which sort all columns. */
    public @Dimensionless int @Dimensionless [] @Dimensionless [] rowSortingPermutations(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int @Dimensionless [] @Dimensionless [] result = new @Dimensionless int @Dimensionless [rows] @Dimensionless [];

        @Dimensionless
        DoubleMatrix temp = new @Dimensionless DoubleMatrix(columns);
        for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
            result[r] = getRow(r, temp).sortingPermutation();
        }

        return result;
    }

    /** Return a vector containing the sums of the columns (having number of columns many entries) */
    public DoubleMatrix columnSums(@Dimensionless DoubleMatrix this) {
        if (rows == ((@Dimensionless int) (1))) {
            return dup();
        } else {
            @Dimensionless
            DoubleMatrix v = new @Dimensionless DoubleMatrix(((@Dimensionless int) (1)), columns);

            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
                for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                    v.put(c, v.get(c) + get(r, c));
                }
            }

            return v;
        }
    }

    /** Return a vector containing the means of all columns. */
    public @Dimensionless DoubleMatrix columnMeans(@Dimensionless DoubleMatrix this) {
        return columnSums().divi(rows);
    }

    /** Return a vector containing the sum of the rows. */
    public @Dimensionless DoubleMatrix rowSums(@Dimensionless DoubleMatrix this) {
        if (columns == ((@Dimensionless int) (1))) {
            return dup();
        } else {
            @Dimensionless
            DoubleMatrix v = new @Dimensionless DoubleMatrix(rows);

            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
                for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                    v.put(r, v.get(r) + get(r, c));
                }
            }

            return v;
        }
    }

    /** Return a vector containing the means of the rows. */
    public @Dimensionless DoubleMatrix rowMeans(@Dimensionless DoubleMatrix this) {
        return rowSums().divi(columns);
    }

    /************************************************************************
     * Column and rows access.
     */

    /** Get a copy of a column. */
    public @Dimensionless DoubleMatrix getColumn(@Dimensionless DoubleMatrix this, @Dimensionless int c) {
        return getColumn(c, new @Dimensionless DoubleMatrix(rows, ((@Dimensionless int) (1))));
    }

    /** Copy a column to the given vector. */
    public DoubleMatrix getColumn(@Dimensionless DoubleMatrix this, int c, DoubleMatrix result) {
        result.checkLength(rows);
        JavaBlas.rcopy(rows, data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)), result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
        return result;
    }

    /** Copy a column back into the matrix. */
    public void putColumn(@Dimensionless DoubleMatrix this, int c, DoubleMatrix v) {
        JavaBlas.rcopy(rows, v.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)));
    }

    /** Get a copy of a row. */
    public @Dimensionless DoubleMatrix getRow(@Dimensionless DoubleMatrix this, @Dimensionless int r) {
        return getRow(r, new @Dimensionless DoubleMatrix(((@Dimensionless int) (1)), columns));
    }

    /** Copy a row to a given vector. */
    public DoubleMatrix getRow(@Dimensionless DoubleMatrix this, int r, DoubleMatrix result) {
        result.checkLength(columns);
        JavaBlas.rcopy(columns, data, index(r, ((@Dimensionless int) (0))), rows, result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
        return result;
    }

    /** Copy a row back into the matrix. */
    public void putRow(@Dimensionless DoubleMatrix this, int r, DoubleMatrix v) {
        JavaBlas.rcopy(columns, v.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(r, ((@Dimensionless int) (0))), rows);
    }

    /** Return column-wise minimums. */
    public @Dimensionless DoubleMatrix columnMins(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        DoubleMatrix mins = new @Dimensionless DoubleMatrix(((@Dimensionless int) (1)), columns);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            mins.put(c, getColumn(c).min());
        }
        return mins;
    }

    /** Return index of minimal element per column. */
    public @Dimensionless int @Dimensionless [] columnArgmins(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int @Dimensionless [] argmins = new @Dimensionless int @Dimensionless [columns];
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            argmins[c] = getColumn(c).argmin();
        }
        return argmins;
    }

    /** Return column-wise maximums. */
    public @Dimensionless DoubleMatrix columnMaxs(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        DoubleMatrix maxs = new @Dimensionless DoubleMatrix(((@Dimensionless int) (1)), columns);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            maxs.put(c, getColumn(c).max());
        }
        return maxs;
    }

    /** Return index of minimal element per column. */
    public @Dimensionless int @Dimensionless [] columnArgmaxs(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int @Dimensionless [] argmaxs = new @Dimensionless int @Dimensionless [columns];
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            argmaxs[c] = getColumn(c).argmax();
        }
        return argmaxs;
    }

    /** Return row-wise minimums. */
    public @Dimensionless DoubleMatrix rowMins(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        DoubleMatrix mins = new @Dimensionless DoubleMatrix(rows);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < rows; c++) {
            mins.put(c, getRow(c).min());
        }
        return mins;
    }

    /** Return index of minimal element per row. */
    public @Dimensionless int @Dimensionless [] rowArgmins(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int @Dimensionless [] argmins = new @Dimensionless int @Dimensionless [rows];
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < rows; c++) {
            argmins[c] = getRow(c).argmin();
        }
        return argmins;
    }

    /** Return row-wise maximums. */
    public @Dimensionless DoubleMatrix rowMaxs(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        DoubleMatrix maxs = new @Dimensionless DoubleMatrix(rows);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < rows; c++) {
            maxs.put(c, getRow(c).max());
        }
        return maxs;
    }

    /** Return index of minimal element per row. */
    public @Dimensionless int @Dimensionless [] rowArgmaxs(@Dimensionless DoubleMatrix this) {
        @Dimensionless
        int @Dimensionless [] argmaxs = new @Dimensionless int @Dimensionless [rows];
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < rows; c++) {
            argmaxs[c] = getRow(c).argmax();
        }
        return argmaxs;
    }

    /**************************************************************************
     * Elementwise Functions
     */
    /** Add a row vector to all rows of the matrix (in place). */
    public DoubleMatrix addiRowVector(@Dimensionless DoubleMatrix this, DoubleMatrix x) {
        x.checkLength(columns);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                put(r, c, get(r, c) + x.get(c));
            }
        }
        return this;
    }

    /** Add a row to all rows of the matrix. */
    public @Dimensionless DoubleMatrix addRowVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        return dup().addiRowVector(x);
    }

    /** Add a vector to all columns of the matrix (in-place). */
    public DoubleMatrix addiColumnVector(@Dimensionless DoubleMatrix this, DoubleMatrix x) {
        x.checkLength(rows);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                put(r, c, get(r, c) + x.get(r));
            }
        }
        return this;
    }

    /** Add a vector to all columns of the matrix. */
    public @Dimensionless DoubleMatrix addColumnVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        return dup().addiColumnVector(x);
    }

    /** Subtract a row vector from all rows of the matrix (in-place). */
    public @Dimensionless DoubleMatrix subiRowVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        // This is a bit crazy, but a row vector must have as length as the columns of the matrix.
        x.checkLength(columns);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                put(r, c, get(r, c) - x.get(c));
            }
        }
        return this;
    }

    /** Subtract a row vector from all rows of the matrix. */
    public @Dimensionless DoubleMatrix subRowVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        return dup().subiRowVector(x);
    }

    /** Subtract a column vector from all columns of the matrix (in-place). */
    public @Dimensionless DoubleMatrix subiColumnVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        x.checkLength(rows);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                put(r, c, get(r, c) - x.get(r));
            }
        }
        return this;
    }

    /** Subtract a vector from all columns of the matrix. */
    public @Dimensionless DoubleMatrix subColumnVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        return dup().subiColumnVector(x);
    }

    /** Multiply a row by a scalar. */
    public @Dimensionless DoubleMatrix mulRow(@Dimensionless DoubleMatrix this, @Dimensionless int r, @Dimensionless double scale) {
        NativeBlas.dscal(columns, scale, data, index(r, ((@Dimensionless int) (0))), rows);
        return this;
    }

    /** Multiply a column by a scalar. */
    public @Dimensionless DoubleMatrix mulColumn(@Dimensionless DoubleMatrix this, @Dimensionless int c, @Dimensionless double scale) {
        NativeBlas.dscal(rows, scale, data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)));
        return this;
    }

    /** Multiply all columns with a column vector (in-place). */
    public @Dimensionless DoubleMatrix muliColumnVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        x.checkLength(rows);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                put(r, c, get(r, c) * x.get(r));
            }
        }
        return this;
    }

    /** Multiply all columns with a column vector. */
    public @Dimensionless DoubleMatrix mulColumnVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        return dup().muliColumnVector(x);
    }

    /** Multiply all rows with a row vector (in-place). */
    public @Dimensionless DoubleMatrix muliRowVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        x.checkLength(columns);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                put(r, c, get(r, c) * x.get(c));
            }
        }
        return this;
    }

    /** Multiply all rows with a row vector. */
    public @Dimensionless DoubleMatrix mulRowVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        return dup().muliRowVector(x);
    }

    public @Dimensionless DoubleMatrix diviRowVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        x.checkLength(columns);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                put(r, c, get(r, c) / x.get(c));
            }
        }
        return this;
    }

    public @Dimensionless DoubleMatrix divRowVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        return dup().diviRowVector(x);
    }

    public @Dimensionless DoubleMatrix diviColumnVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        x.checkLength(rows);
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
            for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
                put(r, c, get(r, c) / x.get(r));
            }
        }
        return this;
    }

    public @Dimensionless DoubleMatrix divColumnVector(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix x) {
        return dup().diviColumnVector(x);
    }

    /**
     * Writes out this matrix to the given data stream.
     * @param dos the data output stream to write to.
     * @throws IOException
     */
    public void out(@Dimensionless DoubleMatrix this, @Dimensionless DataOutputStream dos) throws IOException {
        dos.writeUTF("double");
        dos.writeInt(columns);
        dos.writeInt(rows);

        dos.writeInt(data.length);
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < data.length; i++) {
            dos.writeDouble(data[i]);
        }
    }

    /**
     * Reads in a matrix from the given data stream. Note
     * that the old data of this matrix will be discarded.
     * @param dis the data input stream to read from.
     * @throws IOException
     */
    public void in(@Dimensionless DataInputStream dis) throws IOException {
        if (!dis.readUTF().equals("double")) {
            throw new @Dimensionless IllegalStateException("The matrix in the specified file is not of the correct type!");
        }

        this.columns = dis.readInt();
        this.rows = dis.readInt();

        final @Dimensionless int MAX = dis.readInt();
        data = new @Dimensionless double @Dimensionless [MAX];
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < MAX; i++) {
            data[i] = dis.readDouble();
        }
    }

    /**
     * Saves this matrix to the specified file.
     * @param filename the file to write the matrix in.
     * @throws IOException thrown on errors while writing the matrix to the file
     */
    public void save(@Dimensionless DoubleMatrix this, @Dimensionless String filename) throws IOException {
        @Dimensionless
        FileOutputStream fos = new @Dimensionless FileOutputStream(filename, false);
        @Dimensionless
        DataOutputStream dos = new @Dimensionless DataOutputStream(fos);
        try {
            this.out(dos);
        } finally {
            dos.close();
            fos.close();
        }
    }

    /**
     * Loads a matrix from a file into this matrix. Note that the old data
     * of this matrix will be discarded.
     * @param filename the file to read the matrix from
     * @throws IOException thrown on errors while reading the matrix
     */
    public void load(@Dimensionless DoubleMatrix this, @Dimensionless String filename) throws IOException {
        @Dimensionless
        FileInputStream fis = new @Dimensionless FileInputStream(filename);
        @Dimensionless
        DataInputStream dis = new @Dimensionless DataInputStream(fis);
        try {
            this.in(dis);
        }
        finally {
            dis.close();
            fis.close();
        }
    }

    public static @Dimensionless DoubleMatrix loadAsciiFile(@Dimensionless String filename) throws IOException {
        @Dimensionless
        BufferedReader is = new @Dimensionless BufferedReader(new @Dimensionless InputStreamReader(new @Dimensionless FileInputStream(filename)));

        // Go through file and count columns and rows. What makes this endeavour a bit difficult is
        // that files can have leading or trailing spaces leading to spurious fields
        // after String.split().
        @Dimensionless
        String line;
        @Dimensionless
        int rows = ((@Dimensionless int) (0));
        @Dimensionless
        int columns = ((@Dimensionless int) (-1));
        while ((line = is.readLine()) != null) {
            @Dimensionless
            String @Dimensionless [] elements = WHITESPACES.split(line);
            @Dimensionless
            int numElements = elements.length;
            if (elements[((@Dimensionless int) (0))].length() == ((@Dimensionless int) (0))) {
                numElements--;
            }
            if (elements[elements.length - ((@Dimensionless int) (1))].length() == ((@Dimensionless int) (0))) {
                numElements--;
            }

            if (columns == ((@Dimensionless int) (-1))) {
                columns = numElements;
            } else {
                if (columns != numElements) {
                    throw new @Dimensionless IOException("Number of elements changes in line " + line + ".");
                }
            }

            rows++;
        }
        is.close();

        @Dimensionless
        FileInputStream fis = new @Dimensionless FileInputStream(filename);
        try {
            // Go through file a second time process the actual data.
            is = new @Dimensionless BufferedReader(new @Dimensionless InputStreamReader(fis));
            @Dimensionless
            DoubleMatrix result = new @Dimensionless DoubleMatrix(rows, columns);
            @Dimensionless
            int r = ((@Dimensionless int) (0));
            while ((line = is.readLine()) != null) {
                @Dimensionless
                String @Dimensionless [] elements = WHITESPACES.split(line);
                @Dimensionless
                int firstElement = (elements[((@Dimensionless int) (0))].length() == ((@Dimensionless int) (0))) ? ((@Dimensionless int) (1)) : ((@Dimensionless int) (0));
                for (int c = ((@Dimensionless int) (0)), cc = firstElement; c < columns; c++, cc++) {
                    result.put(r, c, Double.valueOf(elements[cc]));
                }
                r++;
            }
            return result;
        } finally {
            fis.close();
        }
    }

    public static @Dimensionless DoubleMatrix loadCSVFile(@Dimensionless String filename) throws IOException {
        @Dimensionless
        BufferedReader is = new @Dimensionless BufferedReader(new @Dimensionless InputStreamReader(new @Dimensionless FileInputStream(filename)));

        @Dimensionless
        List<@Dimensionless DoubleMatrix> rows = new @Dimensionless LinkedList<@Dimensionless DoubleMatrix>();
        @Dimensionless
        String line;
        @Dimensionless
        int columns = ((@Dimensionless int) (-1));
        while ((line = is.readLine()) != null) {
            @Dimensionless
            String @Dimensionless [] elements = COMMA.split(line);
            @Dimensionless
            int numElements = elements.length;
            if (elements[((@Dimensionless int) (0))].length() == ((@Dimensionless int) (0))) {
                numElements--;
            }
            if (elements[elements.length - ((@Dimensionless int) (1))].length() == ((@Dimensionless int) (0))) {
                numElements--;
            }

            if (columns == ((@Dimensionless int) (-1))) {
                columns = numElements;
            } else {
                if (columns != numElements) {
                    throw new @Dimensionless IOException("Number of elements changes in line " + line + ".");
                }
            }

            @Dimensionless
            DoubleMatrix row = new @Dimensionless DoubleMatrix(columns);
            for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
                row.put(c, Double.valueOf(elements[c]));
            }
            rows.add(row);
        }
        is.close();

        System.out.println("Done reading file");

        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(rows.size(), columns);
        @Dimensionless
        int r = ((@Dimensionless int) (0));
        @Dimensionless
        Iterator<@Dimensionless DoubleMatrix> ri = rows.iterator();
        while (ri.hasNext()) {
            result.putRow(r, ri.next());
            r++;
        }
        return result;
    }

    /****************************************************************
     * Autogenerated code
     */
    /***** Code for operators ***************************************/

    /* Overloads for the usual arithmetic operations */
    /*#
    def gen_overloads(base, result_rows, result_cols, verb=''); <<-EOS
    #{doc verb.capitalize + " a matrix (in place)."}
    public DoubleMatrix #{base}i(DoubleMatrix other) {
      return #{base}i(other, this);
    }

    #{doc verb.capitalize + " a matrix."}
    public DoubleMatrix #{base}(DoubleMatrix other) {
      return #{base}i(other, new DoubleMatrix(#{result_rows}, #{result_cols}));
    }

    #{doc verb.capitalize + " a scalar (in place)."}
    public DoubleMatrix #{base}i(double v) {
      return #{base}i(v, this);
    }

    #{doc verb.capitalize + " a scalar."}
    public DoubleMatrix #{base}(double v) {
      return #{base}i(v, new DoubleMatrix(rows, columns));
    }
    EOS
    end
    #*/

    /* Generating code for logical operators. This not only generates the stubs
     * but really all of the code.
     */
    /*#
    def gen_compare(name, op, cmp); <<-EOS
    #{doc 'Test for ' + cmp + ' (in-place).'}
    public DoubleMatrix #{name}i(DoubleMatrix other, DoubleMatrix result) {
      if (other.isScalar())
        return #{name}i(other.scalar(), result);

      assertSameLength(other);
      ensureResultLength(other, result);

      for (int i = 0; i < length; i++)
        result.put(i, get(i) #{op} other.get(i) ? 1.0 : 0.0);
      return result;
    }

    #{doc 'Test for ' + cmp + ' (in-place).'}
    public DoubleMatrix #{name}i(DoubleMatrix other) {
      return #{name}i(other, this);
    }

    #{doc 'Test for ' + cmp + '.'}
    public DoubleMatrix #{name}(DoubleMatrix other) {
      return #{name}i(other, new DoubleMatrix(rows, columns));
    }

    #{doc 'Test for ' + cmp + ' against a scalar (in-place).'}
    public DoubleMatrix #{name}i(double value, DoubleMatrix result) {
      ensureResultLength(null, result);
      for (int i = 0; i < length; i++)
        result.put(i, get(i) #{op} value ? 1.0 : 0.0);
      return result;
    }

    #{doc 'Test for ' + cmp + ' against a scalar (in-place).'}
    public DoubleMatrix #{name}i(double value) {
      return #{name}i(value, this);
    }

    #{doc 'test for ' + cmp + ' against a scalar.'}
    public DoubleMatrix #{name}(double value) {
      return #{name}i(value, new DoubleMatrix(rows, columns));
    }
    EOS
    end
    #*/
    /*#
    def gen_logical(name, op, cmp); <<-EOS
    #{doc 'Compute elementwise ' + cmp + ' (in-place).'}
    public DoubleMatrix #{name}i(DoubleMatrix other, DoubleMatrix result) {
      assertSameLength(other);
      ensureResultLength(other, result);

      for (int i = 0; i < length; i++)
        result.put(i, (get(i) != 0.0) #{op} (other.get(i) != 0.0) ? 1.0 : 0.0);
      return result;
    }

    #{doc 'Compute elementwise ' + cmp + ' (in-place).'}
    public DoubleMatrix #{name}i(DoubleMatrix other) {
      return #{name}i(other, this);
    }

    #{doc 'Compute elementwise ' + cmp + '.'}
    public DoubleMatrix #{name}(DoubleMatrix other) {
      return #{name}i(other, new DoubleMatrix(rows, columns));
    }

    #{doc 'Compute elementwise ' + cmp + ' against a scalar (in-place).'}
    public DoubleMatrix #{name}i(double value, DoubleMatrix result) {
      ensureResultLength(null, result);
      boolean val = (value != 0.0);
      for (int i = 0; i < length; i++)
        result.put(i, (get(i) != 0.0) #{op} val ? 1.0 : 0.0);
      return result;
    }

    #{doc 'Compute elementwise ' + cmp + ' against a scalar (in-place).'}
    public DoubleMatrix #{name}i(double value) {
      return #{name}i(value, this);
    }

    #{doc 'Compute elementwise ' + cmp + ' against a scalar.'}
    public DoubleMatrix #{name}(double value) {
      return #{name}i(value, new DoubleMatrix(rows, columns));
    }
    EOS
    end
    #*/

    /*# collect(gen_overloads('add', 'rows', 'columns', 'add'),
    gen_overloads('sub', 'rows', 'columns', 'subtract'),
    gen_overloads('rsub', 'rows', 'columns', '(right-)subtract'),
    gen_overloads('div', 'rows', 'columns', 'elementwise divide by'),
    gen_overloads('rdiv', 'rows', 'columns', '(right-)elementwise divide by'),
    gen_overloads('mul', 'rows', 'columns', 'elementwise multiply by'),
    gen_overloads('mmul', 'rows', 'other.columns', 'matrix-multiply by'),
    gen_compare('lt', '<', '"less than"'),
    gen_compare('gt', '>', '"greater than"'),
    gen_compare('le', '<=', '"less than or equal"'),
    gen_compare('ge', '>=', '"greater than or equal"'),
    gen_compare('eq', '==', 'equality'),
    gen_compare('ne', '!=', 'inequality'),
    gen_logical('and', '&', 'logical and'),
    gen_logical('or', '|', 'logical or'),
    gen_logical('xor', '^', 'logical xor'))
    #*/
//RJPP-BEGIN------------------------------------------------------------
    /** Add a matrix (in place). */
    public @Dimensionless DoubleMatrix addi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return addi(other, this);
    }

    /** Add a matrix. */
    public @Dimensionless DoubleMatrix add(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return addi(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Add a scalar (in place). */
    public @Dimensionless DoubleMatrix addi(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return addi(v, this);
    }

    /** Add a scalar. */
    public @Dimensionless DoubleMatrix add(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return addi(v, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Subtract a matrix (in place). */
    public @Dimensionless DoubleMatrix subi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return subi(other, this);
    }

    /** Subtract a matrix. */
    public @Dimensionless DoubleMatrix sub(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return subi(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Subtract a scalar (in place). */
    public DoubleMatrix subi(@Dimensionless DoubleMatrix this, double v) {
      return subi(v, this);
    }

    /** Subtract a scalar. */
    public @Dimensionless DoubleMatrix sub(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return subi(v, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** (right-)subtract a matrix (in place). */
    public @Dimensionless DoubleMatrix rsubi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return rsubi(other, this);
    }

    /** (right-)subtract a matrix. */
    public @Dimensionless DoubleMatrix rsub(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return rsubi(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** (right-)subtract a scalar (in place). */
    public @Dimensionless DoubleMatrix rsubi(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return rsubi(v, this);
    }

    /** (right-)subtract a scalar. */
    public @Dimensionless DoubleMatrix rsub(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return rsubi(v, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Elementwise divide by a matrix (in place). */
    public @Dimensionless DoubleMatrix divi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return divi(other, this);
    }

    /** Elementwise divide by a matrix. */
    public @Dimensionless DoubleMatrix div(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return divi(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Elementwise divide by a scalar (in place). */
    public DoubleMatrix divi(@Dimensionless DoubleMatrix this, double v) {
      return divi(v, this);
    }

    /** Elementwise divide by a scalar. */
    public @Dimensionless DoubleMatrix div(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return divi(v, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** (right-)elementwise divide by a matrix (in place). */
    public @Dimensionless DoubleMatrix rdivi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return rdivi(other, this);
    }

    /** (right-)elementwise divide by a matrix. */
    public @Dimensionless DoubleMatrix rdiv(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return rdivi(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** (right-)elementwise divide by a scalar (in place). */
    public @Dimensionless DoubleMatrix rdivi(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return rdivi(v, this);
    }

    /** (right-)elementwise divide by a scalar. */
    public @Dimensionless DoubleMatrix rdiv(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return rdivi(v, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Elementwise multiply by a matrix (in place). */
    public @Dimensionless DoubleMatrix muli(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return muli(other, this);
    }

    /** Elementwise multiply by a matrix. */
    public DoubleMatrix mul(@Dimensionless DoubleMatrix this, DoubleMatrix other) {
      return muli(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Elementwise multiply by a scalar (in place). */
    public DoubleMatrix muli(@Dimensionless DoubleMatrix this, double v) {
      return muli(v, this);
    }

    /** Elementwise multiply by a scalar. */
    public @Dimensionless DoubleMatrix mul(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return muli(v, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Matrix-multiply by a matrix (in place). */
    public @Dimensionless DoubleMatrix mmuli(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return mmuli(other, this);
    }

    /** Matrix-multiply by a matrix. */
    public DoubleMatrix mmul(@Dimensionless DoubleMatrix this, DoubleMatrix other) {
      return mmuli(other, new @Dimensionless DoubleMatrix(rows, other.columns));
    }

    /** Matrix-multiply by a scalar (in place). */
    public @Dimensionless DoubleMatrix mmuli(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return mmuli(v, this);
    }

    /** Matrix-multiply by a scalar. */
    public @Dimensionless DoubleMatrix mmul(@Dimensionless DoubleMatrix this, @Dimensionless double v) {
      return mmuli(v, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for "less than" (in-place). */
    public @Dimensionless DoubleMatrix lti(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
      if (other.isScalar())
        return lti(other.scalar(), result);

      assertSameLength(other);
      ensureResultLength(other, result);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) < other.get(i) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for "less than" (in-place). */
    public @Dimensionless DoubleMatrix lti(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return lti(other, this);
    }

    /** Test for "less than". */
    public @Dimensionless DoubleMatrix lt(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return lti(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for "less than" against a scalar (in-place). */
    public @Dimensionless DoubleMatrix lti(@Dimensionless DoubleMatrix this, @Dimensionless double value, @Dimensionless DoubleMatrix result) {
      ensureResultLength(null, result);
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) < value ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for "less than" against a scalar (in-place). */
    public @Dimensionless DoubleMatrix lti(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return lti(value, this);
    }

    /** test for "less than" against a scalar. */
    public @Dimensionless DoubleMatrix lt(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return lti(value, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for "greater than" (in-place). */
    public @Dimensionless DoubleMatrix gti(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
      if (other.isScalar())
        return gti(other.scalar(), result);

      assertSameLength(other);
      ensureResultLength(other, result);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) > other.get(i) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for "greater than" (in-place). */
    public @Dimensionless DoubleMatrix gti(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return gti(other, this);
    }

    /** Test for "greater than". */
    public @Dimensionless DoubleMatrix gt(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return gti(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for "greater than" against a scalar (in-place). */
    public @Dimensionless DoubleMatrix gti(@Dimensionless DoubleMatrix this, @Dimensionless double value, @Dimensionless DoubleMatrix result) {
      ensureResultLength(null, result);
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) > value ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for "greater than" against a scalar (in-place). */
    public @Dimensionless DoubleMatrix gti(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return gti(value, this);
    }

    /** test for "greater than" against a scalar. */
    public @Dimensionless DoubleMatrix gt(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return gti(value, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for "less than or equal" (in-place). */
    public @Dimensionless DoubleMatrix lei(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
      if (other.isScalar())
        return lei(other.scalar(), result);

      assertSameLength(other);
      ensureResultLength(other, result);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) <= other.get(i) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for "less than or equal" (in-place). */
    public @Dimensionless DoubleMatrix lei(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return lei(other, this);
    }

    /** Test for "less than or equal". */
    public @Dimensionless DoubleMatrix le(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return lei(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for "less than or equal" against a scalar (in-place). */
    public @Dimensionless DoubleMatrix lei(@Dimensionless DoubleMatrix this, @Dimensionless double value, @Dimensionless DoubleMatrix result) {
      ensureResultLength(null, result);
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) <= value ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for "less than or equal" against a scalar (in-place). */
    public @Dimensionless DoubleMatrix lei(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return lei(value, this);
    }

    /** test for "less than or equal" against a scalar. */
    public @Dimensionless DoubleMatrix le(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return lei(value, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for "greater than or equal" (in-place). */
    public @Dimensionless DoubleMatrix gei(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
      if (other.isScalar())
        return gei(other.scalar(), result);

      assertSameLength(other);
      ensureResultLength(other, result);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) >= other.get(i) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for "greater than or equal" (in-place). */
    public @Dimensionless DoubleMatrix gei(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return gei(other, this);
    }

    /** Test for "greater than or equal". */
    public @Dimensionless DoubleMatrix ge(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return gei(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for "greater than or equal" against a scalar (in-place). */
    public @Dimensionless DoubleMatrix gei(@Dimensionless DoubleMatrix this, @Dimensionless double value, @Dimensionless DoubleMatrix result) {
      ensureResultLength(null, result);
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) >= value ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for "greater than or equal" against a scalar (in-place). */
    public @Dimensionless DoubleMatrix gei(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return gei(value, this);
    }

    /** test for "greater than or equal" against a scalar. */
    public @Dimensionless DoubleMatrix ge(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return gei(value, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for equality (in-place). */
    public @Dimensionless DoubleMatrix eqi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
      if (other.isScalar())
        return eqi(other.scalar(), result);

      assertSameLength(other);
      ensureResultLength(other, result);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) == other.get(i) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for equality (in-place). */
    public @Dimensionless DoubleMatrix eqi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return eqi(other, this);
    }

    /** Test for equality. */
    public @Dimensionless DoubleMatrix eq(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return eqi(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for equality against a scalar (in-place). */
    public @Dimensionless DoubleMatrix eqi(@Dimensionless DoubleMatrix this, @Dimensionless double value, @Dimensionless DoubleMatrix result) {
      ensureResultLength(null, result);
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) == value ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for equality against a scalar (in-place). */
    public @Dimensionless DoubleMatrix eqi(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return eqi(value, this);
    }

    /** test for equality against a scalar. */
    public @Dimensionless DoubleMatrix eq(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return eqi(value, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for inequality (in-place). */
    public @Dimensionless DoubleMatrix nei(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
      if (other.isScalar())
        return nei(other.scalar(), result);

      assertSameLength(other);
      ensureResultLength(other, result);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) != other.get(i) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for inequality (in-place). */
    public @Dimensionless DoubleMatrix nei(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return nei(other, this);
    }

    /** Test for inequality. */
    public @Dimensionless DoubleMatrix ne(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return nei(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Test for inequality against a scalar (in-place). */
    public @Dimensionless DoubleMatrix nei(@Dimensionless DoubleMatrix this, @Dimensionless double value, @Dimensionless DoubleMatrix result) {
      ensureResultLength(null, result);
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, get(i) != value ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Test for inequality against a scalar (in-place). */
    public @Dimensionless DoubleMatrix nei(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return nei(value, this);
    }

    /** test for inequality against a scalar. */
    public @Dimensionless DoubleMatrix ne(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return nei(value, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Compute elementwise logical and (in-place). */
    public @Dimensionless DoubleMatrix andi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
      assertSameLength(other);
      ensureResultLength(other, result);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, (get(i) != ((@Dimensionless double) (0.0))) & (other.get(i) != ((@Dimensionless double) (0.0))) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Compute elementwise logical and (in-place). */
    public @Dimensionless DoubleMatrix andi(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return andi(other, this);
    }

    /** Compute elementwise logical and. */
    public @Dimensionless DoubleMatrix and(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return andi(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Compute elementwise logical and against a scalar (in-place). */
    public @Dimensionless DoubleMatrix andi(@Dimensionless DoubleMatrix this, @Dimensionless double value, @Dimensionless DoubleMatrix result) {
      ensureResultLength(null, result);
      @Dimensionless
      boolean val = (value != ((@Dimensionless double) (0.0)));
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, (get(i) != ((@Dimensionless double) (0.0))) & val ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Compute elementwise logical and against a scalar (in-place). */
    public @Dimensionless DoubleMatrix andi(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return andi(value, this);
    }

    /** Compute elementwise logical and against a scalar. */
    public @Dimensionless DoubleMatrix and(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return andi(value, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Compute elementwise logical or (in-place). */
    public @Dimensionless DoubleMatrix ori(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
      assertSameLength(other);
      ensureResultLength(other, result);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, (get(i) != ((@Dimensionless double) (0.0))) | (other.get(i) != ((@Dimensionless double) (0.0))) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Compute elementwise logical or (in-place). */
    public @Dimensionless DoubleMatrix ori(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return ori(other, this);
    }

    /** Compute elementwise logical or. */
    public @Dimensionless DoubleMatrix or(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return ori(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Compute elementwise logical or against a scalar (in-place). */
    public @Dimensionless DoubleMatrix ori(@Dimensionless DoubleMatrix this, @Dimensionless double value, @Dimensionless DoubleMatrix result) {
      ensureResultLength(null, result);
      @Dimensionless
      boolean val = (value != ((@Dimensionless double) (0.0)));
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, (get(i) != ((@Dimensionless double) (0.0))) | val ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Compute elementwise logical or against a scalar (in-place). */
    public @Dimensionless DoubleMatrix ori(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return ori(value, this);
    }

    /** Compute elementwise logical or against a scalar. */
    public @Dimensionless DoubleMatrix or(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return ori(value, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Compute elementwise logical xor (in-place). */
    public @Dimensionless DoubleMatrix xori(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other, @Dimensionless DoubleMatrix result) {
      assertSameLength(other);
      ensureResultLength(other, result);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, (get(i) != ((@Dimensionless double) (0.0))) ^ (other.get(i) != ((@Dimensionless double) (0.0))) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Compute elementwise logical xor (in-place). */
    public @Dimensionless DoubleMatrix xori(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return xori(other, this);
    }

    /** Compute elementwise logical xor. */
    public @Dimensionless DoubleMatrix xor(@Dimensionless DoubleMatrix this, @Dimensionless DoubleMatrix other) {
      return xori(other, new @Dimensionless DoubleMatrix(rows, columns));
    }

    /** Compute elementwise logical xor against a scalar (in-place). */
    public @Dimensionless DoubleMatrix xori(@Dimensionless DoubleMatrix this, @Dimensionless double value, @Dimensionless DoubleMatrix result) {
      ensureResultLength(null, result);
      @Dimensionless
      boolean val = (value != ((@Dimensionless double) (0.0)));
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
        result.put(i, (get(i) != ((@Dimensionless double) (0.0))) ^ val ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
      return result;
    }

    /** Compute elementwise logical xor against a scalar (in-place). */
    public @Dimensionless DoubleMatrix xori(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return xori(value, this);
    }

    /** Compute elementwise logical xor against a scalar. */
    public @Dimensionless DoubleMatrix xor(@Dimensionless DoubleMatrix this, @Dimensionless double value) {
      return xori(value, new @Dimensionless DoubleMatrix(rows, columns));
    }
//RJPP-END--------------------------------------------------------------

    public @Dimensionless ComplexDoubleMatrix toComplex(@Dimensionless DoubleMatrix this) {
      return new @Dimensionless ComplexDoubleMatrix(this);
    }
}
