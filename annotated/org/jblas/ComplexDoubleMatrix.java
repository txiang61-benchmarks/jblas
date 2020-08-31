// --- BEGIN LICENSE BLOCK ---
/* 
 * Copyright (c) 2009, Mikio L. Braun
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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

public class ComplexDoubleMatrix {
	
	public @Dimensionless int rows;
	public @Dimensionless int columns;
	public @Dimensionless int length;
	public @Dimensionless double @Dimensionless [] data = null; // rows are contiguous

	/**************************************************************************
	 * 
	 * Constructors and factory functions
	 * 
	 **************************************************************************/

	/**
   * Create a new matrix with <i>newRows</i> rows, <i>newColumns</i> columns
	 * using <i>newData></i> as the data.
	 */
	public ComplexDoubleMatrix(@Dimensionless int newRows, @Dimensionless int newColumns, @Dimensionless double @Dimensionless ... newData) {
		rows = newRows;
		columns = newColumns;
		length = rows * columns;

    if (newData.length != ((@Dimensionless int) (2)) * newRows * newColumns)
			throw new @Dimensionless IllegalArgumentException(
					"Passed data must match matrix dimensions.");
                data = newData;
	}
	
	/**
	 * Creates a new <i>n</i> times <i>m</i> <tt>ComplexDoubleMatrix</tt>.
	 * @param newRows the number of rows (<i>n</i>) of the new matrix.
	 * @param newColumns the number of columns (<i>m</i>) of the new matrix.
	 */
	public ComplexDoubleMatrix(@Dimensionless int newRows, @Dimensionless int newColumns) {
		this(newRows, newColumns, new @Dimensionless double @Dimensionless [((@Dimensionless int) (2)) * newRows * newColumns]);
	}
	
	/**
	 * Creates a new <tt>ComplexDoubleMatrix</tt> of size 0 times 0.
	 */
	public ComplexDoubleMatrix() {
		this(((@Dimensionless int) (0)), ((@Dimensionless int) (0)), null);
	}

	/**
	 * Create a Matrix of length <tt>len</tt>. By default, this creates a row vector.
	 * @param len
	 */
	public ComplexDoubleMatrix(@Dimensionless int len) {
		this(len, ((@Dimensionless int) (1)), new @Dimensionless double @Dimensionless [((@Dimensionless int) (2)) * len]);
	}
	
	public ComplexDoubleMatrix(@Dimensionless double @Dimensionless [] newData) {
		this(newData.length/ ((@Dimensionless int) (2)), ((@Dimensionless int) (1)), newData);
	}

	public ComplexDoubleMatrix(@Dimensionless ComplexDouble @Dimensionless [] newData) {
		this(newData.length);
				
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < newData.length; i++)
			put(i, newData[i]);
	}
		
        
  /** Construct a complex matrix from a real matrix. */
  public ComplexDoubleMatrix(DoubleMatrix m) {
    this(m.rows, m.columns);

    NativeBlas.dcopy(m.length, m.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, ((@Dimensionless int) (0)), ((@Dimensionless int) (2)));
  }

  /** Construct a complex matrix from separate real and imaginary parts. Either
   * part can be set to null in which case it will be ignored.
   */
  public ComplexDoubleMatrix(@Dimensionless DoubleMatrix real, @Dimensionless DoubleMatrix imag) {
      this((real != null) ? real.rows : imag.rows, (real != null) ? real.columns : imag.columns);

      if (real != null && imag != null)
      real.assertSameSize(imag);

      if (real != null)
          NativeBlas.dcopy(length, real.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, ((@Dimensionless int) (0)), ((@Dimensionless int) (2)));
      if (imag != null)
          NativeBlas.dcopy(length, imag.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, ((@Dimensionless int) (1)), ((@Dimensionless int) (2)));
  }
        
        /**
	 * Creates a new matrix by reading it from a file.
	 * @param filename the path and name of the file to read the matrix from
	 * @throws IOException 
	 */
	public ComplexDoubleMatrix(@Dimensionless String filename) throws IOException {
		load(filename);
	}
	
	/**
	 * Creates a new <i>n</i> times <i>m</i> <tt>ComplexDoubleMatrix</tt> from
	 * the given <i>n</i> times <i>m</i> 2D data array. The first dimension of the array makes the
	 * rows (<i>n</i>) and the second dimension the columns (<i>m</i>). For example, the
	 * given code <br/><br/>
	 * <code>new ComplexDoubleMatrix(new double[][]{{1d, 2d, 3d}, {4d, 5d, 6d}, {7d, 8d, 9d}}).print();</code><br/><br/>
	 * will constructs the following matrix:
	 * <pre>
	 * 1.0	2.0	3.0
	 * 4.0	5.0	6.0
	 * 7.0	8.0	9.0
	 * </pre>.
	 * @param data <i>n</i> times <i>m</i> data array
	 */ 
	public ComplexDoubleMatrix(@Dimensionless double @Dimensionless [] @Dimensionless [] data) {
		this(data.length, data[((@Dimensionless int) (0))].length);
						
		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++)
			assert(data[r].length == columns);
		
		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++)
			for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++)
				put(r, c, data[r][c]);
	}
	
	/**
	 * Creates a new matrix in which all values are equal 0.
	 * @param rows number of rows
	 * @param columns number of columns
	 * @return new matrix
	 */
	public static @Dimensionless ComplexDoubleMatrix zeros(@Dimensionless int rows, @Dimensionless int columns) {
		return new @Dimensionless ComplexDoubleMatrix(rows, columns);
	}
	
	public static @Dimensionless ComplexDoubleMatrix zeros(@Dimensionless int length) {
		return zeros(length, ((@Dimensionless int) (1)));
	}

	/**
	 * Creates a new matrix in which all values are equal 1.
	 * @param rows number of rows
	 * @param columns number of columns
	 * @return new matrix
	 */
	public static @Dimensionless ComplexDoubleMatrix ones(@Dimensionless int rows, @Dimensionless int columns) {
		@Dimensionless
		ComplexDoubleMatrix m = new @Dimensionless ComplexDoubleMatrix(rows, columns);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows * columns; i++)
			m.put(i, ((@Dimensionless double) (1.0)));
		
		return m;
	}
	
	public static @Dimensionless ComplexDoubleMatrix ones(@Dimensionless int length) {
		return ones(length, ((@Dimensionless int) (1)));
	}
	
	/**
	 * Creates a new matrix where the values of the given vector are the diagonal values of
	 * the matrix.
	 * @param x the diagonal values
	 * @return new matrix
	 */
	public static @Dimensionless ComplexDoubleMatrix diag(@Dimensionless ComplexDoubleMatrix x) {
		@Dimensionless
		ComplexDoubleMatrix m = new @Dimensionless ComplexDoubleMatrix(x.length, x.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
			m.put(i, i, x.get(i));
		
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
  public static @Dimensionless ComplexDoubleMatrix diag(@Dimensionless ComplexDoubleMatrix x, @Dimensionless int rows, @Dimensionless int columns) {
    if (x.length > rows || x.length > columns) {
      throw new @Dimensionless SizeException("Length of diagonal matrix must be larger than both rows and columns.");
    }
    
    @Dimensionless
    ComplexDoubleMatrix m = new @Dimensionless ComplexDoubleMatrix(rows, columns);

    for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
      m.put(i, i, x.get(i));

    return m;
  }
	
	/**
	 * Create a 1 * 1 - matrix. For many operations, this matrix functions like a
	 * normal double
	 * @param s value of the matrix
	 * @return the constructed ComplexDoubleMatrix 
	 */
	public static @Dimensionless ComplexDoubleMatrix scalar(@Dimensionless double s) {
		@Dimensionless
		ComplexDoubleMatrix m = new @Dimensionless ComplexDoubleMatrix(((@Dimensionless int) (1)), ((@Dimensionless int) (1)));
		m.put(((@Dimensionless int) (0)), ((@Dimensionless int) (0)), s);
		return m;
	}
	
	/** Test whether a matrix is scalar */
	public @Dimensionless boolean isScalar(@Dimensionless ComplexDoubleMatrix this) {
		return length == ((@Dimensionless int) (1));
	}
	
	/** Return the first element of the matrix */
	public @Dimensionless ComplexDouble scalar(@Dimensionless ComplexDoubleMatrix this) {
		return get(((@Dimensionless int) (0)));
	}
	
	public static @Dimensionless ComplexDoubleMatrix concatHorizontally(@Dimensionless ComplexDoubleMatrix A, @Dimensionless ComplexDoubleMatrix B) {
		if (A.rows != B.rows)
			throw new @Dimensionless SizeException("Matrices don't have same number of rows.");
		
		@Dimensionless
		ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(A.rows, A.columns + B.columns);
		SimpleBlas.copy(A, result);
		NativeBlas.zcopy(B.length, B.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), result.data, A.length, ((@Dimensionless int) (1)));
		return result;
	}

	public static @Dimensionless ComplexDoubleMatrix concatVertically(@Dimensionless ComplexDoubleMatrix A, @Dimensionless ComplexDoubleMatrix B) {
		if (A.columns != B.columns)
			throw new @Dimensionless SizeException("Matrices don't have same number of columns.");
		
		@Dimensionless
		ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(A.rows + B.rows, A.columns);

		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < A.columns; i++) {
			NativeBlas.zcopy(A.rows, A.data, A.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)), result.data, result.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)));
			NativeBlas.zcopy(B.rows, B.data, B.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)), result.data, result.index(A.rows, i), ((@Dimensionless int) (1)));
		}
		
		return result;
	}
	
	/**************************************************************************
	 * Working with slices (Man! 30+ methods just to make this a bit flexible...) 
	 */

	public @Dimensionless ComplexDoubleMatrix get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices) {
		@Dimensionless
		ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			result.put(i, get(indices[i]));
		
		return result;
	}
	
	public @Dimensionless ComplexDoubleMatrix get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices) {
		@Dimensionless
		ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(((@Dimensionless int) (1)), indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			result.put(i, get(r, indices[i]));
		
		return result;
	}
	
	public @Dimensionless ComplexDoubleMatrix get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c) {
		@Dimensionless
		ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(indices.length, ((@Dimensionless int) (1)));
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			result.put(i, get(indices[i], c));
		
		return result;
	}
	
	public @Dimensionless ComplexDoubleMatrix get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices) {
		@Dimensionless
		ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(rindices.length, cindices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				result.put(i, j, get(rindices[i], cindices[j]));
		
		return result;
	}
	
	public @Dimensionless ComplexDoubleMatrix get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices) {
		return get(indices.findIndices());
	}

	public @Dimensionless ComplexDoubleMatrix get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless ComplexDoubleMatrix indices) {
		return get(r, indices.findIndices());
	}
	
	public @Dimensionless ComplexDoubleMatrix get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless int c) {
		return get(indices.findIndices(), c);
	}

	public @Dimensionless ComplexDoubleMatrix get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix rindices, @Dimensionless ComplexDoubleMatrix cindices) {
		return get(rindices.findIndices(), cindices.findIndices());
	}
	
	private void checkLength(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int l) {
		if (length != l)
			throw new @Dimensionless SizeException("Matrix does not have the necessary length (" + length + " != " + l + ").");
	}

	private void checkRows(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r) {
		if (rows != r)
			throw new @Dimensionless SizeException("Matrix does not have the necessary length (" + length + " != " + r + ").");
	}
	
	private void checkColumns(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int c) {
		if (columns != c)
			throw new @Dimensionless SizeException("Matrix does not have the necessary length (" + length + " != " + c + ").");
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless ComplexDoubleMatrix x) {
		if (x.isScalar())
			return put(indices, x.scalar());
		x.checkLength(indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], x.get(i));
		
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless ComplexDoubleMatrix x) {
		if (x.isScalar())
			return put(r, indices, x.scalar());
		x.checkColumns(indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(r, indices[i], x.get(i));
		
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless ComplexDoubleMatrix x) {
		if (x.isScalar())
			return put(indices, c, x.scalar());		
		x.checkRows(indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], c, x.get(i));
		
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless ComplexDoubleMatrix x) {
		if (x.isScalar())
			return put(rindices, cindices, x.scalar());		
		x.checkRows(rindices.length);
		x.checkColumns(cindices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				put(rindices[i], cindices[j], x.get(i,j));
		
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless double v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless double v) {
		return put(indices, v);
	}

	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless double v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			putImag(indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless ComplexDouble v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless double v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(r, indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless double v) {
		return put(r, indices, v);
	}

	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless double v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			putImag(r, indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless ComplexDouble v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(r, indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless double v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], c, v);
		
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless double v) {
		return put(indices, c, v);
	}
	
	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless double v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			putImag(indices[i], c, v);
		
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless ComplexDouble v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], c, v);
		
		return this;
 	}
	
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless double v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				put(rindices[i], cindices[j], v);
		
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless double v) {
		return put(rindices, cindices, v);
	}
	
	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless double v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				put(rindices[i], cindices[j], v);
		
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless ComplexDouble v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				put(rindices[i], cindices[j], v);
		
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless ComplexDoubleMatrix v) {
		return put(indices.findIndices(), v);
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless ComplexDoubleMatrix v) {
		return put(r, indices.findIndices(), v);
	}
	
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless int c, @Dimensionless ComplexDoubleMatrix v) {
		return put(indices.findIndices(), c, v);
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix rindices, @Dimensionless ComplexDoubleMatrix cindices, @Dimensionless ComplexDoubleMatrix v) {
		return put(rindices.findIndices(), cindices.findIndices(), v);
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless double v) {
		return put(indices.findIndices(), v);
	}

	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless double v) {
		return put(indices, v);
	}

	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless double v) {
		return putImag(indices.findIndices(), v);
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless ComplexDouble v) {
		return put(indices.findIndices(), v);
	}
	
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless double v) {
		return put(r, indices.findIndices(), v);
	}
	
	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless double v) {
		return put(r, indices, v);
	}

	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless double v) {
		return putImag(r, indices.findIndices(), v);
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless ComplexDouble v) {
		return put(r, indices.findIndices(), v);
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless int c, @Dimensionless double v) {
		return put(indices.findIndices(), c, v);
	}

	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless int c, @Dimensionless double v) {
		return put(indices, c, v);
	}

	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless int c, @Dimensionless double v) {
		return putImag(indices.findIndices(), c, v);
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix indices, @Dimensionless int c, @Dimensionless ComplexDouble v) {
		return put(indices.findIndices(), c, v);
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix rindices, @Dimensionless ComplexDoubleMatrix cindices, @Dimensionless double v) {
		return put(rindices.findIndices(), cindices.findIndices(), v);
	}

	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix rindices, @Dimensionless ComplexDoubleMatrix cindices, @Dimensionless double v) {
		return putReal(rindices.findIndices(), cindices.findIndices(), v);
	}

	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix rindices, @Dimensionless ComplexDoubleMatrix cindices, @Dimensionless double v) {
		return putImag(rindices.findIndices(), cindices.findIndices(), v);
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix rindices, @Dimensionless ComplexDoubleMatrix cindices, @Dimensionless ComplexDouble v) {
		return put(rindices.findIndices(), cindices.findIndices(), v);
	}

	
	public @Dimensionless int @Dimensionless [] findIndices(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		int len = ((@Dimensionless int) (0));
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			if (!get(i).isZero())
				len++;
		
		@Dimensionless
		int @Dimensionless [] indices = new @Dimensionless int @Dimensionless [len];
		@Dimensionless
		int c = ((@Dimensionless int) (0));
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			if (!get(i).isZero())
				indices[c++] = i;
		
		return indices;
	}
	
	/**************************************************************************
	 * Basic operations (copying, resizing, element access)
	 */
	
	/** Return transposed copy of this matrix */
	public @Dimensionless ComplexDoubleMatrix transpose(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(columns, rows);

                @Dimensionless
                ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless int) (0)));

		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < columns; j++)
				result.put(j, i, get(i, j, c));
		
		return result;
	}

        public @Dimensionless ComplexDoubleMatrix hermitian(@Dimensionless ComplexDoubleMatrix this) {
            @Dimensionless
            ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(columns, rows);

            @Dimensionless
            ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless int) (0)));

            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows; i++)
                for (@Dimensionless int j = ((@Dimensionless int) (0)); j < columns; j++)
                    result.put(j, i, get(i, j, c).conji());
            return result;
        }

        /**
         * Compute complex conjugate (in-place).
         */
        public @Dimensionless ComplexDoubleMatrix conji(@Dimensionless ComplexDoubleMatrix this) {
            @Dimensionless
            ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                put(i, get(i, c).conji());
            return this;
        }

        /**
         * Compute complex conjugate.
         */
        public @Dimensionless ComplexDoubleMatrix conj(@Dimensionless ComplexDoubleMatrix this) {
            return dup().conji();
        }

		
	/** Compare two matrices.
	 * @param o Object to compare to
	 * @return true if and only if other is also a ComplexDoubleMatrix which has the same size and the
	 * maximal absolute difference in matrix elements is smaller thatn 1e-6.  */
	public @Dimensionless boolean equals(@Dimensionless ComplexDoubleMatrix this, @Dimensionless Object o) {
		if (!(o instanceof ComplexDoubleMatrix))
			return false;

		@Dimensionless
		ComplexDoubleMatrix other = (@Dimensionless ComplexDoubleMatrix) o;

		if (!sameSize(other))
			return false;

    return Arrays.equals(data, other.data);
	}
  
  public @Dimensionless int hashCode(@Dimensionless ComplexDoubleMatrix this) {
    return rows ^ columns ^ Arrays.hashCode(data);
  }

	
	/** Resize the matrix. All elements will be set to zero. */
	public void resize(@Dimensionless int newRows, @Dimensionless int newColumns) {
		rows = newRows;
		columns = newColumns;
		length = newRows * newColumns;
		data = new @Dimensionless double @Dimensionless [((@Dimensionless int) (2)) * rows * columns];
	}

	
	/** Reshape the matrix. Number of elements must not change. */
	public @Dimensionless ComplexDoubleMatrix reshape(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int newRows, @Dimensionless int newColumns) {
		if (length != newRows * newColumns)
			throw new @Dimensionless IllegalArgumentException(
					"Number of elements must not change.");

		rows = newRows;
		columns = newColumns;
		
		return this;
	}

	/** Checks whether two matrices have the same size. */
	public @Dimensionless boolean sameSize(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix a) {
		return rows == a.rows && columns == a.columns;
	}

	/** 
	 * Assert that two matrices have the same size.
	 * 
	 * @param a the other matrix
	 * @throws SizeException if matrix sizes don't match. 
	 * */
	public void assertSameSize(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix a) {
		if (!sameSize(a))
			throw new @Dimensionless SizeException("Matrices must have the same size.");
	}
	
	/** 
	 * Check whether this can be multiplied with a. 
	 * 
	 * @param a right-hand-side of the multiplication.
	 * @return true iff <tt>this.columns == a.rows</tt>
	 */
	public @Dimensionless boolean multipliesWith(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix a) {
		return columns == a.rows;
	}
	
	public void assertMultipliesWith(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix a) {
		if (!multipliesWith(a))
			throw new @Dimensionless SizeException("Number of columns of left matrix must be equal to number of rows of right matrix.");
	}
	
	public @Dimensionless boolean sameLength(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix a) {
		return length == a.length;
	}
	
	public void assertSameLength(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix a) {
		if (!sameLength(a))
			throw new @Dimensionless SizeException("Matrices must have same length (is: " + length + " and " + a.length + ")");
	}
	
	/** Copy ComplexDoubleMatrix a to this. this a is resized if necessary. */
	public @Dimensionless ComplexDoubleMatrix copy(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix a) {
		if (!sameSize(a))
			resize(a.rows, a.columns);
		
		SimpleBlas.copy(a, this);
		return a;
	}
	
	/** Returns a duplicate of this matrix. Geometry is the same (including offsets, transpose, etc.),
	 * but the buffer is not shared.
	 */
	public @Dimensionless ComplexDoubleMatrix dup(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		ComplexDoubleMatrix out = new @Dimensionless ComplexDoubleMatrix(rows, columns);

                JavaBlas.rcopy(((@Dimensionless int) (2))*length, data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), out.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		
		return out;
	}
	
	public @Dimensionless ComplexDoubleMatrix swapColumns(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i, @Dimensionless int j) {
		NativeBlas.zswap(rows, data, index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), j), ((@Dimensionless int) (1)));
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix swapRows(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i, @Dimensionless int j) {
		NativeBlas.zswap(columns, data, index(i, ((@Dimensionless int) (0))), rows, data, index(j, ((@Dimensionless int) (0))), rows);
		return this;
	}
		
	/** Set matrix element */
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless double value) {
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)] =  value;
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless double realValue, @Dimensionless double complexValue) {
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)] =  realValue;
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)+ ((@Dimensionless int) (1))] =  complexValue;
		return this;
	}

        public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless ComplexDouble value) {
		@Dimensionless
		int i = ((@Dimensionless int) (2))*index(rowIndex, columnIndex);
		data[i] = value.real(); data[i+ ((@Dimensionless int) (1))] = value.imag();
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless double value) {
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)] = value;
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless double value) {
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)+ ((@Dimensionless int) (1))] = value;
		return this;
	}
	
	/** Retrieve matrix element */
	public @Dimensionless ComplexDouble get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex) {
            @Dimensionless
            int i = ((@Dimensionless int) (2))*index(rowIndex, columnIndex);
            return new @Dimensionless ComplexDouble(data[i], data[i+ ((@Dimensionless int) (1))]);
	}

        /** Get matrix element, passing the variable to store the result. */
        public @Dimensionless ComplexDouble get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless ComplexDouble result) {
            return get(index(rowIndex, columnIndex), result);
        }
	
	public @Dimensionless DoubleMatrix getReal(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		DoubleMatrix result = new @Dimensionless DoubleMatrix(rows, columns);
		
		NativeBlas.dcopy(length, data, ((@Dimensionless int) (0)), ((@Dimensionless int) (2)), result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		
		return result;
	}

	/** Get index of an element */
	public @Dimensionless int index(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex) {
		return rows * columnIndex + rowIndex;
	}

  /** Compute the row index of a linear index. */
  public @Dimensionless int indexRows(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i) {
    return i - indexColumns(i) * rows;
  }

  /** Compute the column index of a linear index. */
  public @Dimensionless int indexColumns(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i) {
    return i / rows;
  }


	public @Dimensionless ComplexDouble get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i) {
		return new @Dimensionless ComplexDouble(data[i * ((@Dimensionless int) (2))], data[i * ((@Dimensionless int) (2)) + ((@Dimensionless int) (1))]);
	}
	
        public @Dimensionless ComplexDouble get(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i, @Dimensionless ComplexDouble result) {
            return result.set(data[i * ((@Dimensionless int) (2))], data[i* ((@Dimensionless int) (2))+ ((@Dimensionless int) (1))]);
        }
        
	public @Dimensionless double getReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i) {
		return data[((@Dimensionless int) (2))*i];
	}
	
	public @Dimensionless double getImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i) {
		return data[((@Dimensionless int) (2))*i + ((@Dimensionless int) (1))]; 
	}

	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i, @Dimensionless double v) {
		data[((@Dimensionless int) (2))*i] = v;
		return this;
	}

        public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i, @Dimensionless double r, @Dimensionless double c) {
            data[((@Dimensionless int) (2))*i] = r;
            data[((@Dimensionless int) (2))*i+ ((@Dimensionless int) (1))] = c;
            return this;
        }
	
	public @Dimensionless ComplexDoubleMatrix put(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i, @Dimensionless ComplexDouble v) {
		data[((@Dimensionless int) (2))*i] = v.real();
		data[((@Dimensionless int) (2))*i+ ((@Dimensionless int) (1))] = v.imag();
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix putReal(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i, @Dimensionless double v) {
		return put(i, v);
	}
	
	public @Dimensionless ComplexDoubleMatrix putImag(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int i, @Dimensionless double v) {
		data[((@Dimensionless int) (2))*i+ ((@Dimensionless int) (1))] = v;
		return this;
	}

	public @Dimensionless int getRows(@Dimensionless ComplexDoubleMatrix this) {
		return rows;
	}
	
	public @Dimensionless int getColumns(@Dimensionless ComplexDoubleMatrix this) {
		return columns;
	}
	
	public @Dimensionless int getLength(@Dimensionless ComplexDoubleMatrix this) {
		return length;
	}
	
	/** Checks whether the matrix is empty. */
	public @Dimensionless boolean isEmpty(@Dimensionless ComplexDoubleMatrix this) {
		return columns == ((@Dimensionless int) (0)) || rows == ((@Dimensionless int) (0));
	}
	
	/** Checks whether the matrix is square. */
	public @Dimensionless boolean isSquare(@Dimensionless ComplexDoubleMatrix this) {
		return columns == rows;
	}
	
	public void assertSquare(@Dimensionless ComplexDoubleMatrix this) {
		if (!isSquare())
			throw new @Dimensionless SizeException("Matrix must be square!");
	}
	
	/** Checks whether the matrix is a vector. */
	public @Dimensionless boolean isVector(@Dimensionless ComplexDoubleMatrix this) {
		return columns == ((@Dimensionless int) (1)) || rows == ((@Dimensionless int) (1));
	}
	
	public @Dimensionless boolean isRowVector(@Dimensionless ComplexDoubleMatrix this) {
		return columns == ((@Dimensionless int) (1));
	}
	
	public @Dimensionless boolean isColumnVector(@Dimensionless ComplexDoubleMatrix this) {
		return rows == ((@Dimensionless int) (1));
	}
		
        /** Get diagonal of the matrix. */
	public @Dimensionless ComplexDoubleMatrix diag(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		ComplexDoubleMatrix d = new @Dimensionless ComplexDoubleMatrix(rows);
		NativeBlas.zcopy(rows, data, ((@Dimensionless int) (0)), rows + ((@Dimensionless int) (1)), d.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		return d;
	}
        
        /** Get real part of the matrix. */
        public @Dimensionless DoubleMatrix real(@Dimensionless ComplexDoubleMatrix this) {
            @Dimensionless
            DoubleMatrix result = new @Dimensionless DoubleMatrix(rows, columns);
            NativeBlas.dcopy(length, data, ((@Dimensionless int) (0)), ((@Dimensionless int) (2)), result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
            return result;
        }
        
        /** Get imaginary part of the matrix. */
        public @Dimensionless DoubleMatrix imag(@Dimensionless ComplexDoubleMatrix this) {
            @Dimensionless
            DoubleMatrix result = new @Dimensionless DoubleMatrix(rows, columns);
            NativeBlas.dcopy(length, data, ((@Dimensionless int) (1)), ((@Dimensionless int) (2)), result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
            return result;            
        }

	
	/** 
	 * Pretty-print this matrix to <tt>System.out</tt>. 
	 * */
	public void print(@Dimensionless ComplexDoubleMatrix this) {
		System.out.println(toString());
	}

	/** 
	 * Generate string representation of this matrix 
	 * (multi-line).
	 * */
	public @Dimensionless String toString(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		StringBuilder s = new @Dimensionless StringBuilder();

		s.append("[");
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows; i++) {
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < columns; j++) {
				s.append(get(i, j));
				if (j < columns - ((@Dimensionless int) (1)))
					s.append(", ");
			}
			if (i < rows - ((@Dimensionless int) (1)))
				s.append("; ");
		}

		s.append("]");
		
		return s.toString();
	}

	public @Dimensionless double @Dimensionless [] toDoubleArray(@Dimensionless ComplexDoubleMatrix this) {
		return data.clone();
	}
	
	public @Dimensionless ComplexDouble @Dimensionless [] toArray(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		ComplexDouble @Dimensionless [] array = new ComplexDouble @Dimensionless [length];
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			array[i] = get(i);
		
		return array;		
	}
	
	public @Dimensionless ComplexDouble @Dimensionless [] @Dimensionless [] toArray2(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		ComplexDouble @Dimensionless [] @Dimensionless [] array = new ComplexDouble @Dimensionless [rows][columns];
		
		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++)
			for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++)
				array[r][c] = get(r, c);
				
		return array;
	}
	
	public @Dimensionless boolean @Dimensionless [] toBooleanArray(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		boolean @Dimensionless [] array = new @Dimensionless boolean @Dimensionless [length];
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			array[i] = !get(i).isZero();
		
		return array;
	}
	
	public @Dimensionless boolean @Dimensionless [] @Dimensionless [] toBooleanArray2(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		boolean @Dimensionless [] @Dimensionless [] array = new @Dimensionless boolean @Dimensionless [rows] @Dimensionless [columns];
		
		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++)
			for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++)
				array[r][c] = !get(r, c).isZero();
				
		return array;
	}

	/**************************************************************************
	 * Arithmetic Operations
	 */

	/** 
	 * Ensures that the result vector has the same length as this. If not,
	 * resizing result is tried, which fails if result == this or result == other.
	 */
	private void ensureResultLength(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
		if (!sameLength(result)) {
			if (result == this || result == other)
				throw new @Dimensionless SizeException("Cannot resize result matrix because it is used in-place.");
			result.resize(rows, columns);
		}
	}

	/** Add two matrices. */
	public @Dimensionless ComplexDoubleMatrix addi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
		if (other.isScalar())
			return addi(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);
		
		if (result == this)
			SimpleBlas.axpy(ComplexDouble.UNIT, other, result);
		else if (result == other)
			SimpleBlas.axpy(ComplexDouble.UNIT, this, result);
		else {
			SimpleBlas.copy(this, result);
			SimpleBlas.axpy(ComplexDouble.UNIT, other, result);
		}

		return result;
	}
	
	/** Add a scalar to a matrix. */
	public @Dimensionless ComplexDoubleMatrix addi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v, @Dimensionless ComplexDoubleMatrix result) {
		ensureResultLength(null, result);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i).add(v));
		return result;
	}
	
	public @Dimensionless ComplexDoubleMatrix addi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v, @Dimensionless ComplexDoubleMatrix result) {
		return addi(new @Dimensionless ComplexDouble(v), result);
	}

	/** Subtract two matrices. */
	public @Dimensionless ComplexDoubleMatrix subi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
		if (other.isScalar())
			return subi(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);
		
		if (result == this)
			SimpleBlas.axpy(ComplexDouble.NEG_UNIT, other, result);
		else if (result == other) {
			SimpleBlas.scal(ComplexDouble.NEG_UNIT, result);
			SimpleBlas.axpy(ComplexDouble.UNIT, this, result);
		}
		else {
			SimpleBlas.copy(this, result);
			SimpleBlas.axpy(ComplexDouble.NEG_UNIT, other, result);
		}
		return result;
	}
	
	/** Subtract a scalar from a matrix */
	public @Dimensionless ComplexDoubleMatrix subi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v, @Dimensionless ComplexDoubleMatrix result) {
		ensureResultLength(null, result);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i).sub(v));
		return result;
	}
	
	public @Dimensionless ComplexDoubleMatrix subi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v, @Dimensionless ComplexDoubleMatrix result) {
		return subi(new @Dimensionless ComplexDouble(v), result);
	}

	/** 
	 * Subtract two matrices, but subtract first from second matrix, that is, 
	 * compute <em>result = other - this</em>. 
	 * */
	public @Dimensionless ComplexDoubleMatrix rsubi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
		return other.subi(this, result);
	}
	
	/** Subtract a matrix from a scalar */
	public @Dimensionless ComplexDoubleMatrix rsubi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble a, @Dimensionless ComplexDoubleMatrix result) {
		ensureResultLength(null, result);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, a.sub(get(i)));
		return result;
	}

	public @Dimensionless ComplexDoubleMatrix rsubi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double a, @Dimensionless ComplexDoubleMatrix result) {
		return rsubi(new @Dimensionless ComplexDouble(a), result);
	}

	/** (Elementwise) Multiplication */ 
	public @Dimensionless ComplexDoubleMatrix muli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
		if (other.isScalar())
			return muli(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);
		
                @Dimensionless
                ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble d = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i, c).muli(other.get(i, d)));
		return result;
	}
	
	/** (Elementwise) Multiplication with a scalar */
	public @Dimensionless ComplexDoubleMatrix muli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v, @Dimensionless ComplexDoubleMatrix result) {
		ensureResultLength(null, result);
		
                @Dimensionless
                ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i, c).muli(v));
		return result;
	}

	public @Dimensionless ComplexDoubleMatrix muli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v, @Dimensionless ComplexDoubleMatrix result) {
		return muli(new @Dimensionless ComplexDouble(v), result);
	}

	/** Matrix-Matrix Multiplication */
	public @Dimensionless ComplexDoubleMatrix mmuli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
		if (other.isScalar())
			return muli(other.scalar(), result);

		/* check sizes and resize if necessary */
		assertMultipliesWith(other);
		if (result.rows != rows || result.columns != other.columns) {
			if (result != this && result != other)
				result.resize(rows, other.columns);
			else
				throw new @Dimensionless SizeException("Cannot resize result matrix because it is used in-place.");
		}
		
		if (result == this || result == other) {
			/* actually, blas cannot do multiplications in-place. Therefore, we will fake by
			 * allocating a temporary object on the side and copy the result later.
			 */
			@Dimensionless
			ComplexDoubleMatrix temp = new @Dimensionless ComplexDoubleMatrix(result.rows, result.columns);
			SimpleBlas.gemm(ComplexDouble.UNIT, this, other, ComplexDouble.ZERO, temp);
			SimpleBlas.copy(temp, result);
		}
		else {
			SimpleBlas.gemm(ComplexDouble.UNIT, this, other, ComplexDouble.ZERO, result);
		}		
		return result;
	}
	
	/** Matrix-Matrix Multiplication with a scalar (for symmetry, does the
	 * same as muli(scalar)
	 */
	public @Dimensionless ComplexDoubleMatrix mmuli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v, @Dimensionless ComplexDoubleMatrix result) {
		return muli(v, result);
	}

	public @Dimensionless ComplexDoubleMatrix mmuli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v, @Dimensionless ComplexDoubleMatrix result) {
		return muli(v, result);
	}
	
	/** (Elementwise) division */
	public @Dimensionless ComplexDoubleMatrix divi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
		if (other.isScalar())
			return divi(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);
		
                @Dimensionless
                ComplexDouble c1 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble c2 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i, c1).divi(other.get(i, c2)));
		return result;
	}
		
	/** (Elementwise) division with a scalar */
	public @Dimensionless ComplexDoubleMatrix divi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble a, @Dimensionless ComplexDoubleMatrix result) {
		ensureResultLength(null, result);
		
                @Dimensionless
                ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i, c).divi(a));
		return result;
	}	

	public @Dimensionless ComplexDoubleMatrix divi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double a, @Dimensionless ComplexDoubleMatrix result) {
		return divi(new @Dimensionless ComplexDouble(a), result);
	}

	/** 
	 * (Elementwise) division, with operands switched. Computes
	 * <em>result = other / this</em>. */
	public @Dimensionless ComplexDoubleMatrix rdivi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
		if (other.isScalar())
			return divi(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);

                @Dimensionless
                ComplexDouble c1 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble c2 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));

		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, other.get(i, c1).divi(get(i, c2)));
		return result;
	}
		
	/** (Elementwise) division with a scalar, with operands switched. Computes
	 * <em>result = a / this</em>.*/
	public @Dimensionless ComplexDoubleMatrix rdivi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble a, @Dimensionless ComplexDoubleMatrix result) {
		ensureResultLength(null, result);

                @Dimensionless
                ComplexDouble c1 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble c2 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));

		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                    c1.copy(a);
                    result.put(i, c1.divi(get(i, c2)));                    
                }
		return result;
	}

	public @Dimensionless ComplexDoubleMatrix rdivi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double a, @Dimensionless ComplexDoubleMatrix result) {
		return rdivi(new @Dimensionless ComplexDouble(a), result);
	}
	
	public @Dimensionless ComplexDoubleMatrix negi(@Dimensionless ComplexDoubleMatrix this) {
                @Dimensionless
                ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			put(i, get(i, c).negi());
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix neg(@Dimensionless ComplexDoubleMatrix this) {
		return dup().negi();
	}

	public @Dimensionless ComplexDoubleMatrix noti(@Dimensionless ComplexDoubleMatrix this) {
                @Dimensionless
                ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			put(i, get(i, c).isZero() ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix not(@Dimensionless ComplexDoubleMatrix this) {
		return dup().noti();
	}
	
	public @Dimensionless ComplexDoubleMatrix truthi(@Dimensionless ComplexDoubleMatrix this) {
                @Dimensionless
                ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			put(i, get(i, c).isZero() ? ((@Dimensionless double) (0.0)) : ((@Dimensionless double) (1.0)));
		return this;
	}
	
	public @Dimensionless ComplexDoubleMatrix truth(@Dimensionless ComplexDoubleMatrix this) {
		return dup().truthi();
	}

	/****************************************************************
	 * Rank one-updates
	 */
	
	/** Computes a rank-1-update A = A + alpha * x * y'. */ 
	public @Dimensionless ComplexDoubleMatrix rankOneUpdate(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble alpha, @Dimensionless ComplexDoubleMatrix x, @Dimensionless ComplexDoubleMatrix y) {
		if (rows != x.length)
			throw new @Dimensionless SizeException("Vector x has wrong length (" + x.length + " != " + rows + ").");
		if (columns != y.length)
			throw new @Dimensionless SizeException("Vector y has wrong length (" + x.length + " != " + columns + ").");			
		
		SimpleBlas.gerc(alpha, x, y, this);
		return this;
	}

	public @Dimensionless ComplexDoubleMatrix rankOneUpdate(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double alpha, @Dimensionless ComplexDoubleMatrix x, @Dimensionless ComplexDoubleMatrix y) {
		return rankOneUpdate(new @Dimensionless ComplexDouble(alpha), x, y);
	}

	/** Computes a rank-1-update A = A + alpha * x * x'. */ 
	public @Dimensionless ComplexDoubleMatrix rankOneUpdate(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double alpha, @Dimensionless ComplexDoubleMatrix x) {
		return rankOneUpdate(new @Dimensionless ComplexDouble(alpha), x, x);
	}

	/** Computes a rank-1-update A = A + alpha * x * x'. */ 
	public @Dimensionless ComplexDoubleMatrix rankOneUpdate(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble alpha, @Dimensionless ComplexDoubleMatrix x) {
		return rankOneUpdate(alpha, x, x);
	}

	/** Computes a rank-1-update A = A + x * x'. */ 
	public @Dimensionless ComplexDoubleMatrix rankOneUpdate(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix x) {
		return rankOneUpdate(((@Dimensionless double) (1.0)), x, x);
	}

	/** Computes a rank-1-update A = A + x * y'. */ 
	public @Dimensionless ComplexDoubleMatrix rankOneUpdate(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix x, @Dimensionless ComplexDoubleMatrix y) {
		return rankOneUpdate(((@Dimensionless double) (1.0)), x, y);
	}

	/****************************************************************
	 * Logical operations
	 */
	
	public @Dimensionless ComplexDouble sum(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		ComplexDouble s = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			s.addi(get(i, c));
		return s;
	}
	
	public @Dimensionless ComplexDouble mean(@Dimensionless ComplexDoubleMatrix this) {
		return sum().div((@Dimensionless double)length);
	}
	
	/** Computes this^T * other */
	public @Dimensionless ComplexDouble dotc(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
		return SimpleBlas.dotc(this, other);
	}
	
	/** Computes this^H * other */
	public @Dimensionless ComplexDouble dotu(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
		return SimpleBlas.dotu(this, other);
	}

	public @Dimensionless double norm2(@Dimensionless ComplexDoubleMatrix this) {
		return SimpleBlas.nrm2(this);
	}
	
	public @Dimensionless double normmax(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		int i = SimpleBlas.iamax(this);
		return get(i).abs();
	}

	public @Dimensionless double norm1(@Dimensionless ComplexDoubleMatrix this) {
		return SimpleBlas.asum(this);
	}
		
	/** Return a vector containing the sums of the columns (having number of columns many entries) */
	public @Dimensionless ComplexDoubleMatrix columnSums(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		ComplexDoubleMatrix v =
                        new @Dimensionless ComplexDoubleMatrix(((@Dimensionless int) (1)), columns);

		for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++)
			v.put(c, getColumn(c).sum());

		return v;
	}

	public @Dimensionless ComplexDoubleMatrix columnMeans(@Dimensionless ComplexDoubleMatrix this) {
		return columnSums().divi(rows);
	}
	
	public @Dimensionless ComplexDoubleMatrix rowSums(@Dimensionless ComplexDoubleMatrix this) {
		@Dimensionless
		ComplexDoubleMatrix v = new @Dimensionless ComplexDoubleMatrix(rows);

		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++)
			v.put(r, getRow(r).sum());

		return v;
	}

	public @Dimensionless ComplexDoubleMatrix rowMeans(@Dimensionless ComplexDoubleMatrix this) {
		return rowSums().divi(columns);
	}

	public @Dimensionless ComplexDoubleMatrix getColumn(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int c) {
		@Dimensionless
		ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(rows, ((@Dimensionless int) (1)));
		NativeBlas.zcopy(rows, data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)), result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		return result;
	}
	
	public void putColumn(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int c, @Dimensionless ComplexDoubleMatrix v) {
		NativeBlas.zcopy(rows, v.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)));
	}

	public @Dimensionless ComplexDoubleMatrix getRow(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r) {
		@Dimensionless
		ComplexDoubleMatrix result = new @Dimensionless ComplexDoubleMatrix(((@Dimensionless int) (1)), columns);
		NativeBlas.zcopy(columns, data, index(r, ((@Dimensionless int) (0))), rows, result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		return result;
	}
	
	public void putRow(@Dimensionless ComplexDoubleMatrix this, @Dimensionless int r, @Dimensionless ComplexDoubleMatrix v) {
		NativeBlas.zcopy(columns, v.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(r, ((@Dimensionless int) (0))), rows);
	}

	/**************************************************************************
	 * Elementwise Functions
	 */

	/** Add a row vector to all rows of the matrix */
	public void addRowVector(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix x) {
		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
			NativeBlas.zaxpy(columns, ComplexDouble.UNIT, x.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(r, ((@Dimensionless int) (0))), rows);
		}
	}

	/** Add a vector to all columns of the matrix */
	public void addColumnVector(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix x) {
		for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
			NativeBlas.zaxpy(rows, ComplexDouble.UNIT, x.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)));
		}
	}

       	/** Add a row vector to all rows of the matrix */
	public void subRowVector(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix x) {
		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
			NativeBlas.zaxpy(columns, ComplexDouble.NEG_UNIT, x.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(r, ((@Dimensionless int) (0))), rows);
		}
	}

	/** Add a vector to all columns of the matrix */
	public void subColumnVector(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix x) {
		for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
			NativeBlas.zaxpy(rows, ComplexDouble.NEG_UNIT, x.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)));
		}
	}

	/**
	 * Writes out this matrix to the given data stream.
	 * @param dos the data output stream to write to.
	 * @throws IOException 
	 */
	public void out(@Dimensionless ComplexDoubleMatrix this, @Dimensionless DataOutputStream dos) throws IOException {
		dos.writeUTF("double");
		dos.writeInt(columns);
		dos.writeInt(rows);
		
		dos.writeInt(data.length);
		for(@Dimensionless int i= ((@Dimensionless int) (0)); i < data.length;i++)
			dos.writeDouble(data[i]);
	}
	
	/**
	 * Reads in a matrix from the given data stream. Note
	 * that the old data of this matrix will be discarded.
	 * @param dis the data input stream to read from.
	 * @throws IOException 
	 */
	public void in(@Dimensionless DataInputStream dis) throws IOException {
		if(!dis.readUTF().equals("double")) 
			throw new @Dimensionless IllegalStateException("The matrix in the specified file is not of the correct type!");
		
		this.columns	= dis.readInt();
		this.rows		= dis.readInt();

		final @Dimensionless int MAX = dis.readInt();
		data = new @Dimensionless double @Dimensionless [MAX];
		for(@Dimensionless int i= ((@Dimensionless int) (0)); i < MAX;i++)
			data[i] = dis.readDouble();
	}	
	
	/**
	 * Saves this matrix to the specified file.
	 * @param filename the file to write the matrix in.
	 * @throws IOException thrown on errors while writing the matrix to the file
	 */
	public void save(@Dimensionless ComplexDoubleMatrix this, @Dimensionless String filename) throws IOException {
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
	public void load(@Dimensionless ComplexDoubleMatrix this, @Dimensionless String filename) throws IOException {
            @Dimensionless
            FileInputStream fis = new @Dimensionless FileInputStream(filename);
            @Dimensionless
            DataInputStream dis = new @Dimensionless DataInputStream(fis);
            try {
		this.in(dis);
            } finally {
                dis.close();
                fis.close();
            }
	}

	/****************************************************************
	 * Autogenerated code
	 */
	
	/***** Code for operators ***************************************/ 

	/* Overloads for the usual arithmetic operations */
	/*#
	 def gen_overloads(base, result_rows, result_cols); <<-EOS
	public ComplexDoubleMatrix #{base}i(ComplexDoubleMatrix other) {
		return #{base}i(other, this);
	}
	 	
	public ComplexDoubleMatrix #{base}(ComplexDoubleMatrix other) {
	  	return #{base}i(other, new ComplexDoubleMatrix(#{result_rows}, #{result_cols}));
	}

	public ComplexDoubleMatrix #{base}i(ComplexDouble v) {
		return #{base}i(v, this);
	}
	
	public ComplexDoubleMatrix #{base}i(double v) {
		return #{base}i(new ComplexDouble(v), this);
	}

	public ComplexDoubleMatrix #{base}(ComplexDouble v) {
		return #{base}i(v, new ComplexDoubleMatrix(rows, columns));
	} 	

	public ComplexDoubleMatrix #{base}(double v) {
		return #{base}i(new ComplexDouble(v), new ComplexDoubleMatrix(rows, columns));
	} 	
	 	
	 	EOS
	  end
	#*/

	/* Generating code for logical operators. This not only generates the stubs 
	 * but really all of the code.
	 */
	
	/*#
	 def gen_compare(name, op); <<-EOS
	 public ComplexDoubleMatrix #{name}i(ComplexDoubleMatrix other, ComplexDoubleMatrix result) {
	    if (other.isScalar())
	       return #{name}i(other.scalar(), result);
	       
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                ComplexDouble c1 = new ComplexDouble(0.0);
                ComplexDouble c2 = new ComplexDouble(0.0);
          
                for (int i = 0; i < length; i++)
                    result.put(i, get(i, c1).#{op}(other.get(i, c2)) ? 1.0 : 0.0);
	   return result;
	 }
	 
	 public ComplexDoubleMatrix #{name}i(ComplexDoubleMatrix other) {
	   return #{name}i(other, this);
	 }
	 
	 public ComplexDoubleMatrix #{name}(ComplexDoubleMatrix other) {
	   return #{name}i(other, new ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public ComplexDoubleMatrix #{name}i(ComplexDouble value, ComplexDoubleMatrix result) {
	   ensureResultLength(null, result);
           ComplexDouble c = new ComplexDouble(0.0);
	   for (int i = 0; i < length; i++)
	     result.put(i, get(i, c).#{op}(value) ? 1.0 : 0.0);
	   return result;
	 }

	 public ComplexDoubleMatrix #{name}i(double value, ComplexDoubleMatrix result) {
	   return #{name}i(new ComplexDouble(value), result);
	 }

	 public ComplexDoubleMatrix #{name}i(ComplexDouble value) {
	   return #{name}i(value, this);
	 }
	 
	 public ComplexDoubleMatrix #{name}i(double value) {
	   return #{name}i(new ComplexDouble(value));
	 }
	 
	 public ComplexDoubleMatrix #{name}(ComplexDouble value) {
	   return #{name}i(value, new ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public ComplexDoubleMatrix #{name}(double value) {
	   return #{name}i(new ComplexDouble(value));
	 }

	 EOS
	 end
	 #*/
	
	/*#
	 def gen_logical(name, op); <<-EOS
	 public ComplexDoubleMatrix #{name}i(ComplexDoubleMatrix other, ComplexDoubleMatrix result) {
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                ComplexDouble t1 = new ComplexDouble(0.0);
                ComplexDouble t2 = new ComplexDouble(0.0);
         
               for (int i = 0; i < length; i++)
                  result.put(i, (!get(i, t1).isZero()) #{op} (!other.get(i, t2).isZero()) ? 1.0 : 0.0);
	   return result;
	 }
	 
	 public ComplexDoubleMatrix #{name}i(ComplexDoubleMatrix other) {
	   return #{name}i(other, this);
	 }
	 
	 public ComplexDoubleMatrix #{name}(ComplexDoubleMatrix other) {
	   return #{name}i(other, new ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public ComplexDoubleMatrix #{name}i(ComplexDouble value, ComplexDoubleMatrix result) {
	 	ensureResultLength(null, result);
	 	boolean val = !value.isZero();
                ComplexDouble t = new ComplexDouble(0.0);
                for (int i = 0; i < length; i++)
                     result.put(i, !get(i, t).isZero() #{op} val ? 1.0 : 0.0);
	   return result;
	 }

 	 public ComplexDoubleMatrix #{name}i(double value, ComplexDoubleMatrix result) {
 	   return #{name}i(new ComplexDouble(value), result);
 	 }

	 public ComplexDoubleMatrix #{name}i(ComplexDouble value) {
	   return #{name}i(value, this);
	 }

 	 public ComplexDoubleMatrix #{name}i(double value) {
 	   return #{name}i(new ComplexDouble(value), this);
 	 }

	 public ComplexDoubleMatrix #{name}(ComplexDouble value) {
	   return #{name}i(value, new ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public ComplexDoubleMatrix #{name}(double value) {
	   return #{name}i(new ComplexDouble(value));
	 }
	 EOS
	 end
	 #*/

	/*# collect(gen_overloads('add', 'rows', 'columns'),
	  gen_overloads('sub', 'rows', 'columns'),
	  gen_overloads('rsub', 'rows', 'columns'),
	  gen_overloads('div', 'rows', 'columns'),
	  gen_overloads('rdiv', 'rows', 'columns'),
	  gen_overloads('mul', 'rows', 'columns'),
	  gen_overloads('mmul', 'rows', 'other.columns'),
	  gen_compare('eq', 'eq'),
	  gen_compare('ne', 'eq'),
	  gen_logical('and', '&'),
	  gen_logical('or', '|'),
	  gen_logical('xor', '^'))
	 #*/
//RJPP-BEGIN------------------------------------------------------------
	public @Dimensionless ComplexDoubleMatrix addi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
		return addi(other, this);
	}
	 	
	public @Dimensionless ComplexDoubleMatrix add(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	  	return addi(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	}

	public @Dimensionless ComplexDoubleMatrix addi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return addi(v, this);
	}
	
	public @Dimensionless ComplexDoubleMatrix addi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return addi(new @Dimensionless ComplexDouble(v), this);
	}

	public @Dimensionless ComplexDoubleMatrix add(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return addi(v, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexDoubleMatrix add(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return addi(new @Dimensionless ComplexDouble(v), new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexDoubleMatrix subi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
		return subi(other, this);
	}
	 	
	public @Dimensionless ComplexDoubleMatrix sub(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	  	return subi(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	}

	public @Dimensionless ComplexDoubleMatrix subi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return subi(v, this);
	}
	
	public @Dimensionless ComplexDoubleMatrix subi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return subi(new @Dimensionless ComplexDouble(v), this);
	}

	public @Dimensionless ComplexDoubleMatrix sub(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return subi(v, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexDoubleMatrix sub(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return subi(new @Dimensionless ComplexDouble(v), new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexDoubleMatrix rsubi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
		return rsubi(other, this);
	}
	 	
	public @Dimensionless ComplexDoubleMatrix rsub(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	  	return rsubi(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	}

	public @Dimensionless ComplexDoubleMatrix rsubi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return rsubi(v, this);
	}
	
	public @Dimensionless ComplexDoubleMatrix rsubi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return rsubi(new @Dimensionless ComplexDouble(v), this);
	}

	public @Dimensionless ComplexDoubleMatrix rsub(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return rsubi(v, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexDoubleMatrix rsub(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return rsubi(new @Dimensionless ComplexDouble(v), new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexDoubleMatrix divi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
		return divi(other, this);
	}
	 	
	public @Dimensionless ComplexDoubleMatrix div(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	  	return divi(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	}

	public @Dimensionless ComplexDoubleMatrix divi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return divi(v, this);
	}
	
	public @Dimensionless ComplexDoubleMatrix divi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return divi(new @Dimensionless ComplexDouble(v), this);
	}

	public @Dimensionless ComplexDoubleMatrix div(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return divi(v, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexDoubleMatrix div(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return divi(new @Dimensionless ComplexDouble(v), new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexDoubleMatrix rdivi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
		return rdivi(other, this);
	}
	 	
	public @Dimensionless ComplexDoubleMatrix rdiv(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	  	return rdivi(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	}

	public @Dimensionless ComplexDoubleMatrix rdivi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return rdivi(v, this);
	}
	
	public @Dimensionless ComplexDoubleMatrix rdivi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return rdivi(new @Dimensionless ComplexDouble(v), this);
	}

	public @Dimensionless ComplexDoubleMatrix rdiv(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return rdivi(v, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexDoubleMatrix rdiv(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return rdivi(new @Dimensionless ComplexDouble(v), new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexDoubleMatrix muli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
		return muli(other, this);
	}
	 	
	public @Dimensionless ComplexDoubleMatrix mul(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	  	return muli(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	}

	public @Dimensionless ComplexDoubleMatrix muli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return muli(v, this);
	}
	
	public @Dimensionless ComplexDoubleMatrix muli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return muli(new @Dimensionless ComplexDouble(v), this);
	}

	public @Dimensionless ComplexDoubleMatrix mul(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return muli(v, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexDoubleMatrix mul(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return muli(new @Dimensionless ComplexDouble(v), new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexDoubleMatrix mmuli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
		return mmuli(other, this);
	}
	 	
	public @Dimensionless ComplexDoubleMatrix mmul(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	  	return mmuli(other, new @Dimensionless ComplexDoubleMatrix(rows, other.columns));
	}

	public @Dimensionless ComplexDoubleMatrix mmuli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return mmuli(v, this);
	}
	
	public @Dimensionless ComplexDoubleMatrix mmuli(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return mmuli(new @Dimensionless ComplexDouble(v), this);
	}

	public @Dimensionless ComplexDoubleMatrix mmul(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble v) {
		return mmuli(v, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexDoubleMatrix mmul(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double v) {
		return mmuli(new @Dimensionless ComplexDouble(v), new @Dimensionless ComplexDoubleMatrix(rows, columns));
	} 	
	 	

	 public @Dimensionless ComplexDoubleMatrix eqi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
	    if (other.isScalar())
	       return eqi(other.scalar(), result);
	       
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexDouble c1 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble c2 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
          
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                    result.put(i, get(i, c1).eq(other.get(i, c2)) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix eqi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return eqi(other, this);
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix eq(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return eqi(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix eqi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value, @Dimensionless ComplexDoubleMatrix result) {
	   ensureResultLength(null, result);
           @Dimensionless
           ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
	     result.put(i, get(i, c).eq(value) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }

	 public @Dimensionless ComplexDoubleMatrix eqi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value, @Dimensionless ComplexDoubleMatrix result) {
	   return eqi(new @Dimensionless ComplexDouble(value), result);
	 }

	 public @Dimensionless ComplexDoubleMatrix eqi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return eqi(value, this);
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix eqi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
	   return eqi(new @Dimensionless ComplexDouble(value));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix eq(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return eqi(value, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix eq(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
	   return eqi(new @Dimensionless ComplexDouble(value));
	 }


	 public @Dimensionless ComplexDoubleMatrix nei(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
	    if (other.isScalar())
	       return nei(other.scalar(), result);
	       
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexDouble c1 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble c2 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
          
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                    result.put(i, get(i, c1).eq(other.get(i, c2)) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix nei(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return nei(other, this);
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix ne(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return nei(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix nei(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value, @Dimensionless ComplexDoubleMatrix result) {
	   ensureResultLength(null, result);
           @Dimensionless
           ComplexDouble c = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
	     result.put(i, get(i, c).eq(value) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }

	 public @Dimensionless ComplexDoubleMatrix nei(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value, @Dimensionless ComplexDoubleMatrix result) {
	   return nei(new @Dimensionless ComplexDouble(value), result);
	 }

	 public @Dimensionless ComplexDoubleMatrix nei(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return nei(value, this);
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix nei(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
	   return nei(new @Dimensionless ComplexDouble(value));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix ne(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return nei(value, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix ne(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
	   return nei(new @Dimensionless ComplexDouble(value));
	 }


	 public @Dimensionless ComplexDoubleMatrix andi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexDouble t1 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble t2 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
         
               for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                  result.put(i, (!get(i, t1).isZero()) & (!other.get(i, t2).isZero()) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix andi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return andi(other, this);
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix and(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return andi(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix andi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value, @Dimensionless ComplexDoubleMatrix result) {
	 	ensureResultLength(null, result);
	 	@Dimensionless
	 	boolean val = !value.isZero();
                @Dimensionless
                ComplexDouble t = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                     result.put(i, !get(i, t).isZero() & val ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }

 	 public @Dimensionless ComplexDoubleMatrix andi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value, @Dimensionless ComplexDoubleMatrix result) {
 	   return andi(new @Dimensionless ComplexDouble(value), result);
 	 }

	 public @Dimensionless ComplexDoubleMatrix andi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return andi(value, this);
	 }

 	 public @Dimensionless ComplexDoubleMatrix andi(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
 	   return andi(new @Dimensionless ComplexDouble(value), this);
 	 }

	 public @Dimensionless ComplexDoubleMatrix and(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return andi(value, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix and(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
	   return andi(new @Dimensionless ComplexDouble(value));
	 }

	 public @Dimensionless ComplexDoubleMatrix ori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexDouble t1 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble t2 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
         
               for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                  result.put(i, (!get(i, t1).isZero()) | (!other.get(i, t2).isZero()) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix ori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return ori(other, this);
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix or(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return ori(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix ori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value, @Dimensionless ComplexDoubleMatrix result) {
	 	ensureResultLength(null, result);
	 	@Dimensionless
	 	boolean val = !value.isZero();
                @Dimensionless
                ComplexDouble t = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                     result.put(i, !get(i, t).isZero() | val ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }

 	 public @Dimensionless ComplexDoubleMatrix ori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value, @Dimensionless ComplexDoubleMatrix result) {
 	   return ori(new @Dimensionless ComplexDouble(value), result);
 	 }

	 public @Dimensionless ComplexDoubleMatrix ori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return ori(value, this);
	 }

 	 public @Dimensionless ComplexDoubleMatrix ori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
 	   return ori(new @Dimensionless ComplexDouble(value), this);
 	 }

	 public @Dimensionless ComplexDoubleMatrix or(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return ori(value, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix or(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
	   return ori(new @Dimensionless ComplexDouble(value));
	 }

	 public @Dimensionless ComplexDoubleMatrix xori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other, @Dimensionless ComplexDoubleMatrix result) {
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexDouble t1 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                @Dimensionless
                ComplexDouble t2 = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
         
               for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                  result.put(i, (!get(i, t1).isZero()) ^ (!other.get(i, t2).isZero()) ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix xori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return xori(other, this);
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix xor(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDoubleMatrix other) {
	   return xori(other, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix xori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value, @Dimensionless ComplexDoubleMatrix result) {
	 	ensureResultLength(null, result);
	 	@Dimensionless
	 	boolean val = !value.isZero();
                @Dimensionless
                ComplexDouble t = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                     result.put(i, !get(i, t).isZero() ^ val ? ((@Dimensionless double) (1.0)) : ((@Dimensionless double) (0.0)));
	   return result;
	 }

 	 public @Dimensionless ComplexDoubleMatrix xori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value, @Dimensionless ComplexDoubleMatrix result) {
 	   return xori(new @Dimensionless ComplexDouble(value), result);
 	 }

	 public @Dimensionless ComplexDoubleMatrix xori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return xori(value, this);
	 }

 	 public @Dimensionless ComplexDoubleMatrix xori(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
 	   return xori(new @Dimensionless ComplexDouble(value), this);
 	 }

	 public @Dimensionless ComplexDoubleMatrix xor(@Dimensionless ComplexDoubleMatrix this, @Dimensionless ComplexDouble value) {
	   return xori(value, new @Dimensionless ComplexDoubleMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexDoubleMatrix xor(@Dimensionless ComplexDoubleMatrix this, @Dimensionless double value) {
	   return xori(new @Dimensionless ComplexDouble(value));
	 }
//RJPP-END--------------------------------------------------------------
}
