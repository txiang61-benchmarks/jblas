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

public class ComplexFloatMatrix {
	
	public @Dimensionless int rows;
	public @Dimensionless int columns;
	public @Dimensionless int length;
	public @Dimensionless float @Dimensionless [] data = null; // rows are contiguous

	/**************************************************************************
	 * 
	 * Constructors and factory functions
	 * 
	 **************************************************************************/

	/**
   * Create a new matrix with <i>newRows</i> rows, <i>newColumns</i> columns
	 * using <i>newData></i> as the data.
	 */
	public ComplexFloatMatrix(@Dimensionless int newRows, @Dimensionless int newColumns, @Dimensionless float @Dimensionless ... newData) {
		rows = newRows;
		columns = newColumns;
		length = rows * columns;

    if (newData.length != ((@Dimensionless int) (2)) * newRows * newColumns)
			throw new @Dimensionless IllegalArgumentException(
					"Passed data must match matrix dimensions.");
                data = newData;
	}
	
	/**
	 * Creates a new <i>n</i> times <i>m</i> <tt>ComplexFloatMatrix</tt>.
	 * @param newRows the number of rows (<i>n</i>) of the new matrix.
	 * @param newColumns the number of columns (<i>m</i>) of the new matrix.
	 */
	public ComplexFloatMatrix(@Dimensionless int newRows, @Dimensionless int newColumns) {
		this(newRows, newColumns, new @Dimensionless float @Dimensionless [((@Dimensionless int) (2)) * newRows * newColumns]);
	}
	
	/**
	 * Creates a new <tt>ComplexFloatMatrix</tt> of size 0 times 0.
	 */
	public ComplexFloatMatrix() {
		this(((@Dimensionless int) (0)), ((@Dimensionless int) (0)), null);
	}

	/**
	 * Create a Matrix of length <tt>len</tt>. By default, this creates a row vector.
	 * @param len
	 */
	public ComplexFloatMatrix(@Dimensionless int len) {
		this(len, ((@Dimensionless int) (1)), new @Dimensionless float @Dimensionless [((@Dimensionless int) (2)) * len]);
	}
	
	public ComplexFloatMatrix(@Dimensionless float @Dimensionless [] newData) {
		this(newData.length/ ((@Dimensionless int) (2)), ((@Dimensionless int) (1)), newData);
	}

	public ComplexFloatMatrix(@Dimensionless ComplexFloat @Dimensionless [] newData) {
		this(newData.length);
				
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < newData.length; i++)
			put(i, newData[i]);
	}
		
        
  /** Construct a complex matrix from a real matrix. */
  public ComplexFloatMatrix(FloatMatrix m) {
    this(m.rows, m.columns);

    NativeBlas.scopy(m.length, m.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, ((@Dimensionless int) (0)), ((@Dimensionless int) (2)));
  }

  /** Construct a complex matrix from separate real and imaginary parts. Either
   * part can be set to null in which case it will be ignored.
   */
  public ComplexFloatMatrix(@Dimensionless FloatMatrix real, @Dimensionless FloatMatrix imag) {
      this((real != null) ? real.rows : imag.rows, (real != null) ? real.columns : imag.columns);

      if (real != null && imag != null)
      real.assertSameSize(imag);

      if (real != null)
          NativeBlas.scopy(length, real.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, ((@Dimensionless int) (0)), ((@Dimensionless int) (2)));
      if (imag != null)
          NativeBlas.scopy(length, imag.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, ((@Dimensionless int) (1)), ((@Dimensionless int) (2)));
  }
        
        /**
	 * Creates a new matrix by reading it from a file.
	 * @param filename the path and name of the file to read the matrix from
	 * @throws IOException 
	 */
	public ComplexFloatMatrix(@Dimensionless String filename) throws IOException {
		load(filename);
	}
	
	/**
	 * Creates a new <i>n</i> times <i>m</i> <tt>ComplexFloatMatrix</tt> from
	 * the given <i>n</i> times <i>m</i> 2D data array. The first dimension of the array makes the
	 * rows (<i>n</i>) and the second dimension the columns (<i>m</i>). For example, the
	 * given code <br/><br/>
	 * <code>new ComplexFloatMatrix(new float[][]{{1d, 2d, 3d}, {4d, 5d, 6d}, {7d, 8d, 9d}}).print();</code><br/><br/>
	 * will constructs the following matrix:
	 * <pre>
	 * 1.0f	2.0f	3.0f
	 * 4.0f	5.0f	6.0f
	 * 7.0f	8.0f	9.0f
	 * </pre>.
	 * @param data <i>n</i> times <i>m</i> data array
	 */ 
	public ComplexFloatMatrix(@Dimensionless float @Dimensionless [] @Dimensionless [] data) {
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
	public static @Dimensionless ComplexFloatMatrix zeros(@Dimensionless int rows, @Dimensionless int columns) {
		return new @Dimensionless ComplexFloatMatrix(rows, columns);
	}
	
	public static @Dimensionless ComplexFloatMatrix zeros(@Dimensionless int length) {
		return zeros(length, ((@Dimensionless int) (1)));
	}

	/**
	 * Creates a new matrix in which all values are equal 1.
	 * @param rows number of rows
	 * @param columns number of columns
	 * @return new matrix
	 */
	public static @Dimensionless ComplexFloatMatrix ones(@Dimensionless int rows, @Dimensionless int columns) {
		@Dimensionless
		ComplexFloatMatrix m = new @Dimensionless ComplexFloatMatrix(rows, columns);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows * columns; i++)
			m.put(i, ((@Dimensionless float) (1.0f)));
		
		return m;
	}
	
	public static @Dimensionless ComplexFloatMatrix ones(@Dimensionless int length) {
		return ones(length, ((@Dimensionless int) (1)));
	}
	
	/**
	 * Creates a new matrix where the values of the given vector are the diagonal values of
	 * the matrix.
	 * @param x the diagonal values
	 * @return new matrix
	 */
	public static @Dimensionless ComplexFloatMatrix diag(@Dimensionless ComplexFloatMatrix x) {
		@Dimensionless
		ComplexFloatMatrix m = new @Dimensionless ComplexFloatMatrix(x.length, x.length);
		
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
  public static @Dimensionless ComplexFloatMatrix diag(@Dimensionless ComplexFloatMatrix x, @Dimensionless int rows, @Dimensionless int columns) {
    if (x.length > rows || x.length > columns) {
      throw new @Dimensionless SizeException("Length of diagonal matrix must be larger than both rows and columns.");
    }
    
    @Dimensionless
    ComplexFloatMatrix m = new @Dimensionless ComplexFloatMatrix(rows, columns);

    for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
      m.put(i, i, x.get(i));

    return m;
  }
	
	/**
	 * Create a 1 * 1 - matrix. For many operations, this matrix functions like a
	 * normal float
	 * @param s value of the matrix
	 * @return the constructed ComplexFloatMatrix 
	 */
	public static @Dimensionless ComplexFloatMatrix scalar(@Dimensionless float s) {
		@Dimensionless
		ComplexFloatMatrix m = new @Dimensionless ComplexFloatMatrix(((@Dimensionless int) (1)), ((@Dimensionless int) (1)));
		m.put(((@Dimensionless int) (0)), ((@Dimensionless int) (0)), s);
		return m;
	}
	
	/** Test whether a matrix is scalar */
	public @Dimensionless boolean isScalar(@Dimensionless ComplexFloatMatrix this) {
		return length == ((@Dimensionless int) (1));
	}
	
	/** Return the first element of the matrix */
	public @Dimensionless ComplexFloat scalar(@Dimensionless ComplexFloatMatrix this) {
		return get(((@Dimensionless int) (0)));
	}
	
	public static @Dimensionless ComplexFloatMatrix concatHorizontally(@Dimensionless ComplexFloatMatrix A, @Dimensionless ComplexFloatMatrix B) {
		if (A.rows != B.rows)
			throw new @Dimensionless SizeException("Matrices don't have same number of rows.");
		
		@Dimensionless
		ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(A.rows, A.columns + B.columns);
		SimpleBlas.copy(A, result);
		NativeBlas.ccopy(B.length, B.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), result.data, A.length, ((@Dimensionless int) (1)));
		return result;
	}

	public static @Dimensionless ComplexFloatMatrix concatVertically(@Dimensionless ComplexFloatMatrix A, @Dimensionless ComplexFloatMatrix B) {
		if (A.columns != B.columns)
			throw new @Dimensionless SizeException("Matrices don't have same number of columns.");
		
		@Dimensionless
		ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(A.rows + B.rows, A.columns);

		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < A.columns; i++) {
			NativeBlas.ccopy(A.rows, A.data, A.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)), result.data, result.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)));
			NativeBlas.ccopy(B.rows, B.data, B.index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)), result.data, result.index(A.rows, i), ((@Dimensionless int) (1)));
		}
		
		return result;
	}
	
	/**************************************************************************
	 * Working with slices (Man! 30+ methods just to make this a bit flexible...) 
	 */

	public @Dimensionless ComplexFloatMatrix get(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices) {
		@Dimensionless
		ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			result.put(i, get(indices[i]));
		
		return result;
	}
	
	public @Dimensionless ComplexFloatMatrix get(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices) {
		@Dimensionless
		ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(((@Dimensionless int) (1)), indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			result.put(i, get(r, indices[i]));
		
		return result;
	}
	
	public @Dimensionless ComplexFloatMatrix get(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c) {
		@Dimensionless
		ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(indices.length, ((@Dimensionless int) (1)));
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			result.put(i, get(indices[i], c));
		
		return result;
	}
	
	public @Dimensionless ComplexFloatMatrix get(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices) {
		@Dimensionless
		ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(rindices.length, cindices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				result.put(i, j, get(rindices[i], cindices[j]));
		
		return result;
	}
	
	public @Dimensionless ComplexFloatMatrix get(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices) {
		return get(indices.findIndices());
	}

	public @Dimensionless ComplexFloatMatrix get(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless ComplexFloatMatrix indices) {
		return get(r, indices.findIndices());
	}
	
	public @Dimensionless ComplexFloatMatrix get(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless int c) {
		return get(indices.findIndices(), c);
	}

	public @Dimensionless ComplexFloatMatrix get(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix rindices, @Dimensionless ComplexFloatMatrix cindices) {
		return get(rindices.findIndices(), cindices.findIndices());
	}
	
	private void checkLength(@Dimensionless ComplexFloatMatrix this, @Dimensionless int l) {
		if (length != l)
			throw new @Dimensionless SizeException("Matrix does not have the necessary length (" + length + " != " + l + ").");
	}

	private void checkRows(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r) {
		if (rows != r)
			throw new @Dimensionless SizeException("Matrix does not have the necessary length (" + length + " != " + r + ").");
	}
	
	private void checkColumns(@Dimensionless ComplexFloatMatrix this, @Dimensionless int c) {
		if (columns != c)
			throw new @Dimensionless SizeException("Matrix does not have the necessary length (" + length + " != " + c + ").");
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless ComplexFloatMatrix x) {
		if (x.isScalar())
			return put(indices, x.scalar());
		x.checkLength(indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], x.get(i));
		
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless ComplexFloatMatrix x) {
		if (x.isScalar())
			return put(r, indices, x.scalar());
		x.checkColumns(indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(r, indices[i], x.get(i));
		
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless ComplexFloatMatrix x) {
		if (x.isScalar())
			return put(indices, c, x.scalar());		
		x.checkRows(indices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], c, x.get(i));
		
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless ComplexFloatMatrix x) {
		if (x.isScalar())
			return put(rindices, cindices, x.scalar());		
		x.checkRows(rindices.length);
		x.checkColumns(cindices.length);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				put(rindices[i], cindices[j], x.get(i,j));
		
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless float v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless float v) {
		return put(indices, v);
	}

	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless float v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			putImag(indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless ComplexFloat v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless float v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(r, indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless float v) {
		return put(r, indices, v);
	}

	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless float v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			putImag(r, indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless int @Dimensionless [] indices, @Dimensionless ComplexFloat v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(r, indices[i], v);
		
		return this;
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless float v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], c, v);
		
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless float v) {
		return put(indices, c, v);
	}
	
	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless float v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			putImag(indices[i], c, v);
		
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] indices, @Dimensionless int c, @Dimensionless ComplexFloat v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < indices.length; i++)
			put(indices[i], c, v);
		
		return this;
 	}
	
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless float v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				put(rindices[i], cindices[j], v);
		
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless float v) {
		return put(rindices, cindices, v);
	}
	
	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless float v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				put(rindices[i], cindices[j], v);
		
		return this;
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int @Dimensionless [] rindices, @Dimensionless int @Dimensionless [] cindices, @Dimensionless ComplexFloat v) {
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rindices.length; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < cindices.length; j++)
				put(rindices[i], cindices[j], v);
		
		return this;
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless ComplexFloatMatrix v) {
		return put(indices.findIndices(), v);
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless ComplexFloatMatrix indices, @Dimensionless ComplexFloatMatrix v) {
		return put(r, indices.findIndices(), v);
	}
	
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless int c, @Dimensionless ComplexFloatMatrix v) {
		return put(indices.findIndices(), c, v);
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix rindices, @Dimensionless ComplexFloatMatrix cindices, @Dimensionless ComplexFloatMatrix v) {
		return put(rindices.findIndices(), cindices.findIndices(), v);
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless float v) {
		return put(indices.findIndices(), v);
	}

	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless float v) {
		return put(indices, v);
	}

	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless float v) {
		return putImag(indices.findIndices(), v);
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless ComplexFloat v) {
		return put(indices.findIndices(), v);
	}
	
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless ComplexFloatMatrix indices, @Dimensionless float v) {
		return put(r, indices.findIndices(), v);
	}
	
	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless ComplexFloatMatrix indices, @Dimensionless float v) {
		return put(r, indices, v);
	}

	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless ComplexFloatMatrix indices, @Dimensionless float v) {
		return putImag(r, indices.findIndices(), v);
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless ComplexFloatMatrix indices, @Dimensionless ComplexFloat v) {
		return put(r, indices.findIndices(), v);
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless int c, @Dimensionless float v) {
		return put(indices.findIndices(), c, v);
	}

	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless int c, @Dimensionless float v) {
		return put(indices, c, v);
	}

	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless int c, @Dimensionless float v) {
		return putImag(indices.findIndices(), c, v);
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix indices, @Dimensionless int c, @Dimensionless ComplexFloat v) {
		return put(indices.findIndices(), c, v);
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix rindices, @Dimensionless ComplexFloatMatrix cindices, @Dimensionless float v) {
		return put(rindices.findIndices(), cindices.findIndices(), v);
	}

	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix rindices, @Dimensionless ComplexFloatMatrix cindices, @Dimensionless float v) {
		return putReal(rindices.findIndices(), cindices.findIndices(), v);
	}

	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix rindices, @Dimensionless ComplexFloatMatrix cindices, @Dimensionless float v) {
		return putImag(rindices.findIndices(), cindices.findIndices(), v);
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix rindices, @Dimensionless ComplexFloatMatrix cindices, @Dimensionless ComplexFloat v) {
		return put(rindices.findIndices(), cindices.findIndices(), v);
	}

	
	public @Dimensionless int @Dimensionless [] findIndices(@Dimensionless ComplexFloatMatrix this) {
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
	public @Dimensionless ComplexFloatMatrix transpose(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(columns, rows);

                @Dimensionless
                ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless int) (0)));

		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows; i++)
			for (@Dimensionless int j = ((@Dimensionless int) (0)); j < columns; j++)
				result.put(j, i, get(i, j, c));
		
		return result;
	}

        public @Dimensionless ComplexFloatMatrix hermitian(@Dimensionless ComplexFloatMatrix this) {
            @Dimensionless
            ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(columns, rows);

            @Dimensionless
            ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless int) (0)));

            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < rows; i++)
                for (@Dimensionless int j = ((@Dimensionless int) (0)); j < columns; j++)
                    result.put(j, i, get(i, j, c).conji());
            return result;
        }

        /**
         * Compute complex conjugate (in-place).
         */
        public @Dimensionless ComplexFloatMatrix conji(@Dimensionless ComplexFloatMatrix this) {
            @Dimensionless
            ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                put(i, get(i, c).conji());
            return this;
        }

        /**
         * Compute complex conjugate.
         */
        public @Dimensionless ComplexFloatMatrix conj(@Dimensionless ComplexFloatMatrix this) {
            return dup().conji();
        }

		
	/** Compare two matrices.
	 * @param o Object to compare to
	 * @return true if and only if other is also a ComplexFloatMatrix which has the same size and the
	 * maximal absolute difference in matrix elements is smaller thatn 1e-6.  */
	public @Dimensionless boolean equals(@Dimensionless ComplexFloatMatrix this, @Dimensionless Object o) {
		if (!(o instanceof ComplexFloatMatrix))
			return false;

		@Dimensionless
		ComplexFloatMatrix other = (@Dimensionless ComplexFloatMatrix) o;

		if (!sameSize(other))
			return false;

    return Arrays.equals(data, other.data);
	}
  
  public @Dimensionless int hashCode(@Dimensionless ComplexFloatMatrix this) {
    return rows ^ columns ^ Arrays.hashCode(data);
  }

	
	/** Resize the matrix. All elements will be set to zero. */
	public void resize(@Dimensionless int newRows, @Dimensionless int newColumns) {
		rows = newRows;
		columns = newColumns;
		length = newRows * newColumns;
		data = new @Dimensionless float @Dimensionless [((@Dimensionless int) (2)) * rows * columns];
	}

	
	/** Reshape the matrix. Number of elements must not change. */
	public @Dimensionless ComplexFloatMatrix reshape(@Dimensionless ComplexFloatMatrix this, @Dimensionless int newRows, @Dimensionless int newColumns) {
		if (length != newRows * newColumns)
			throw new @Dimensionless IllegalArgumentException(
					"Number of elements must not change.");

		rows = newRows;
		columns = newColumns;
		
		return this;
	}

	/** Checks whether two matrices have the same size. */
	public @Dimensionless boolean sameSize(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix a) {
		return rows == a.rows && columns == a.columns;
	}

	/** 
	 * Assert that two matrices have the same size.
	 * 
	 * @param a the other matrix
	 * @throws SizeException if matrix sizes don't match. 
	 * */
	public void assertSameSize(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix a) {
		if (!sameSize(a))
			throw new @Dimensionless SizeException("Matrices must have the same size.");
	}
	
	/** 
	 * Check whether this can be multiplied with a. 
	 * 
	 * @param a right-hand-side of the multiplication.
	 * @return true iff <tt>this.columns == a.rows</tt>
	 */
	public @Dimensionless boolean multipliesWith(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix a) {
		return columns == a.rows;
	}
	
	public void assertMultipliesWith(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix a) {
		if (!multipliesWith(a))
			throw new @Dimensionless SizeException("Number of columns of left matrix must be equal to number of rows of right matrix.");
	}
	
	public @Dimensionless boolean sameLength(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix a) {
		return length == a.length;
	}
	
	public void assertSameLength(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix a) {
		if (!sameLength(a))
			throw new @Dimensionless SizeException("Matrices must have same length (is: " + length + " and " + a.length + ")");
	}
	
	/** Copy ComplexFloatMatrix a to this. this a is resized if necessary. */
	public @Dimensionless ComplexFloatMatrix copy(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix a) {
		if (!sameSize(a))
			resize(a.rows, a.columns);
		
		SimpleBlas.copy(a, this);
		return a;
	}
	
	/** Returns a duplicate of this matrix. Geometry is the same (including offsets, transpose, etc.),
	 * but the buffer is not shared.
	 */
	public @Dimensionless ComplexFloatMatrix dup(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		ComplexFloatMatrix out = new @Dimensionless ComplexFloatMatrix(rows, columns);

                JavaBlas.rcopy(((@Dimensionless int) (2))*length, data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), out.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		
		return out;
	}
	
	public @Dimensionless ComplexFloatMatrix swapColumns(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i, @Dimensionless int j) {
		NativeBlas.cswap(rows, data, index(((@Dimensionless int) (0)), i), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), j), ((@Dimensionless int) (1)));
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix swapRows(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i, @Dimensionless int j) {
		NativeBlas.cswap(columns, data, index(i, ((@Dimensionless int) (0))), rows, data, index(j, ((@Dimensionless int) (0))), rows);
		return this;
	}
		
	/** Set matrix element */
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless float value) {
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)] =  value;
		return this;
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless float realValue, @Dimensionless float complexValue) {
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)] =  realValue;
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)+ ((@Dimensionless int) (1))] =  complexValue;
		return this;
	}

        public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless ComplexFloat value) {
		@Dimensionless
		int i = ((@Dimensionless int) (2))*index(rowIndex, columnIndex);
		data[i] = value.real(); data[i+ ((@Dimensionless int) (1))] = value.imag();
		return this;
	}

	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless float value) {
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)] = value;
		return this;
	}

	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless float value) {
		data[((@Dimensionless int) (2))*index(rowIndex, columnIndex)+ ((@Dimensionless int) (1))] = value;
		return this;
	}
	
	/** Retrieve matrix element */
	public @Dimensionless ComplexFloat get(@Dimensionless ComplexFloatMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex) {
            @Dimensionless
            int i = ((@Dimensionless int) (2))*index(rowIndex, columnIndex);
            return new @Dimensionless ComplexFloat(data[i], data[i+ ((@Dimensionless int) (1))]);
	}

        /** Get matrix element, passing the variable to store the result. */
        public @Dimensionless ComplexFloat get(@Dimensionless ComplexFloatMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex, @Dimensionless ComplexFloat result) {
            return get(index(rowIndex, columnIndex), result);
        }
	
	public @Dimensionless FloatMatrix getReal(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		FloatMatrix result = new @Dimensionless FloatMatrix(rows, columns);
		
		NativeBlas.scopy(length, data, ((@Dimensionless int) (0)), ((@Dimensionless int) (2)), result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		
		return result;
	}

	/** Get index of an element */
	public @Dimensionless int index(@Dimensionless ComplexFloatMatrix this, @Dimensionless int rowIndex, @Dimensionless int columnIndex) {
		return rows * columnIndex + rowIndex;
	}

  /** Compute the row index of a linear index. */
  public @Dimensionless int indexRows(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i) {
    return i - indexColumns(i) * rows;
  }

  /** Compute the column index of a linear index. */
  public @Dimensionless int indexColumns(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i) {
    return i / rows;
  }


	public @Dimensionless ComplexFloat get(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i) {
		return new @Dimensionless ComplexFloat(data[i * ((@Dimensionless int) (2))], data[i * ((@Dimensionless int) (2)) + ((@Dimensionless int) (1))]);
	}
	
        public @Dimensionless ComplexFloat get(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i, @Dimensionless ComplexFloat result) {
            return result.set(data[i * ((@Dimensionless int) (2))], data[i* ((@Dimensionless int) (2))+ ((@Dimensionless int) (1))]);
        }
        
	public @Dimensionless float getReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i) {
		return data[((@Dimensionless int) (2))*i];
	}
	
	public @Dimensionless float getImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i) {
		return data[((@Dimensionless int) (2))*i + ((@Dimensionless int) (1))]; 
	}

	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i, @Dimensionless float v) {
		data[((@Dimensionless int) (2))*i] = v;
		return this;
	}

        public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i, @Dimensionless float r, @Dimensionless float c) {
            data[((@Dimensionless int) (2))*i] = r;
            data[((@Dimensionless int) (2))*i+ ((@Dimensionless int) (1))] = c;
            return this;
        }
	
	public @Dimensionless ComplexFloatMatrix put(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i, @Dimensionless ComplexFloat v) {
		data[((@Dimensionless int) (2))*i] = v.real();
		data[((@Dimensionless int) (2))*i+ ((@Dimensionless int) (1))] = v.imag();
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix putReal(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i, @Dimensionless float v) {
		return put(i, v);
	}
	
	public @Dimensionless ComplexFloatMatrix putImag(@Dimensionless ComplexFloatMatrix this, @Dimensionless int i, @Dimensionless float v) {
		data[((@Dimensionless int) (2))*i+ ((@Dimensionless int) (1))] = v;
		return this;
	}

	public @Dimensionless int getRows(@Dimensionless ComplexFloatMatrix this) {
		return rows;
	}
	
	public @Dimensionless int getColumns(@Dimensionless ComplexFloatMatrix this) {
		return columns;
	}
	
	public @Dimensionless int getLength(@Dimensionless ComplexFloatMatrix this) {
		return length;
	}
	
	/** Checks whether the matrix is empty. */
	public @Dimensionless boolean isEmpty(@Dimensionless ComplexFloatMatrix this) {
		return columns == ((@Dimensionless int) (0)) || rows == ((@Dimensionless int) (0));
	}
	
	/** Checks whether the matrix is square. */
	public @Dimensionless boolean isSquare(@Dimensionless ComplexFloatMatrix this) {
		return columns == rows;
	}
	
	public void assertSquare(@Dimensionless ComplexFloatMatrix this) {
		if (!isSquare())
			throw new @Dimensionless SizeException("Matrix must be square!");
	}
	
	/** Checks whether the matrix is a vector. */
	public @Dimensionless boolean isVector(@Dimensionless ComplexFloatMatrix this) {
		return columns == ((@Dimensionless int) (1)) || rows == ((@Dimensionless int) (1));
	}
	
	public @Dimensionless boolean isRowVector(@Dimensionless ComplexFloatMatrix this) {
		return columns == ((@Dimensionless int) (1));
	}
	
	public @Dimensionless boolean isColumnVector(@Dimensionless ComplexFloatMatrix this) {
		return rows == ((@Dimensionless int) (1));
	}
		
        /** Get diagonal of the matrix. */
	public @Dimensionless ComplexFloatMatrix diag(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		ComplexFloatMatrix d = new @Dimensionless ComplexFloatMatrix(rows);
		NativeBlas.ccopy(rows, data, ((@Dimensionless int) (0)), rows + ((@Dimensionless int) (1)), d.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		return d;
	}
        
        /** Get real part of the matrix. */
        public @Dimensionless FloatMatrix real(@Dimensionless ComplexFloatMatrix this) {
            @Dimensionless
            FloatMatrix result = new @Dimensionless FloatMatrix(rows, columns);
            NativeBlas.scopy(length, data, ((@Dimensionless int) (0)), ((@Dimensionless int) (2)), result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
            return result;
        }
        
        /** Get imaginary part of the matrix. */
        public @Dimensionless FloatMatrix imag(@Dimensionless ComplexFloatMatrix this) {
            @Dimensionless
            FloatMatrix result = new @Dimensionless FloatMatrix(rows, columns);
            NativeBlas.scopy(length, data, ((@Dimensionless int) (1)), ((@Dimensionless int) (2)), result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
            return result;            
        }

	
	/** 
	 * Pretty-print this matrix to <tt>System.out</tt>. 
	 * */
	public void print(@Dimensionless ComplexFloatMatrix this) {
		System.out.println(toString());
	}

	/** 
	 * Generate string representation of this matrix 
	 * (multi-line).
	 * */
	public @Dimensionless String toString(@Dimensionless ComplexFloatMatrix this) {
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

	public @Dimensionless float @Dimensionless [] toDoubleArray(@Dimensionless ComplexFloatMatrix this) {
		return data.clone();
	}
	
	public @Dimensionless ComplexFloat @Dimensionless [] toArray(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		ComplexFloat @Dimensionless [] array = new ComplexFloat @Dimensionless [length];
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			array[i] = get(i);
		
		return array;		
	}
	
	public @Dimensionless ComplexFloat @Dimensionless [] @Dimensionless [] toArray2(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		ComplexFloat @Dimensionless [] @Dimensionless [] array = new ComplexFloat @Dimensionless [rows][columns];
		
		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++)
			for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++)
				array[r][c] = get(r, c);
				
		return array;
	}
	
	public @Dimensionless boolean @Dimensionless [] toBooleanArray(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		boolean @Dimensionless [] array = new @Dimensionless boolean @Dimensionless [length];
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			array[i] = !get(i).isZero();
		
		return array;
	}
	
	public @Dimensionless boolean @Dimensionless [] @Dimensionless [] toBooleanArray2(@Dimensionless ComplexFloatMatrix this) {
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
	private void ensureResultLength(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
		if (!sameLength(result)) {
			if (result == this || result == other)
				throw new @Dimensionless SizeException("Cannot resize result matrix because it is used in-place.");
			result.resize(rows, columns);
		}
	}

	/** Add two matrices. */
	public @Dimensionless ComplexFloatMatrix addi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
		if (other.isScalar())
			return addi(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);
		
		if (result == this)
			SimpleBlas.axpy(ComplexFloat.UNIT, other, result);
		else if (result == other)
			SimpleBlas.axpy(ComplexFloat.UNIT, this, result);
		else {
			SimpleBlas.copy(this, result);
			SimpleBlas.axpy(ComplexFloat.UNIT, other, result);
		}

		return result;
	}
	
	/** Add a scalar to a matrix. */
	public @Dimensionless ComplexFloatMatrix addi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v, @Dimensionless ComplexFloatMatrix result) {
		ensureResultLength(null, result);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i).add(v));
		return result;
	}
	
	public @Dimensionless ComplexFloatMatrix addi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v, @Dimensionless ComplexFloatMatrix result) {
		return addi(new @Dimensionless ComplexFloat(v), result);
	}

	/** Subtract two matrices. */
	public @Dimensionless ComplexFloatMatrix subi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
		if (other.isScalar())
			return subi(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);
		
		if (result == this)
			SimpleBlas.axpy(ComplexFloat.NEG_UNIT, other, result);
		else if (result == other) {
			SimpleBlas.scal(ComplexFloat.NEG_UNIT, result);
			SimpleBlas.axpy(ComplexFloat.UNIT, this, result);
		}
		else {
			SimpleBlas.copy(this, result);
			SimpleBlas.axpy(ComplexFloat.NEG_UNIT, other, result);
		}
		return result;
	}
	
	/** Subtract a scalar from a matrix */
	public @Dimensionless ComplexFloatMatrix subi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v, @Dimensionless ComplexFloatMatrix result) {
		ensureResultLength(null, result);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i).sub(v));
		return result;
	}
	
	public @Dimensionless ComplexFloatMatrix subi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v, @Dimensionless ComplexFloatMatrix result) {
		return subi(new @Dimensionless ComplexFloat(v), result);
	}

	/** 
	 * Subtract two matrices, but subtract first from second matrix, that is, 
	 * compute <em>result = other - this</em>. 
	 * */
	public @Dimensionless ComplexFloatMatrix rsubi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
		return other.subi(this, result);
	}
	
	/** Subtract a matrix from a scalar */
	public @Dimensionless ComplexFloatMatrix rsubi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat a, @Dimensionless ComplexFloatMatrix result) {
		ensureResultLength(null, result);
		
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, a.sub(get(i)));
		return result;
	}

	public @Dimensionless ComplexFloatMatrix rsubi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float a, @Dimensionless ComplexFloatMatrix result) {
		return rsubi(new @Dimensionless ComplexFloat(a), result);
	}

	/** (Elementwise) Multiplication */ 
	public @Dimensionless ComplexFloatMatrix muli(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
		if (other.isScalar())
			return muli(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);
		
                @Dimensionless
                ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat d = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i, c).muli(other.get(i, d)));
		return result;
	}
	
	/** (Elementwise) Multiplication with a scalar */
	public @Dimensionless ComplexFloatMatrix muli(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v, @Dimensionless ComplexFloatMatrix result) {
		ensureResultLength(null, result);
		
                @Dimensionless
                ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i, c).muli(v));
		return result;
	}

	public @Dimensionless ComplexFloatMatrix muli(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v, @Dimensionless ComplexFloatMatrix result) {
		return muli(new @Dimensionless ComplexFloat(v), result);
	}

	/** Matrix-Matrix Multiplication */
	public @Dimensionless ComplexFloatMatrix mmuli(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
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
			ComplexFloatMatrix temp = new @Dimensionless ComplexFloatMatrix(result.rows, result.columns);
			SimpleBlas.gemm(ComplexFloat.UNIT, this, other, ComplexFloat.ZERO, temp);
			SimpleBlas.copy(temp, result);
		}
		else {
			SimpleBlas.gemm(ComplexFloat.UNIT, this, other, ComplexFloat.ZERO, result);
		}		
		return result;
	}
	
	/** Matrix-Matrix Multiplication with a scalar (for symmetry, does the
	 * same as muli(scalar)
	 */
	public @Dimensionless ComplexFloatMatrix mmuli(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v, @Dimensionless ComplexFloatMatrix result) {
		return muli(v, result);
	}

	public @Dimensionless ComplexFloatMatrix mmuli(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v, @Dimensionless ComplexFloatMatrix result) {
		return muli(v, result);
	}
	
	/** (Elementwise) division */
	public @Dimensionless ComplexFloatMatrix divi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
		if (other.isScalar())
			return divi(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);
		
                @Dimensionless
                ComplexFloat c1 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat c2 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i, c1).divi(other.get(i, c2)));
		return result;
	}
		
	/** (Elementwise) division with a scalar */
	public @Dimensionless ComplexFloatMatrix divi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat a, @Dimensionless ComplexFloatMatrix result) {
		ensureResultLength(null, result);
		
                @Dimensionless
                ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, get(i, c).divi(a));
		return result;
	}	

	public @Dimensionless ComplexFloatMatrix divi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float a, @Dimensionless ComplexFloatMatrix result) {
		return divi(new @Dimensionless ComplexFloat(a), result);
	}

	/** 
	 * (Elementwise) division, with operands switched. Computes
	 * <em>result = other / this</em>. */
	public @Dimensionless ComplexFloatMatrix rdivi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
		if (other.isScalar())
			return divi(other.scalar(), result);
		
		assertSameLength(other);
		ensureResultLength(other, result);

                @Dimensionless
                ComplexFloat c1 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat c2 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));

		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			result.put(i, other.get(i, c1).divi(get(i, c2)));
		return result;
	}
		
	/** (Elementwise) division with a scalar, with operands switched. Computes
	 * <em>result = a / this</em>.*/
	public @Dimensionless ComplexFloatMatrix rdivi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat a, @Dimensionless ComplexFloatMatrix result) {
		ensureResultLength(null, result);

                @Dimensionless
                ComplexFloat c1 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat c2 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));

		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++) {
                    c1.copy(a);
                    result.put(i, c1.divi(get(i, c2)));                    
                }
		return result;
	}

	public @Dimensionless ComplexFloatMatrix rdivi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float a, @Dimensionless ComplexFloatMatrix result) {
		return rdivi(new @Dimensionless ComplexFloat(a), result);
	}
	
	public @Dimensionless ComplexFloatMatrix negi(@Dimensionless ComplexFloatMatrix this) {
                @Dimensionless
                ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			put(i, get(i, c).negi());
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix neg(@Dimensionless ComplexFloatMatrix this) {
		return dup().negi();
	}

	public @Dimensionless ComplexFloatMatrix noti(@Dimensionless ComplexFloatMatrix this) {
                @Dimensionless
                ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			put(i, get(i, c).isZero() ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix not(@Dimensionless ComplexFloatMatrix this) {
		return dup().noti();
	}
	
	public @Dimensionless ComplexFloatMatrix truthi(@Dimensionless ComplexFloatMatrix this) {
                @Dimensionless
                ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			put(i, get(i, c).isZero() ? ((@Dimensionless float) (0.0f)) : ((@Dimensionless float) (1.0f)));
		return this;
	}
	
	public @Dimensionless ComplexFloatMatrix truth(@Dimensionless ComplexFloatMatrix this) {
		return dup().truthi();
	}

	/****************************************************************
	 * Rank one-updates
	 */
	
	/** Computes a rank-1-update A = A + alpha * x * y'. */ 
	public @Dimensionless ComplexFloatMatrix rankOneUpdate(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat alpha, @Dimensionless ComplexFloatMatrix x, @Dimensionless ComplexFloatMatrix y) {
		if (rows != x.length)
			throw new @Dimensionless SizeException("Vector x has wrong length (" + x.length + " != " + rows + ").");
		if (columns != y.length)
			throw new @Dimensionless SizeException("Vector y has wrong length (" + x.length + " != " + columns + ").");			
		
		SimpleBlas.gerc(alpha, x, y, this);
		return this;
	}

	public @Dimensionless ComplexFloatMatrix rankOneUpdate(@Dimensionless ComplexFloatMatrix this, @Dimensionless float alpha, @Dimensionless ComplexFloatMatrix x, @Dimensionless ComplexFloatMatrix y) {
		return rankOneUpdate(new @Dimensionless ComplexFloat(alpha), x, y);
	}

	/** Computes a rank-1-update A = A + alpha * x * x'. */ 
	public @Dimensionless ComplexFloatMatrix rankOneUpdate(@Dimensionless ComplexFloatMatrix this, @Dimensionless float alpha, @Dimensionless ComplexFloatMatrix x) {
		return rankOneUpdate(new @Dimensionless ComplexFloat(alpha), x, x);
	}

	/** Computes a rank-1-update A = A + alpha * x * x'. */ 
	public @Dimensionless ComplexFloatMatrix rankOneUpdate(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat alpha, @Dimensionless ComplexFloatMatrix x) {
		return rankOneUpdate(alpha, x, x);
	}

	/** Computes a rank-1-update A = A + x * x'. */ 
	public @Dimensionless ComplexFloatMatrix rankOneUpdate(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix x) {
		return rankOneUpdate(((@Dimensionless float) (1.0f)), x, x);
	}

	/** Computes a rank-1-update A = A + x * y'. */ 
	public @Dimensionless ComplexFloatMatrix rankOneUpdate(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix x, @Dimensionless ComplexFloatMatrix y) {
		return rankOneUpdate(((@Dimensionless float) (1.0f)), x, y);
	}

	/****************************************************************
	 * Logical operations
	 */
	
	public @Dimensionless ComplexFloat sum(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		ComplexFloat s = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
		for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
			s.addi(get(i, c));
		return s;
	}
	
	public @Dimensionless ComplexFloat mean(@Dimensionless ComplexFloatMatrix this) {
		return sum().div((@Dimensionless float)length);
	}
	
	/** Computes this^T * other */
	public @Dimensionless ComplexFloat dotc(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
		return SimpleBlas.dotc(this, other);
	}
	
	/** Computes this^H * other */
	public @Dimensionless ComplexFloat dotu(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
		return SimpleBlas.dotu(this, other);
	}

	public @Dimensionless float norm2(@Dimensionless ComplexFloatMatrix this) {
		return SimpleBlas.nrm2(this);
	}
	
	public @Dimensionless float normmax(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		int i = SimpleBlas.iamax(this);
		return get(i).abs();
	}

	public @Dimensionless float norm1(@Dimensionless ComplexFloatMatrix this) {
		return SimpleBlas.asum(this);
	}
		
	/** Return a vector containing the sums of the columns (having number of columns many entries) */
	public @Dimensionless ComplexFloatMatrix columnSums(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		ComplexFloatMatrix v =
                        new @Dimensionless ComplexFloatMatrix(((@Dimensionless int) (1)), columns);

		for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++)
			v.put(c, getColumn(c).sum());

		return v;
	}

	public @Dimensionless ComplexFloatMatrix columnMeans(@Dimensionless ComplexFloatMatrix this) {
		return columnSums().divi(rows);
	}
	
	public @Dimensionless ComplexFloatMatrix rowSums(@Dimensionless ComplexFloatMatrix this) {
		@Dimensionless
		ComplexFloatMatrix v = new @Dimensionless ComplexFloatMatrix(rows);

		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++)
			v.put(r, getRow(r).sum());

		return v;
	}

	public @Dimensionless ComplexFloatMatrix rowMeans(@Dimensionless ComplexFloatMatrix this) {
		return rowSums().divi(columns);
	}

	public @Dimensionless ComplexFloatMatrix getColumn(@Dimensionless ComplexFloatMatrix this, @Dimensionless int c) {
		@Dimensionless
		ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(rows, ((@Dimensionless int) (1)));
		NativeBlas.ccopy(rows, data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)), result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		return result;
	}
	
	public void putColumn(@Dimensionless ComplexFloatMatrix this, @Dimensionless int c, @Dimensionless ComplexFloatMatrix v) {
		NativeBlas.ccopy(rows, v.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)));
	}

	public @Dimensionless ComplexFloatMatrix getRow(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r) {
		@Dimensionless
		ComplexFloatMatrix result = new @Dimensionless ComplexFloatMatrix(((@Dimensionless int) (1)), columns);
		NativeBlas.ccopy(columns, data, index(r, ((@Dimensionless int) (0))), rows, result.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
		return result;
	}
	
	public void putRow(@Dimensionless ComplexFloatMatrix this, @Dimensionless int r, @Dimensionless ComplexFloatMatrix v) {
		NativeBlas.ccopy(columns, v.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(r, ((@Dimensionless int) (0))), rows);
	}

	/**************************************************************************
	 * Elementwise Functions
	 */

	/** Add a row vector to all rows of the matrix */
	public void addRowVector(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix x) {
		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
			NativeBlas.caxpy(columns, ComplexFloat.UNIT, x.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(r, ((@Dimensionless int) (0))), rows);
		}
	}

	/** Add a vector to all columns of the matrix */
	public void addColumnVector(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix x) {
		for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
			NativeBlas.caxpy(rows, ComplexFloat.UNIT, x.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)));
		}
	}

       	/** Add a row vector to all rows of the matrix */
	public void subRowVector(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix x) {
		for (@Dimensionless int r = ((@Dimensionless int) (0)); r < rows; r++) {
			NativeBlas.caxpy(columns, ComplexFloat.NEG_UNIT, x.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(r, ((@Dimensionless int) (0))), rows);
		}
	}

	/** Add a vector to all columns of the matrix */
	public void subColumnVector(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix x) {
		for (@Dimensionless int c = ((@Dimensionless int) (0)); c < columns; c++) {
			NativeBlas.caxpy(rows, ComplexFloat.NEG_UNIT, x.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), data, index(((@Dimensionless int) (0)), c), ((@Dimensionless int) (1)));
		}
	}

	/**
	 * Writes out this matrix to the given data stream.
	 * @param dos the data output stream to write to.
	 * @throws IOException 
	 */
	public void out(@Dimensionless ComplexFloatMatrix this, @Dimensionless DataOutputStream dos) throws IOException {
		dos.writeUTF("float");
		dos.writeInt(columns);
		dos.writeInt(rows);
		
		dos.writeInt(data.length);
		for(@Dimensionless int i= ((@Dimensionless int) (0)); i < data.length;i++)
			dos.writeFloat(data[i]);
	}
	
	/**
	 * Reads in a matrix from the given data stream. Note
	 * that the old data of this matrix will be discarded.
	 * @param dis the data input stream to read from.
	 * @throws IOException 
	 */
	public void in(@Dimensionless DataInputStream dis) throws IOException {
		if(!dis.readUTF().equals("float")) 
			throw new @Dimensionless IllegalStateException("The matrix in the specified file is not of the correct type!");
		
		this.columns	= dis.readInt();
		this.rows		= dis.readInt();

		final @Dimensionless int MAX = dis.readInt();
		data = new @Dimensionless float @Dimensionless [MAX];
		for(@Dimensionless int i= ((@Dimensionless int) (0)); i < MAX;i++)
			data[i] = dis.readFloat();
	}	
	
	/**
	 * Saves this matrix to the specified file.
	 * @param filename the file to write the matrix in.
	 * @throws IOException thrown on errors while writing the matrix to the file
	 */
	public void save(@Dimensionless ComplexFloatMatrix this, @Dimensionless String filename) throws IOException {
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
	public void load(@Dimensionless ComplexFloatMatrix this, @Dimensionless String filename) throws IOException {
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
	public ComplexFloatMatrix #{base}i(ComplexFloatMatrix other) {
		return #{base}i(other, this);
	}
	 	
	public ComplexFloatMatrix #{base}(ComplexFloatMatrix other) {
	  	return #{base}i(other, new ComplexFloatMatrix(#{result_rows}, #{result_cols}));
	}

	public ComplexFloatMatrix #{base}i(ComplexFloat v) {
		return #{base}i(v, this);
	}
	
	public ComplexFloatMatrix #{base}i(float v) {
		return #{base}i(new ComplexFloat(v), this);
	}

	public ComplexFloatMatrix #{base}(ComplexFloat v) {
		return #{base}i(v, new ComplexFloatMatrix(rows, columns));
	} 	

	public ComplexFloatMatrix #{base}(float v) {
		return #{base}i(new ComplexFloat(v), new ComplexFloatMatrix(rows, columns));
	} 	
	 	
	 	EOS
	  end
	#*/

	/* Generating code for logical operators. This not only generates the stubs 
	 * but really all of the code.
	 */
	
	/*#
	 def gen_compare(name, op); <<-EOS
	 public ComplexFloatMatrix #{name}i(ComplexFloatMatrix other, ComplexFloatMatrix result) {
	    if (other.isScalar())
	       return #{name}i(other.scalar(), result);
	       
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                ComplexFloat c1 = new ComplexFloat(0.0f);
                ComplexFloat c2 = new ComplexFloat(0.0f);
          
                for (int i = 0; i < length; i++)
                    result.put(i, get(i, c1).#{op}(other.get(i, c2)) ? 1.0f : 0.0f);
	   return result;
	 }
	 
	 public ComplexFloatMatrix #{name}i(ComplexFloatMatrix other) {
	   return #{name}i(other, this);
	 }
	 
	 public ComplexFloatMatrix #{name}(ComplexFloatMatrix other) {
	   return #{name}i(other, new ComplexFloatMatrix(rows, columns));
	 }
	 
	 public ComplexFloatMatrix #{name}i(ComplexFloat value, ComplexFloatMatrix result) {
	   ensureResultLength(null, result);
           ComplexFloat c = new ComplexFloat(0.0f);
	   for (int i = 0; i < length; i++)
	     result.put(i, get(i, c).#{op}(value) ? 1.0f : 0.0f);
	   return result;
	 }

	 public ComplexFloatMatrix #{name}i(float value, ComplexFloatMatrix result) {
	   return #{name}i(new ComplexFloat(value), result);
	 }

	 public ComplexFloatMatrix #{name}i(ComplexFloat value) {
	   return #{name}i(value, this);
	 }
	 
	 public ComplexFloatMatrix #{name}i(float value) {
	   return #{name}i(new ComplexFloat(value));
	 }
	 
	 public ComplexFloatMatrix #{name}(ComplexFloat value) {
	   return #{name}i(value, new ComplexFloatMatrix(rows, columns));
	 }
	 
	 public ComplexFloatMatrix #{name}(float value) {
	   return #{name}i(new ComplexFloat(value));
	 }

	 EOS
	 end
	 #*/
	
	/*#
	 def gen_logical(name, op); <<-EOS
	 public ComplexFloatMatrix #{name}i(ComplexFloatMatrix other, ComplexFloatMatrix result) {
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                ComplexFloat t1 = new ComplexFloat(0.0f);
                ComplexFloat t2 = new ComplexFloat(0.0f);
         
               for (int i = 0; i < length; i++)
                  result.put(i, (!get(i, t1).isZero()) #{op} (!other.get(i, t2).isZero()) ? 1.0f : 0.0f);
	   return result;
	 }
	 
	 public ComplexFloatMatrix #{name}i(ComplexFloatMatrix other) {
	   return #{name}i(other, this);
	 }
	 
	 public ComplexFloatMatrix #{name}(ComplexFloatMatrix other) {
	   return #{name}i(other, new ComplexFloatMatrix(rows, columns));
	 }
	 
	 public ComplexFloatMatrix #{name}i(ComplexFloat value, ComplexFloatMatrix result) {
	 	ensureResultLength(null, result);
	 	boolean val = !value.isZero();
                ComplexFloat t = new ComplexFloat(0.0f);
                for (int i = 0; i < length; i++)
                     result.put(i, !get(i, t).isZero() #{op} val ? 1.0f : 0.0f);
	   return result;
	 }

 	 public ComplexFloatMatrix #{name}i(float value, ComplexFloatMatrix result) {
 	   return #{name}i(new ComplexFloat(value), result);
 	 }

	 public ComplexFloatMatrix #{name}i(ComplexFloat value) {
	   return #{name}i(value, this);
	 }

 	 public ComplexFloatMatrix #{name}i(float value) {
 	   return #{name}i(new ComplexFloat(value), this);
 	 }

	 public ComplexFloatMatrix #{name}(ComplexFloat value) {
	   return #{name}i(value, new ComplexFloatMatrix(rows, columns));
	 }
	 
	 public ComplexFloatMatrix #{name}(float value) {
	   return #{name}i(new ComplexFloat(value));
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
	public @Dimensionless ComplexFloatMatrix addi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
		return addi(other, this);
	}
	 	
	public @Dimensionless ComplexFloatMatrix add(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	  	return addi(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	}

	public @Dimensionless ComplexFloatMatrix addi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return addi(v, this);
	}
	
	public @Dimensionless ComplexFloatMatrix addi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return addi(new @Dimensionless ComplexFloat(v), this);
	}

	public @Dimensionless ComplexFloatMatrix add(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return addi(v, new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexFloatMatrix add(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return addi(new @Dimensionless ComplexFloat(v), new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexFloatMatrix subi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
		return subi(other, this);
	}
	 	
	public @Dimensionless ComplexFloatMatrix sub(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	  	return subi(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	}

	public @Dimensionless ComplexFloatMatrix subi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return subi(v, this);
	}
	
	public @Dimensionless ComplexFloatMatrix subi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return subi(new @Dimensionless ComplexFloat(v), this);
	}

	public @Dimensionless ComplexFloatMatrix sub(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return subi(v, new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexFloatMatrix sub(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return subi(new @Dimensionless ComplexFloat(v), new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexFloatMatrix rsubi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
		return rsubi(other, this);
	}
	 	
	public @Dimensionless ComplexFloatMatrix rsub(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	  	return rsubi(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	}

	public @Dimensionless ComplexFloatMatrix rsubi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return rsubi(v, this);
	}
	
	public @Dimensionless ComplexFloatMatrix rsubi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return rsubi(new @Dimensionless ComplexFloat(v), this);
	}

	public @Dimensionless ComplexFloatMatrix rsub(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return rsubi(v, new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexFloatMatrix rsub(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return rsubi(new @Dimensionless ComplexFloat(v), new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexFloatMatrix divi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
		return divi(other, this);
	}
	 	
	public @Dimensionless ComplexFloatMatrix div(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	  	return divi(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	}

	public @Dimensionless ComplexFloatMatrix divi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return divi(v, this);
	}
	
	public @Dimensionless ComplexFloatMatrix divi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return divi(new @Dimensionless ComplexFloat(v), this);
	}

	public @Dimensionless ComplexFloatMatrix div(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return divi(v, new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexFloatMatrix div(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return divi(new @Dimensionless ComplexFloat(v), new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexFloatMatrix rdivi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
		return rdivi(other, this);
	}
	 	
	public @Dimensionless ComplexFloatMatrix rdiv(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	  	return rdivi(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	}

	public @Dimensionless ComplexFloatMatrix rdivi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return rdivi(v, this);
	}
	
	public @Dimensionless ComplexFloatMatrix rdivi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return rdivi(new @Dimensionless ComplexFloat(v), this);
	}

	public @Dimensionless ComplexFloatMatrix rdiv(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return rdivi(v, new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexFloatMatrix rdiv(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return rdivi(new @Dimensionless ComplexFloat(v), new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexFloatMatrix muli(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
		return muli(other, this);
	}
	 	
	public @Dimensionless ComplexFloatMatrix mul(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	  	return muli(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	}

	public @Dimensionless ComplexFloatMatrix muli(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return muli(v, this);
	}
	
	public @Dimensionless ComplexFloatMatrix muli(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return muli(new @Dimensionless ComplexFloat(v), this);
	}

	public @Dimensionless ComplexFloatMatrix mul(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return muli(v, new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexFloatMatrix mul(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return muli(new @Dimensionless ComplexFloat(v), new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	
	 	

	public @Dimensionless ComplexFloatMatrix mmuli(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
		return mmuli(other, this);
	}
	 	
	public @Dimensionless ComplexFloatMatrix mmul(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	  	return mmuli(other, new @Dimensionless ComplexFloatMatrix(rows, other.columns));
	}

	public @Dimensionless ComplexFloatMatrix mmuli(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return mmuli(v, this);
	}
	
	public @Dimensionless ComplexFloatMatrix mmuli(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return mmuli(new @Dimensionless ComplexFloat(v), this);
	}

	public @Dimensionless ComplexFloatMatrix mmul(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat v) {
		return mmuli(v, new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	

	public @Dimensionless ComplexFloatMatrix mmul(@Dimensionless ComplexFloatMatrix this, @Dimensionless float v) {
		return mmuli(new @Dimensionless ComplexFloat(v), new @Dimensionless ComplexFloatMatrix(rows, columns));
	} 	
	 	

	 public @Dimensionless ComplexFloatMatrix eqi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
	    if (other.isScalar())
	       return eqi(other.scalar(), result);
	       
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexFloat c1 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat c2 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
          
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                    result.put(i, get(i, c1).eq(other.get(i, c2)) ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexFloatMatrix eqi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return eqi(other, this);
	 }
	 
	 public @Dimensionless ComplexFloatMatrix eq(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return eqi(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix eqi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value, @Dimensionless ComplexFloatMatrix result) {
	   ensureResultLength(null, result);
           @Dimensionless
           ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
	     result.put(i, get(i, c).eq(value) ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }

	 public @Dimensionless ComplexFloatMatrix eqi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value, @Dimensionless ComplexFloatMatrix result) {
	   return eqi(new @Dimensionless ComplexFloat(value), result);
	 }

	 public @Dimensionless ComplexFloatMatrix eqi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return eqi(value, this);
	 }
	 
	 public @Dimensionless ComplexFloatMatrix eqi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
	   return eqi(new @Dimensionless ComplexFloat(value));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix eq(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return eqi(value, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix eq(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
	   return eqi(new @Dimensionless ComplexFloat(value));
	 }


	 public @Dimensionless ComplexFloatMatrix nei(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
	    if (other.isScalar())
	       return nei(other.scalar(), result);
	       
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexFloat c1 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat c2 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
          
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                    result.put(i, get(i, c1).eq(other.get(i, c2)) ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexFloatMatrix nei(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return nei(other, this);
	 }
	 
	 public @Dimensionless ComplexFloatMatrix ne(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return nei(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix nei(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value, @Dimensionless ComplexFloatMatrix result) {
	   ensureResultLength(null, result);
           @Dimensionless
           ComplexFloat c = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
	     result.put(i, get(i, c).eq(value) ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }

	 public @Dimensionless ComplexFloatMatrix nei(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value, @Dimensionless ComplexFloatMatrix result) {
	   return nei(new @Dimensionless ComplexFloat(value), result);
	 }

	 public @Dimensionless ComplexFloatMatrix nei(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return nei(value, this);
	 }
	 
	 public @Dimensionless ComplexFloatMatrix nei(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
	   return nei(new @Dimensionless ComplexFloat(value));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix ne(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return nei(value, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix ne(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
	   return nei(new @Dimensionless ComplexFloat(value));
	 }


	 public @Dimensionless ComplexFloatMatrix andi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexFloat t1 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat t2 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
         
               for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                  result.put(i, (!get(i, t1).isZero()) & (!other.get(i, t2).isZero()) ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexFloatMatrix andi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return andi(other, this);
	 }
	 
	 public @Dimensionless ComplexFloatMatrix and(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return andi(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix andi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value, @Dimensionless ComplexFloatMatrix result) {
	 	ensureResultLength(null, result);
	 	@Dimensionless
	 	boolean val = !value.isZero();
                @Dimensionless
                ComplexFloat t = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                     result.put(i, !get(i, t).isZero() & val ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }

 	 public @Dimensionless ComplexFloatMatrix andi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value, @Dimensionless ComplexFloatMatrix result) {
 	   return andi(new @Dimensionless ComplexFloat(value), result);
 	 }

	 public @Dimensionless ComplexFloatMatrix andi(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return andi(value, this);
	 }

 	 public @Dimensionless ComplexFloatMatrix andi(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
 	   return andi(new @Dimensionless ComplexFloat(value), this);
 	 }

	 public @Dimensionless ComplexFloatMatrix and(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return andi(value, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix and(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
	   return andi(new @Dimensionless ComplexFloat(value));
	 }

	 public @Dimensionless ComplexFloatMatrix ori(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexFloat t1 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat t2 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
         
               for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                  result.put(i, (!get(i, t1).isZero()) | (!other.get(i, t2).isZero()) ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexFloatMatrix ori(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return ori(other, this);
	 }
	 
	 public @Dimensionless ComplexFloatMatrix or(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return ori(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix ori(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value, @Dimensionless ComplexFloatMatrix result) {
	 	ensureResultLength(null, result);
	 	@Dimensionless
	 	boolean val = !value.isZero();
                @Dimensionless
                ComplexFloat t = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                     result.put(i, !get(i, t).isZero() | val ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }

 	 public @Dimensionless ComplexFloatMatrix ori(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value, @Dimensionless ComplexFloatMatrix result) {
 	   return ori(new @Dimensionless ComplexFloat(value), result);
 	 }

	 public @Dimensionless ComplexFloatMatrix ori(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return ori(value, this);
	 }

 	 public @Dimensionless ComplexFloatMatrix ori(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
 	   return ori(new @Dimensionless ComplexFloat(value), this);
 	 }

	 public @Dimensionless ComplexFloatMatrix or(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return ori(value, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix or(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
	   return ori(new @Dimensionless ComplexFloat(value));
	 }

	 public @Dimensionless ComplexFloatMatrix xori(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other, @Dimensionless ComplexFloatMatrix result) {
	 	assertSameLength(other);
	 	ensureResultLength(other, result);
	 	
                @Dimensionless
                ComplexFloat t1 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                @Dimensionless
                ComplexFloat t2 = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
         
               for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                  result.put(i, (!get(i, t1).isZero()) ^ (!other.get(i, t2).isZero()) ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }
	 
	 public @Dimensionless ComplexFloatMatrix xori(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return xori(other, this);
	 }
	 
	 public @Dimensionless ComplexFloatMatrix xor(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloatMatrix other) {
	   return xori(other, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix xori(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value, @Dimensionless ComplexFloatMatrix result) {
	 	ensureResultLength(null, result);
	 	@Dimensionless
	 	boolean val = !value.isZero();
                @Dimensionless
                ComplexFloat t = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < length; i++)
                     result.put(i, !get(i, t).isZero() ^ val ? ((@Dimensionless float) (1.0f)) : ((@Dimensionless float) (0.0f)));
	   return result;
	 }

 	 public @Dimensionless ComplexFloatMatrix xori(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value, @Dimensionless ComplexFloatMatrix result) {
 	   return xori(new @Dimensionless ComplexFloat(value), result);
 	 }

	 public @Dimensionless ComplexFloatMatrix xori(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return xori(value, this);
	 }

 	 public @Dimensionless ComplexFloatMatrix xori(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
 	   return xori(new @Dimensionless ComplexFloat(value), this);
 	 }

	 public @Dimensionless ComplexFloatMatrix xor(@Dimensionless ComplexFloatMatrix this, @Dimensionless ComplexFloat value) {
	   return xori(value, new @Dimensionless ComplexFloatMatrix(rows, columns));
	 }
	 
	 public @Dimensionless ComplexFloatMatrix xor(@Dimensionless ComplexFloatMatrix this, @Dimensionless float value) {
	   return xori(new @Dimensionless ComplexFloat(value));
	 }
//RJPP-END--------------------------------------------------------------
}