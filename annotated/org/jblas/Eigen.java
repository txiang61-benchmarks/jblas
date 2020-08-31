// --- BEGIN LICENSE BLOCK ---
/* 
 * Copyright (c) 2009-13, Mikio L. Braun
 *               2011, Nicola Oury
 *               2013, Alexander Sehlström
 *
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
import org.jblas.exceptions.NoEigenResultException;
import org.jblas.ranges.IntervalRange;
import org.jblas.ranges.Range;

/**
 * <p>Eigenvalue and Eigenvector related functions.</p>
 * <p/>
 * <p>Methods exist for working with symmetric matrices or general eigenvalues.
 * The symmetric versions are usually much faster on symmetric matrices.</p>
 */
@Dimensionless
public class Eigen {
    private static final @Dimensionless DoubleMatrix dummyDouble = new @Dimensionless DoubleMatrix(((@Dimensionless int) (1)));

    /**
     * Compute the eigenvalues for a symmetric matrix.
     */
    public static @Dimensionless DoubleMatrix symmetricEigenvalues(@Dimensionless DoubleMatrix A) {
        A.assertSquare();
        @Dimensionless
        DoubleMatrix eigenvalues = new @Dimensionless DoubleMatrix(A.rows);
        @Dimensionless
        int isuppz @Dimensionless [] = new int @Dimensionless [((@Dimensionless int) (2)) * A.rows];
        SimpleBlas.syevr('N', 'A', 'U', A.dup(), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), eigenvalues, dummyDouble, isuppz);
        return eigenvalues;
    }

    /**
     * Computes the eigenvalues and eigenvectors for a symmetric matrix.
     *
     * @return an array of DoubleMatrix objects containing the eigenvectors
     *         stored as the columns of the first matrix, and the eigenvalues as
     *         diagonal elements of the second matrix.
     */
    public static @Dimensionless DoubleMatrix @Dimensionless [] symmetricEigenvectors(@Dimensionless DoubleMatrix A) {
        A.assertSquare();
        @Dimensionless
        DoubleMatrix eigenvalues = new @Dimensionless DoubleMatrix(A.rows);
        @Dimensionless
        DoubleMatrix eigenvectors = A.dup();
        @Dimensionless
        int isuppz @Dimensionless [] = new int @Dimensionless [((@Dimensionless int) (2)) * A.rows];
        SimpleBlas.syevr('V', 'A', 'U', A.dup(), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), eigenvalues, eigenvectors, isuppz);
        return new DoubleMatrix @Dimensionless []{eigenvectors, DoubleMatrix.diag(eigenvalues)};
    }

    /**
     * Computes the eigenvalues of a general matrix.
     */
    public static @Dimensionless ComplexDoubleMatrix eigenvalues(@Dimensionless DoubleMatrix A) {
        A.assertSquare();
        @Dimensionless
        DoubleMatrix WR = new @Dimensionless DoubleMatrix(A.rows);
        @Dimensionless
        DoubleMatrix WI = WR.dup();
        SimpleBlas.geev('N', 'N', A.dup(), WR, WI, dummyDouble, dummyDouble);

        return new @Dimensionless ComplexDoubleMatrix(WR, WI);
    }

    /**
     * Computes the eigenvalues and eigenvectors of a general matrix.
     *
     * @return an array of ComplexDoubleMatrix objects containing the eigenvectors
     *         stored as the columns of the first matrix, and the eigenvalues as the
     *         diagonal elements of the second matrix.
     */
    public static @Dimensionless ComplexDoubleMatrix @Dimensionless [] eigenvectors(@Dimensionless DoubleMatrix A) {
        A.assertSquare();
        // setting up result arrays
        @Dimensionless
        DoubleMatrix WR = new @Dimensionless DoubleMatrix(A.rows);
        @Dimensionless
        DoubleMatrix WI = WR.dup();
        @Dimensionless
        DoubleMatrix VR = new @Dimensionless DoubleMatrix(A.rows, A.rows);

        SimpleBlas.geev('N', 'V', A.dup(), WR, WI, dummyDouble, VR);

        // transferring the result
        @Dimensionless
        ComplexDoubleMatrix E = new @Dimensionless ComplexDoubleMatrix(WR, WI);
        @Dimensionless
        ComplexDoubleMatrix V = new @Dimensionless ComplexDoubleMatrix(A.rows, A.rows);
        //System.err.printf("VR = %s\n", VR.toString());
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < A.rows; i++) {
            if (E.get(i).isReal()) {
                V.putColumn(i, new @Dimensionless ComplexDoubleMatrix(VR.getColumn(i)));
            } else {
                @Dimensionless
                ComplexDoubleMatrix v = new @Dimensionless ComplexDoubleMatrix(VR.getColumn(i), VR.getColumn(i + ((@Dimensionless int) (1))));
                V.putColumn(i, v);
                V.putColumn(i + ((@Dimensionless int) (1)), v.conji());
                i += ((@Dimensionless int) (1));
            }
        }
        return new ComplexDoubleMatrix @Dimensionless []{V, ComplexDoubleMatrix.diag(E)};
    }

    /**
     * Compute generalized eigenvalues of the problem A x = L B x.
     *
     * @param A symmetric Matrix A. Only the upper triangle will be considered.
     * @param B symmetric Matrix B. Only the upper triangle will be considered.
     * @return a vector of eigenvalues L.
     */
    public static @Dimensionless DoubleMatrix symmetricGeneralizedEigenvalues(@Dimensionless DoubleMatrix A, @Dimensionless DoubleMatrix B) {
        A.assertSquare();
        B.assertSquare();
        @Dimensionless
        DoubleMatrix W = new @Dimensionless DoubleMatrix(A.rows);
        SimpleBlas.sygvd(((@Dimensionless int) (1)), 'N', 'U', A.dup(), B.dup(), W);
        return W;
    }

    /**
     * Solve a general problem A x = L B x.
     *
     * @param A symmetric matrix A
     * @param B symmetric matrix B
     * @return an array of matrices of length two. The first one is an array of the eigenvectors X
     *         The second one is A vector containing the corresponding eigenvalues L.
     */
    public static @Dimensionless DoubleMatrix @Dimensionless [] symmetricGeneralizedEigenvectors(@Dimensionless DoubleMatrix A, @Dimensionless DoubleMatrix B) {
        A.assertSquare();
        B.assertSquare();
        @Dimensionless
        DoubleMatrix @Dimensionless [] result = new DoubleMatrix @Dimensionless [((@Dimensionless int) (2))];
        @Dimensionless
        DoubleMatrix dA = A.dup();
        @Dimensionless
        DoubleMatrix dB = B.dup();
        @Dimensionless
        DoubleMatrix W = new @Dimensionless DoubleMatrix(dA.rows);
        SimpleBlas.sygvd(((@Dimensionless int) (1)), 'V', 'U', dA, dB, W);
        result[((@Dimensionless int) (0))] = dA;
        result[((@Dimensionless int) (1))] = W;
        return result;
    }

  /**
   * Computes selected eigenvalues of the real generalized symmetric-definite eigenproblem of the form A x = L B x
   * or, equivalently, (A - L B)x = 0. Here A and B are assumed to be symmetric and B is also positive definite.
   * The selection is based on the given range of values for the desired eigenvalues.
   * <p/>
   * The range is half open: (vl,vu].
   *
   * @param A  symmetric Matrix A. Only the upper triangle will be considered.
   * @param B  symmetric Matrix B. Only the upper triangle will be considered.
   * @param vl lower bound of the smallest eigenvalue to return
   * @param vu upper bound of the largest eigenvalue to return
   * @throws IllegalArgumentException if <code>vl &gt; vu</code>
   * @throws NoEigenResultException   if no eigevalues are found for the selected range: (vl,vu]
   * @return a vector of eigenvalues L
   */
  public static @Dimensionless DoubleMatrix symmetricGeneralizedEigenvalues(@Dimensionless DoubleMatrix A, @Dimensionless DoubleMatrix B, @Dimensionless double vl, @Dimensionless double vu) {
    A.assertSquare();
    B.assertSquare();
    A.assertSameSize(B);
    if (vu <= vl) {
      throw new @Dimensionless IllegalArgumentException("Bound exception: make sure vu > vl");
    }
    @Dimensionless
    double abstol = (@Dimensionless double) ((@Dimensionless double) (1e-9));    // What is a good tolerance?
    @Dimensionless
    int @Dimensionless [] m = new int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    DoubleMatrix W = new @Dimensionless DoubleMatrix(A.rows);
    @Dimensionless
    DoubleMatrix Z = new @Dimensionless DoubleMatrix(A.rows, A.rows);
    SimpleBlas.sygvx(((@Dimensionless int) (1)), 'N', 'V', 'U', A.dup(), B.dup(), vl, vu, ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), abstol, m, W, Z);
    if (m[((@Dimensionless int) (0))] == ((@Dimensionless int) (0))) {
      throw new @Dimensionless NoEigenResultException("No eigenvalues found for selected range");
    }
    return W.get(new @Dimensionless IntervalRange(((@Dimensionless int) (0)), m[((@Dimensionless int) (0))]), ((@Dimensionless int) (0)));
  }

  /**
   * Computes selected eigenvalues of the real generalized symmetric-definite eigenproblem of the form A x = L B x
   * or, equivalently, (A - L B)x = 0. Here A and B are assumed to be symmetric and B is also positive definite.
   * The selection is based on the given range of indices for the desired eigenvalues.
   *
   * @param A  symmetric Matrix A. Only the upper triangle will be considered.
   * @param B  symmetric Matrix B. Only the upper triangle will be considered.
   * @param il lower index (in ascending order) of the smallest eigenvalue to return (index is 0-based)
   * @param iu upper index (in ascending order) of the largest eigenvalue to return (index is 0-based)
   * @throws IllegalArgumentException if <code>il &gt; iu</code> or <code>il &lt; 0 </code> or <code>iu &gt; A.rows - 1</code>
   * @return a vector of eigenvalues L
   */
  public static @Dimensionless DoubleMatrix symmetricGeneralizedEigenvalues(@Dimensionless DoubleMatrix A, @Dimensionless DoubleMatrix B, @Dimensionless int il, @Dimensionless int iu) {
    A.assertSquare();
    B.assertSquare();
    A.assertSameSize(B);
    if (iu < il) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure iu >= il");
    }
    if (il < ((@Dimensionless int) (0))) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure il >= 0");
    }
    if (iu > A.rows - ((@Dimensionless int) (1))) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure iu <= A.rows - 1");
    }
    @Dimensionless
    double abstol = (@Dimensionless double) ((@Dimensionless double) (1e-9));    // What is a good tolerance?
    @Dimensionless
    int @Dimensionless [] m = new int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    DoubleMatrix W = new @Dimensionless DoubleMatrix(A.rows);
    @Dimensionless
    DoubleMatrix Z = new @Dimensionless DoubleMatrix(A.rows, A.columns);
    SimpleBlas.sygvx(((@Dimensionless int) (1)), 'N', 'I', 'U', A.dup(), B.dup(), ((@Dimensionless double) (0.0)), ((@Dimensionless double) (0.0)), il + ((@Dimensionless int) (1)), iu + ((@Dimensionless int) (1)), abstol, m, W, Z);
    return W.get(new @Dimensionless IntervalRange(((@Dimensionless int) (0)), m[((@Dimensionless int) (0))]), ((@Dimensionless int) (0)));
  }

  /**
   * Computes selected eigenvalues and their corresponding eigenvectors of the real generalized symmetric-definite
   * eigenproblem of the form A x = L B x or, equivalently, (A - L B)x = 0. Here A and B are assumed to be symmetric
   * and B is also positive definite. The selection is based on the given range of values for the desired eigenvalues.
   *
   * The range is half open: (vl,vu].
   *
   * @param A  symmetric Matrix A. Only the upper triangle will be considered.
   * @param B  symmetric Matrix B. Only the upper triangle will be considered.
   * @param vl lower bound of the smallest eigenvalue to return
   * @param vu upper bound of the largest eigenvalue to return
   * @return an array of matrices of length two. The first one is an array of the eigenvectors x.
   *         The second one is a vector containing the corresponding eigenvalues L.
   * @throws IllegalArgumentException if <code>vl &gt; vu</code>
   * @throws NoEigenResultException   if no eigevalues are found for the selected range: (vl,vu]
   */
  public static @Dimensionless DoubleMatrix @Dimensionless [] symmetricGeneralizedEigenvectors(@Dimensionless DoubleMatrix A, @Dimensionless DoubleMatrix B, @Dimensionless double vl, @Dimensionless double vu) {
    A.assertSquare();
    B.assertSquare();
    A.assertSameSize(B);
    if (vu <= vl) {
      throw new @Dimensionless IllegalArgumentException("Bound exception: make sure vu > vl");
    }
    @Dimensionless
    double abstol = (@Dimensionless double) ((@Dimensionless double) (1e-9));    // What is a good tolerance?
    @Dimensionless
    int @Dimensionless [] m = new int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    DoubleMatrix W = new @Dimensionless DoubleMatrix(A.rows);
    @Dimensionless
    DoubleMatrix Z = new @Dimensionless DoubleMatrix(A.rows, A.columns);
    SimpleBlas.sygvx(((@Dimensionless int) (1)), 'V', 'V', 'U', A.dup(), B.dup(), vl, vu, ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), abstol, m, W, Z);
    if (m[((@Dimensionless int) (0))] == ((@Dimensionless int) (0))) {
      throw new @Dimensionless NoEigenResultException("No eigenvalues found for selected range");
    }
    @Dimensionless
    DoubleMatrix @Dimensionless [] result = new DoubleMatrix @Dimensionless [((@Dimensionless int) (2))];
    @Dimensionless
    Range r = new @Dimensionless IntervalRange(((@Dimensionless int) (0)), m[((@Dimensionless int) (0))]);
    result[((@Dimensionless int) (0))] = Z.get(new @Dimensionless IntervalRange(((@Dimensionless int) (0)), A.rows), r);
    result[((@Dimensionless int) (1))] = W.get(r, ((@Dimensionless int) (0)));
    return result;
  }

  /**
   * Computes selected eigenvalues and their corresponding
   * eigenvectors of the real generalized symmetric-definite
   * eigenproblem of the form A x = L B x or, equivalently,
   * (A - L B)x = 0. Here A and B are assumed to be symmetric
   * and B is also positive definite. The selection is based
   * on the given range of values for the desired eigenvalues.
   *
   * @param A  symmetric Matrix A. Only the upper triangle will be considered.
   * @param B  symmetric Matrix B. Only the upper triangle will be considered.
   * @param il lower index (in ascending order) of the smallest eigenvalue to return (index is 0-based)
   * @param iu upper index (in ascending order) of the largest eigenvalue to return (index is 0-based)
   * @throws IllegalArgumentException if <code>il &gt; iu</code> or <code>il &lt; 0 </code>
   *         or <code>iu &gt; A.rows - 1</code>
   * @return an array of matrices of length two. The first one is an array of the eigenvectors x.
   * The second one is a vector containing the corresponding eigenvalues L.
   */
  public static @Dimensionless DoubleMatrix @Dimensionless [] symmetricGeneralizedEigenvectors(@Dimensionless DoubleMatrix A, @Dimensionless DoubleMatrix B, @Dimensionless int il, @Dimensionless int iu) {
    A.assertSquare();
    B.assertSquare();
    A.assertSameSize(B);
    if (iu < il) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure iu >= il");
    }
    if (il < ((@Dimensionless int) (0))) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure il >= 0");
    }
    if (iu > A.rows - ((@Dimensionless int) (1))) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure iu <= A.rows - 1");
    }
    @Dimensionless
    double abstol = (@Dimensionless double) ((@Dimensionless double) (1e-9));    // What is a good tolerance?
    @Dimensionless
    int @Dimensionless [] m = new int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    DoubleMatrix W = new @Dimensionless DoubleMatrix(A.rows);
    @Dimensionless
    DoubleMatrix Z = new @Dimensionless DoubleMatrix(A.rows, A.columns);
    SimpleBlas.sygvx(((@Dimensionless int) (1)), 'V', 'I', 'U', A.dup(), B.dup(), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), il + ((@Dimensionless int) (1)), iu + ((@Dimensionless int) (1)), abstol, m, W, Z);
    @Dimensionless
    DoubleMatrix @Dimensionless [] result = new DoubleMatrix @Dimensionless [((@Dimensionless int) (2))];
    @Dimensionless
    Range r = new @Dimensionless IntervalRange(((@Dimensionless int) (0)), m[((@Dimensionless int) (0))]);
    result[((@Dimensionless int) (0))] = Z.get(new @Dimensionless IntervalRange(((@Dimensionless int) (0)), A.rows), r);
    result[((@Dimensionless int) (1))] = W.get(r, ((@Dimensionless int) (0)));
    return result;
  }

//BEGIN
  // The code below has been automatically generated.
  // DO NOT EDIT!
    private static final @Dimensionless FloatMatrix dummyFloat = new @Dimensionless FloatMatrix(((@Dimensionless int) (1)));

    /**
     * Compute the eigenvalues for a symmetric matrix.
     */
    public static @Dimensionless FloatMatrix symmetricEigenvalues(@Dimensionless FloatMatrix A) {
        A.assertSquare();
        @Dimensionless
        FloatMatrix eigenvalues = new @Dimensionless FloatMatrix(A.rows);
        @Dimensionless
        int isuppz @Dimensionless [] = new int @Dimensionless [((@Dimensionless int) (2)) * A.rows];
        SimpleBlas.syevr('N', 'A', 'U', A.dup(), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), eigenvalues, dummyFloat, isuppz);
        return eigenvalues;
    }

    /**
     * Computes the eigenvalues and eigenvectors for a symmetric matrix.
     *
     * @return an array of FloatMatrix objects containing the eigenvectors
     *         stored as the columns of the first matrix, and the eigenvalues as
     *         diagonal elements of the second matrix.
     */
    public static @Dimensionless FloatMatrix @Dimensionless [] symmetricEigenvectors(@Dimensionless FloatMatrix A) {
        A.assertSquare();
        @Dimensionless
        FloatMatrix eigenvalues = new @Dimensionless FloatMatrix(A.rows);
        @Dimensionless
        FloatMatrix eigenvectors = A.dup();
        @Dimensionless
        int isuppz @Dimensionless [] = new int @Dimensionless [((@Dimensionless int) (2)) * A.rows];
        SimpleBlas.syevr('V', 'A', 'U', A.dup(), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), eigenvalues, eigenvectors, isuppz);
        return new FloatMatrix @Dimensionless []{eigenvectors, FloatMatrix.diag(eigenvalues)};
    }

    /**
     * Computes the eigenvalues of a general matrix.
     */
    public static @Dimensionless ComplexFloatMatrix eigenvalues(@Dimensionless FloatMatrix A) {
        A.assertSquare();
        @Dimensionless
        FloatMatrix WR = new @Dimensionless FloatMatrix(A.rows);
        @Dimensionless
        FloatMatrix WI = WR.dup();
        SimpleBlas.geev('N', 'N', A.dup(), WR, WI, dummyFloat, dummyFloat);

        return new @Dimensionless ComplexFloatMatrix(WR, WI);
    }

    /**
     * Computes the eigenvalues and eigenvectors of a general matrix.
     *
     * @return an array of ComplexFloatMatrix objects containing the eigenvectors
     *         stored as the columns of the first matrix, and the eigenvalues as the
     *         diagonal elements of the second matrix.
     */
    public static @Dimensionless ComplexFloatMatrix @Dimensionless [] eigenvectors(@Dimensionless FloatMatrix A) {
        A.assertSquare();
        // setting up result arrays
        @Dimensionless
        FloatMatrix WR = new @Dimensionless FloatMatrix(A.rows);
        @Dimensionless
        FloatMatrix WI = WR.dup();
        @Dimensionless
        FloatMatrix VR = new @Dimensionless FloatMatrix(A.rows, A.rows);

        SimpleBlas.geev('N', 'V', A.dup(), WR, WI, dummyFloat, VR);

        // transferring the result
        @Dimensionless
        ComplexFloatMatrix E = new @Dimensionless ComplexFloatMatrix(WR, WI);
        @Dimensionless
        ComplexFloatMatrix V = new @Dimensionless ComplexFloatMatrix(A.rows, A.rows);
        //System.err.printf("VR = %s\n", VR.toString());
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < A.rows; i++) {
            if (E.get(i).isReal()) {
                V.putColumn(i, new @Dimensionless ComplexFloatMatrix(VR.getColumn(i)));
            } else {
                @Dimensionless
                ComplexFloatMatrix v = new @Dimensionless ComplexFloatMatrix(VR.getColumn(i), VR.getColumn(i + ((@Dimensionless int) (1))));
                V.putColumn(i, v);
                V.putColumn(i + ((@Dimensionless int) (1)), v.conji());
                i += ((@Dimensionless int) (1));
            }
        }
        return new ComplexFloatMatrix @Dimensionless []{V, ComplexFloatMatrix.diag(E)};
    }

    /**
     * Compute generalized eigenvalues of the problem A x = L B x.
     *
     * @param A symmetric Matrix A. Only the upper triangle will be considered.
     * @param B symmetric Matrix B. Only the upper triangle will be considered.
     * @return a vector of eigenvalues L.
     */
    public static @Dimensionless FloatMatrix symmetricGeneralizedEigenvalues(@Dimensionless FloatMatrix A, @Dimensionless FloatMatrix B) {
        A.assertSquare();
        B.assertSquare();
        @Dimensionless
        FloatMatrix W = new @Dimensionless FloatMatrix(A.rows);
        SimpleBlas.sygvd(((@Dimensionless int) (1)), 'N', 'U', A.dup(), B.dup(), W);
        return W;
    }

    /**
     * Solve a general problem A x = L B x.
     *
     * @param A symmetric matrix A
     * @param B symmetric matrix B
     * @return an array of matrices of length two. The first one is an array of the eigenvectors X
     *         The second one is A vector containing the corresponding eigenvalues L.
     */
    public static @Dimensionless FloatMatrix @Dimensionless [] symmetricGeneralizedEigenvectors(@Dimensionless FloatMatrix A, @Dimensionless FloatMatrix B) {
        A.assertSquare();
        B.assertSquare();
        @Dimensionless
        FloatMatrix @Dimensionless [] result = new FloatMatrix @Dimensionless [((@Dimensionless int) (2))];
        @Dimensionless
        FloatMatrix dA = A.dup();
        @Dimensionless
        FloatMatrix dB = B.dup();
        @Dimensionless
        FloatMatrix W = new @Dimensionless FloatMatrix(dA.rows);
        SimpleBlas.sygvd(((@Dimensionless int) (1)), 'V', 'U', dA, dB, W);
        result[((@Dimensionless int) (0))] = dA;
        result[((@Dimensionless int) (1))] = W;
        return result;
    }

  /**
   * Computes selected eigenvalues of the real generalized symmetric-definite eigenproblem of the form A x = L B x
   * or, equivalently, (A - L B)x = 0. Here A and B are assumed to be symmetric and B is also positive definite.
   * The selection is based on the given range of values for the desired eigenvalues.
   * <p/>
   * The range is half open: (vl,vu].
   *
   * @param A  symmetric Matrix A. Only the upper triangle will be considered.
   * @param B  symmetric Matrix B. Only the upper triangle will be considered.
   * @param vl lower bound of the smallest eigenvalue to return
   * @param vu upper bound of the largest eigenvalue to return
   * @throws IllegalArgumentException if <code>vl &gt; vu</code>
   * @throws NoEigenResultException   if no eigevalues are found for the selected range: (vl,vu]
   * @return a vector of eigenvalues L
   */
  public static @Dimensionless FloatMatrix symmetricGeneralizedEigenvalues(@Dimensionless FloatMatrix A, @Dimensionless FloatMatrix B, @Dimensionless float vl, @Dimensionless float vu) {
    A.assertSquare();
    B.assertSquare();
    A.assertSameSize(B);
    if (vu <= vl) {
      throw new @Dimensionless IllegalArgumentException("Bound exception: make sure vu > vl");
    }
    @Dimensionless
    float abstol = (@Dimensionless float) ((@Dimensionless double) (1e-9));    // What is a good tolerance?
    @Dimensionless
    int @Dimensionless [] m = new int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    FloatMatrix W = new @Dimensionless FloatMatrix(A.rows);
    @Dimensionless
    FloatMatrix Z = new @Dimensionless FloatMatrix(A.rows, A.rows);
    SimpleBlas.sygvx(((@Dimensionless int) (1)), 'N', 'V', 'U', A.dup(), B.dup(), vl, vu, ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), abstol, m, W, Z);
    if (m[((@Dimensionless int) (0))] == ((@Dimensionless int) (0))) {
      throw new @Dimensionless NoEigenResultException("No eigenvalues found for selected range");
    }
    return W.get(new @Dimensionless IntervalRange(((@Dimensionless int) (0)), m[((@Dimensionless int) (0))]), ((@Dimensionless int) (0)));
  }

  /**
   * Computes selected eigenvalues of the real generalized symmetric-definite eigenproblem of the form A x = L B x
   * or, equivalently, (A - L B)x = 0. Here A and B are assumed to be symmetric and B is also positive definite.
   * The selection is based on the given range of indices for the desired eigenvalues.
   *
   * @param A  symmetric Matrix A. Only the upper triangle will be considered.
   * @param B  symmetric Matrix B. Only the upper triangle will be considered.
   * @param il lower index (in ascending order) of the smallest eigenvalue to return (index is 0-based)
   * @param iu upper index (in ascending order) of the largest eigenvalue to return (index is 0-based)
   * @throws IllegalArgumentException if <code>il &gt; iu</code> or <code>il &lt; 0 </code> or <code>iu &gt; A.rows - 1</code>
   * @return a vector of eigenvalues L
   */
  public static @Dimensionless FloatMatrix symmetricGeneralizedEigenvalues(@Dimensionless FloatMatrix A, @Dimensionless FloatMatrix B, @Dimensionless int il, @Dimensionless int iu) {
    A.assertSquare();
    B.assertSquare();
    A.assertSameSize(B);
    if (iu < il) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure iu >= il");
    }
    if (il < ((@Dimensionless int) (0))) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure il >= 0");
    }
    if (iu > A.rows - ((@Dimensionless int) (1))) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure iu <= A.rows - 1");
    }
    @Dimensionless
    float abstol = (@Dimensionless float) ((@Dimensionless double) (1e-9));    // What is a good tolerance?
    @Dimensionless
    int @Dimensionless [] m = new int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    FloatMatrix W = new @Dimensionless FloatMatrix(A.rows);
    @Dimensionless
    FloatMatrix Z = new @Dimensionless FloatMatrix(A.rows, A.columns);
    SimpleBlas.sygvx(((@Dimensionless int) (1)), 'N', 'I', 'U', A.dup(), B.dup(), ((@Dimensionless float) (0.0f)), ((@Dimensionless float) (0.0f)), il + ((@Dimensionless int) (1)), iu + ((@Dimensionless int) (1)), abstol, m, W, Z);
    return W.get(new @Dimensionless IntervalRange(((@Dimensionless int) (0)), m[((@Dimensionless int) (0))]), ((@Dimensionless int) (0)));
  }

  /**
   * Computes selected eigenvalues and their corresponding eigenvectors of the real generalized symmetric-definite
   * eigenproblem of the form A x = L B x or, equivalently, (A - L B)x = 0. Here A and B are assumed to be symmetric
   * and B is also positive definite. The selection is based on the given range of values for the desired eigenvalues.
   *
   * The range is half open: (vl,vu].
   *
   * @param A  symmetric Matrix A. Only the upper triangle will be considered.
   * @param B  symmetric Matrix B. Only the upper triangle will be considered.
   * @param vl lower bound of the smallest eigenvalue to return
   * @param vu upper bound of the largest eigenvalue to return
   * @return an array of matrices of length two. The first one is an array of the eigenvectors x.
   *         The second one is a vector containing the corresponding eigenvalues L.
   * @throws IllegalArgumentException if <code>vl &gt; vu</code>
   * @throws NoEigenResultException   if no eigevalues are found for the selected range: (vl,vu]
   */
  public static @Dimensionless FloatMatrix @Dimensionless [] symmetricGeneralizedEigenvectors(@Dimensionless FloatMatrix A, @Dimensionless FloatMatrix B, @Dimensionless float vl, @Dimensionless float vu) {
    A.assertSquare();
    B.assertSquare();
    A.assertSameSize(B);
    if (vu <= vl) {
      throw new @Dimensionless IllegalArgumentException("Bound exception: make sure vu > vl");
    }
    @Dimensionless
    float abstol = (@Dimensionless float) ((@Dimensionless double) (1e-9));    // What is a good tolerance?
    @Dimensionless
    int @Dimensionless [] m = new int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    FloatMatrix W = new @Dimensionless FloatMatrix(A.rows);
    @Dimensionless
    FloatMatrix Z = new @Dimensionless FloatMatrix(A.rows, A.columns);
    SimpleBlas.sygvx(((@Dimensionless int) (1)), 'V', 'V', 'U', A.dup(), B.dup(), vl, vu, ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), abstol, m, W, Z);
    if (m[((@Dimensionless int) (0))] == ((@Dimensionless int) (0))) {
      throw new @Dimensionless NoEigenResultException("No eigenvalues found for selected range");
    }
    @Dimensionless
    FloatMatrix @Dimensionless [] result = new FloatMatrix @Dimensionless [((@Dimensionless int) (2))];
    @Dimensionless
    Range r = new @Dimensionless IntervalRange(((@Dimensionless int) (0)), m[((@Dimensionless int) (0))]);
    result[((@Dimensionless int) (0))] = Z.get(new @Dimensionless IntervalRange(((@Dimensionless int) (0)), A.rows), r);
    result[((@Dimensionless int) (1))] = W.get(r, ((@Dimensionless int) (0)));
    return result;
  }

  /**
   * Computes selected eigenvalues and their corresponding
   * eigenvectors of the real generalized symmetric-definite
   * eigenproblem of the form A x = L B x or, equivalently,
   * (A - L B)x = 0. Here A and B are assumed to be symmetric
   * and B is also positive definite. The selection is based
   * on the given range of values for the desired eigenvalues.
   *
   * @param A  symmetric Matrix A. Only the upper triangle will be considered.
   * @param B  symmetric Matrix B. Only the upper triangle will be considered.
   * @param il lower index (in ascending order) of the smallest eigenvalue to return (index is 0-based)
   * @param iu upper index (in ascending order) of the largest eigenvalue to return (index is 0-based)
   * @throws IllegalArgumentException if <code>il &gt; iu</code> or <code>il &lt; 0 </code>
   *         or <code>iu &gt; A.rows - 1</code>
   * @return an array of matrices of length two. The first one is an array of the eigenvectors x.
   * The second one is a vector containing the corresponding eigenvalues L.
   */
  public static @Dimensionless FloatMatrix @Dimensionless [] symmetricGeneralizedEigenvectors(@Dimensionless FloatMatrix A, @Dimensionless FloatMatrix B, @Dimensionless int il, @Dimensionless int iu) {
    A.assertSquare();
    B.assertSquare();
    A.assertSameSize(B);
    if (iu < il) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure iu >= il");
    }
    if (il < ((@Dimensionless int) (0))) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure il >= 0");
    }
    if (iu > A.rows - ((@Dimensionless int) (1))) {
      throw new @Dimensionless IllegalArgumentException("Index exception: make sure iu <= A.rows - 1");
    }
    @Dimensionless
    float abstol = (@Dimensionless float) ((@Dimensionless double) (1e-9));    // What is a good tolerance?
    @Dimensionless
    int @Dimensionless [] m = new int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    FloatMatrix W = new @Dimensionless FloatMatrix(A.rows);
    @Dimensionless
    FloatMatrix Z = new @Dimensionless FloatMatrix(A.rows, A.columns);
    SimpleBlas.sygvx(((@Dimensionless int) (1)), 'V', 'I', 'U', A.dup(), B.dup(), ((@Dimensionless int) (0)), ((@Dimensionless int) (0)), il + ((@Dimensionless int) (1)), iu + ((@Dimensionless int) (1)), abstol, m, W, Z);
    @Dimensionless
    FloatMatrix @Dimensionless [] result = new FloatMatrix @Dimensionless [((@Dimensionless int) (2))];
    @Dimensionless
    Range r = new @Dimensionless IntervalRange(((@Dimensionless int) (0)), m[((@Dimensionless int) (0))]);
    result[((@Dimensionless int) (0))] = Z.get(new @Dimensionless IntervalRange(((@Dimensionless int) (0)), A.rows), r);
    result[((@Dimensionless int) (1))] = W.get(r, ((@Dimensionless int) (0)));
    return result;
  }

//END
}
