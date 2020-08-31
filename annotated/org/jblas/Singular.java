/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.jblas;

import units.qual.Dimensionless;
import org.jblas.exceptions.LapackConvergenceException;

import static org.jblas.util.Functions.min;

/**
 *
 */
public class Singular {

    /**
     * Compute a singular-value decomposition of A.
     *
     * @return A DoubleMatrix[3] array of U, S, V such that A = U * diag(S) * V'
     */
    public static @Dimensionless DoubleMatrix @Dimensionless [] fullSVD(@Dimensionless DoubleMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;

        @Dimensionless
        DoubleMatrix U = new @Dimensionless DoubleMatrix(m, m);
        @Dimensionless
        DoubleMatrix S = new @Dimensionless DoubleMatrix(min(m, n));
        @Dimensionless
        DoubleMatrix V = new @Dimensionless DoubleMatrix(n, n);

        @Dimensionless
        int info = NativeBlas.dgesvd('A', 'A', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), U.data, ((@Dimensionless int) (0)), m, V.data, ((@Dimensionless int) (0)), n);

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return new DoubleMatrix @Dimensionless []{U, S, V.transpose()};
    }

    /**
     * Compute a singular-value decomposition of A (sparse variant).
     * Sparse means that the matrices U and V are not square but
     * only have as many columns (or rows) as necessary.
     * 
     * @param A
     * @return A DoubleMatrix[3] array of U, S, V such that A = U * diag(S) * V'
     */
    public static DoubleMatrix[] sparseSVD(DoubleMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;

        @Dimensionless
        DoubleMatrix U = new @Dimensionless DoubleMatrix(m, min(m, n));
        @Dimensionless
        DoubleMatrix S = new @Dimensionless DoubleMatrix(min(m, n));
        @Dimensionless
        DoubleMatrix V = new @Dimensionless DoubleMatrix(min(m, n), n);

        @Dimensionless
        int info = NativeBlas.dgesvd('S', 'S', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), U.data, ((@Dimensionless int) (0)), m, V.data, ((@Dimensionless int) (0)), min(m, n));

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return new DoubleMatrix @Dimensionless []{U, S, V.transpose()};
    }

  /**
   * Compute a singular-value decomposition of A (sparse variant).
   * Sparse means that the matrices U and V are not square but only have
   * as many columns (or rows) as necessary.
   *
   * @param A
   * @return A ComplexDoubleMatrix[3] array of U, S, V such that A = U * diag(S) * V*
   */
    public static @Dimensionless ComplexDoubleMatrix @Dimensionless [] sparseSVD(@Dimensionless ComplexDoubleMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;

        @Dimensionless
        ComplexDoubleMatrix U = new @Dimensionless ComplexDoubleMatrix(m, min(m, n));
        @Dimensionless
        DoubleMatrix S = new @Dimensionless DoubleMatrix(min(m, n));
        @Dimensionless
        ComplexDoubleMatrix V = new @Dimensionless ComplexDoubleMatrix(min(m, n), n);

        @Dimensionless
        double @Dimensionless [] rwork = new double @Dimensionless [((@Dimensionless int) (5))*min(m,n)];

        @Dimensionless
        int info = NativeBlas.zgesvd('S', 'S', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), U.data, ((@Dimensionless int) (0)), m, V.data, ((@Dimensionless int) (0)), min(m, n), rwork, ((@Dimensionless int) (0)));

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return new ComplexDoubleMatrix @Dimensionless []{U, new @Dimensionless ComplexDoubleMatrix(S), V.hermitian()};
    }

    /**
     * Compute a singular-value decomposition of A.
     *
     * @return A ComplexDoubleMatrix[3] array of U, S, V such that A = U * diag(S) * V'
     */
    public static @Dimensionless ComplexDoubleMatrix @Dimensionless [] fullSVD(@Dimensionless ComplexDoubleMatrix A) {
      @Dimensionless
      int m = A.rows;
      @Dimensionless
      int n = A.columns;

      @Dimensionless
      ComplexDoubleMatrix U = new @Dimensionless ComplexDoubleMatrix(m, m);
      @Dimensionless
      DoubleMatrix S = new @Dimensionless DoubleMatrix(min(m, n));
      @Dimensionless
      ComplexDoubleMatrix V = new @Dimensionless ComplexDoubleMatrix(n, n);

      @Dimensionless
      double @Dimensionless [] rwork = new double @Dimensionless [((@Dimensionless int) (5))*min(m,n)];

      @Dimensionless
      int info = NativeBlas.zgesvd('A', 'A', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), U.data, ((@Dimensionless int) (0)), m, V.data, ((@Dimensionless int) (0)), n, rwork, ((@Dimensionless int) (0)));

      if (info > ((@Dimensionless int) (0))) {
        throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
      }

      return new ComplexDoubleMatrix @Dimensionless []{U, new @Dimensionless ComplexDoubleMatrix(S), V.hermitian()};
    }

    /**
     * Compute the singular values of a matrix.
     *
     * @param A DoubleMatrix of dimension m * n
     * @return A min(m, n) vector of singular values.
     */
    public static DoubleMatrix SVDValues(DoubleMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;
        @Dimensionless
        DoubleMatrix S = new @Dimensionless DoubleMatrix(min(m, n));

        @Dimensionless
        int info = NativeBlas.dgesvd('N', 'N', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), null, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), null, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return S;
    }

    /**
     * Compute the singular values of a complex matrix.
     *
     * @param A ComplexDoubleMatrix of dimension m * n
     * @return A real-valued (!) min(m, n) vector of singular values.
     */
    public static @Dimensionless DoubleMatrix SVDValues(@Dimensionless ComplexDoubleMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;
        @Dimensionless
        DoubleMatrix S = new @Dimensionless DoubleMatrix(min(m, n));
        @Dimensionless
        double @Dimensionless [] rwork = new double @Dimensionless [((@Dimensionless int) (5))*min(m,n)];

        @Dimensionless
        int info = NativeBlas.zgesvd('N', 'N', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), null, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), null, ((@Dimensionless int) (0)), min(m,n), rwork, ((@Dimensionless int) (0)));

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return S;
    }

    //BEGIN
  // The code below has been automatically generated.
  // DO NOT EDIT!

    /**
     * Compute a singular-value decomposition of A.
     *
     * @return A FloatMatrix[3] array of U, S, V such that A = U * diag(S) * V'
     */
    public static @Dimensionless FloatMatrix @Dimensionless [] fullSVD(@Dimensionless FloatMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;

        @Dimensionless
        FloatMatrix U = new @Dimensionless FloatMatrix(m, m);
        @Dimensionless
        FloatMatrix S = new @Dimensionless FloatMatrix(min(m, n));
        @Dimensionless
        FloatMatrix V = new @Dimensionless FloatMatrix(n, n);

        @Dimensionless
        int info = NativeBlas.sgesvd('A', 'A', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), U.data, ((@Dimensionless int) (0)), m, V.data, ((@Dimensionless int) (0)), n);

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return new FloatMatrix @Dimensionless []{U, S, V.transpose()};
    }

    /**
     * Compute a singular-value decomposition of A (sparse variant).
     * Sparse means that the matrices U and V are not square but
     * only have as many columns (or rows) as necessary.
     * 
     * @param A
     * @return A FloatMatrix[3] array of U, S, V such that A = U * diag(S) * V'
     */
    public static @Dimensionless FloatMatrix @Dimensionless [] sparseSVD(@Dimensionless FloatMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;

        @Dimensionless
        FloatMatrix U = new @Dimensionless FloatMatrix(m, min(m, n));
        @Dimensionless
        FloatMatrix S = new @Dimensionless FloatMatrix(min(m, n));
        @Dimensionless
        FloatMatrix V = new @Dimensionless FloatMatrix(min(m, n), n);

        @Dimensionless
        int info = NativeBlas.sgesvd('S', 'S', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), U.data, ((@Dimensionless int) (0)), m, V.data, ((@Dimensionless int) (0)), min(m, n));

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return new FloatMatrix @Dimensionless []{U, S, V.transpose()};
    }

  /**
   * Compute a singular-value decomposition of A (sparse variant).
   * Sparse means that the matrices U and V are not square but only have
   * as many columns (or rows) as necessary.
   *
   * @param A
   * @return A ComplexFloatMatrix[3] array of U, S, V such that A = U * diag(S) * V*
   */
    public static @Dimensionless ComplexFloatMatrix @Dimensionless [] sparseSVD(@Dimensionless ComplexFloatMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;

        @Dimensionless
        ComplexFloatMatrix U = new @Dimensionless ComplexFloatMatrix(m, min(m, n));
        @Dimensionless
        FloatMatrix S = new @Dimensionless FloatMatrix(min(m, n));
        @Dimensionless
        ComplexFloatMatrix V = new @Dimensionless ComplexFloatMatrix(min(m, n), n);

        @Dimensionless
        float @Dimensionless [] rwork = new float @Dimensionless [((@Dimensionless int) (5))*min(m,n)];

        @Dimensionless
        int info = NativeBlas.cgesvd('S', 'S', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), U.data, ((@Dimensionless int) (0)), m, V.data, ((@Dimensionless int) (0)), min(m, n), rwork, ((@Dimensionless int) (0)));

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return new ComplexFloatMatrix @Dimensionless []{U, new @Dimensionless ComplexFloatMatrix(S), V.hermitian()};
    }

    /**
     * Compute a singular-value decomposition of A.
     *
     * @return A ComplexFloatMatrix[3] array of U, S, V such that A = U * diag(S) * V'
     */
    public static @Dimensionless ComplexFloatMatrix @Dimensionless [] fullSVD(@Dimensionless ComplexFloatMatrix A) {
      @Dimensionless
      int m = A.rows;
      @Dimensionless
      int n = A.columns;

      @Dimensionless
      ComplexFloatMatrix U = new @Dimensionless ComplexFloatMatrix(m, m);
      @Dimensionless
      FloatMatrix S = new @Dimensionless FloatMatrix(min(m, n));
      @Dimensionless
      ComplexFloatMatrix V = new @Dimensionless ComplexFloatMatrix(n, n);

      @Dimensionless
      float @Dimensionless [] rwork = new float @Dimensionless [((@Dimensionless int) (5))*min(m,n)];

      @Dimensionless
      int info = NativeBlas.cgesvd('A', 'A', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), U.data, ((@Dimensionless int) (0)), m, V.data, ((@Dimensionless int) (0)), n, rwork, ((@Dimensionless int) (0)));

      if (info > ((@Dimensionless int) (0))) {
        throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
      }

      return new ComplexFloatMatrix @Dimensionless []{U, new @Dimensionless ComplexFloatMatrix(S), V.hermitian()};
    }

    /**
     * Compute the singular values of a matrix.
     *
     * @param A FloatMatrix of dimension m * n
     * @return A min(m, n) vector of singular values.
     */
    public static @Dimensionless FloatMatrix SVDValues(@Dimensionless FloatMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;
        @Dimensionless
        FloatMatrix S = new @Dimensionless FloatMatrix(min(m, n));

        @Dimensionless
        int info = NativeBlas.sgesvd('N', 'N', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), null, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), null, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return S;
    }

    /**
     * Compute the singular values of a complex matrix.
     *
     * @param A ComplexFloatMatrix of dimension m * n
     * @return A real-valued (!) min(m, n) vector of singular values.
     */
    public static @Dimensionless FloatMatrix SVDValues(@Dimensionless ComplexFloatMatrix A) {
        @Dimensionless
        int m = A.rows;
        @Dimensionless
        int n = A.columns;
        @Dimensionless
        FloatMatrix S = new @Dimensionless FloatMatrix(min(m, n));
        @Dimensionless
        float @Dimensionless [] rwork = new float @Dimensionless [((@Dimensionless int) (5))*min(m,n)];

        @Dimensionless
        int info = NativeBlas.cgesvd('N', 'N', m, n, A.dup().data, ((@Dimensionless int) (0)), m, S.data, ((@Dimensionless int) (0)), null, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), null, ((@Dimensionless int) (0)), min(m,n), rwork, ((@Dimensionless int) (0)));

        if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackConvergenceException("GESVD", info + " superdiagonals of an intermediate bidiagonal form failed to converge.");
        }

        return S;
    }

    //END
}
