/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.jblas;

import units.qual.Dimensionless;
import org.jblas.exceptions.LapackArgumentException;
import org.jblas.exceptions.LapackPositivityException;
import org.jblas.util.Permutations;
import static org.jblas.util.Functions.min;

/**
 * Matrix which collects all kinds of decompositions.
 */
@Dimensionless
public class Decompose {
    /**
     * Class to hold an LU decomposition result.
     *
     * Contains a lower matrix L, and upper matrix U, and a permutation matrix
     * P such that P*L*U is the original matrix.
     * @param <T>
     */
    @Dimensionless
    public static class LUDecomposition<@Dimensionless T extends @Dimensionless Object> {

        public @Dimensionless T l;
        public @Dimensionless T u;
        public @Dimensionless T p;

        public LUDecomposition(Decompose.@Dimensionless LUDecomposition<T> this, @Dimensionless T l, @Dimensionless T u, @Dimensionless T p) {
            this.l = l;
            this.u = u;
            this.p = p;
        }

        @Override
        public @Dimensionless String toString(Decompose.@Dimensionless LUDecomposition<@Dimensionless T> this) {
          return String.format("<LUDecomposition L=%s U=%s P=%s>", l, u, p);
        }
    }

    /**
     * Compute LU Decomposition of a general matrix.
     *
     * Computes the LU decomposition using GETRF. Returns three matrices L, U, P,
     * where L is lower diagonal, U is upper diagonal, and P is a permutation
     * matrix such that A = P * L * U.
     *
     * @param A general matrix
     * @return An LUDecomposition object.
     */
    public static @Dimensionless LUDecomposition<@Dimensionless DoubleMatrix> lu(@Dimensionless DoubleMatrix A) {
        @Dimensionless
        int @Dimensionless [] ipiv = new int @Dimensionless [min(A.rows, A.columns)];
        @Dimensionless
        DoubleMatrix result = A.dup();
        NativeBlas.dgetrf(A.rows, A.columns, result.data, ((@Dimensionless int) (0)), A.rows, ipiv, ((@Dimensionless int) (0)));

        // collect result
        @Dimensionless
        DoubleMatrix l = new @Dimensionless DoubleMatrix(A.rows, min(A.rows, A.columns));
        @Dimensionless
        DoubleMatrix u = new @Dimensionless DoubleMatrix(min(A.columns, A.rows), A.columns);
        decomposeLowerUpper(result, l, u);
        @Dimensionless
        DoubleMatrix p = Permutations.permutationDoubleMatrixFromPivotIndices(A.rows, ipiv);
        return new @Dimensionless LUDecomposition<@Dimensionless DoubleMatrix>(l, u, p);
    }

    private static void decomposeLowerUpper(@Dimensionless DoubleMatrix A, @Dimensionless DoubleMatrix L, @Dimensionless DoubleMatrix U) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < A.rows; i++) {
            for (@Dimensionless int j = ((@Dimensionless int) (0)); j < A.columns; j++) {
                if (i < j) {
                    U.put(i, j, A.get(i, j));
                } else if (i == j) {
                    U.put(i, i, A.get(i, i));
                    L.put(i, i, ((@Dimensionless double) (1.0)));
                } else {
                    L.put(i, j, A.get(i, j));
                }

            }
        }
    }

    /**if (info )
     * Compute Cholesky decomposition of A
     *
     * @param A symmetric, positive definite matrix (only upper half is used)
     * @return upper triangular matrix U such that  A = U' * U
     */
    public static @Dimensionless FloatMatrix cholesky(@Dimensionless FloatMatrix A) {
        @Dimensionless
        FloatMatrix result = A.dup();
        @Dimensionless
        int info = NativeBlas.spotrf('U', A.rows, result.data, ((@Dimensionless int) (0)), A.rows);
        if (info < ((@Dimensionless int) (0))) {
            throw new @Dimensionless LapackArgumentException("DPOTRF", -info);
        } else if (info > ((@Dimensionless int) (0))) {
            throw new @Dimensionless LapackPositivityException("DPOTRF", "Minor " + info + " was negative. Matrix must be positive definite.");
        }
        clearLower(result);
        return result;
    }

    private static void clearLower(@Dimensionless FloatMatrix A) {
        for (@Dimensionless int j = ((@Dimensionless int) (0)); j < A.columns; j++)
            for (@Dimensionless int i = j + ((@Dimensionless int) (1)); i < A.rows; i++)
                A.put(i, j, ((@Dimensionless float) (0.0f)));
    }

  /**
   * Compute LU Decomposition of a general matrix.
   *
   * Computes the LU decomposition using GETRF. Returns three matrices L, U, P,
   * where L is lower diagonal, U is upper diagonal, and P is a permutation
   * matrix such that A = P * L * U.
   *
   * @param A general matrix
   * @return An LUDecomposition object.
   */
  public static @Dimensionless LUDecomposition<@Dimensionless FloatMatrix> lu(@Dimensionless FloatMatrix A) {
      @Dimensionless
      int @Dimensionless [] ipiv = new int @Dimensionless [min(A.rows, A.columns)];
      @Dimensionless
      FloatMatrix result = A.dup();
      NativeBlas.sgetrf(A.rows, A.columns, result.data, ((@Dimensionless int) (0)), A.rows, ipiv, ((@Dimensionless int) (0)));

      // collect result
      @Dimensionless
      FloatMatrix l = new @Dimensionless FloatMatrix(A.rows, min(A.rows, A.columns));
      @Dimensionless
      FloatMatrix u = new @Dimensionless FloatMatrix(min(A.columns, A.rows), A.columns);
      decomposeLowerUpper(result, l, u);
      @Dimensionless
      FloatMatrix p = Permutations.permutationFloatMatrixFromPivotIndices(A.rows, ipiv);
      return new @Dimensionless LUDecomposition<@Dimensionless FloatMatrix>(l, u, p);
  }

  private static void decomposeLowerUpper(@Dimensionless FloatMatrix A, @Dimensionless FloatMatrix L, @Dimensionless FloatMatrix U) {
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < A.rows; i++) {
          for (@Dimensionless int j = ((@Dimensionless int) (0)); j < A.columns; j++) {
              if (i < j) {
                  U.put(i, j, A.get(i, j));
              } else if (i == j) {
                  U.put(i, i, A.get(i, i));
                  L.put(i, i, ((@Dimensionless float) (1.0f)));
              } else {
                  L.put(i, j, A.get(i, j));
              }

          }
      }
  }

  /**
   * Compute Cholesky decomposition of A
   *
   * @param A symmetric, positive definite matrix (only upper half is used)
   * @return upper triangular matrix U such that  A = U' * U
   */
  public static @Dimensionless DoubleMatrix cholesky(@Dimensionless DoubleMatrix A) {
      @Dimensionless
      DoubleMatrix result = A.dup();
      @Dimensionless
      int info = NativeBlas.dpotrf('U', A.rows, result.data, ((@Dimensionless int) (0)), A.rows);
      if (info < ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackArgumentException("DPOTRF", -info);
      } else if (info > ((@Dimensionless int) (0))) {
          throw new @Dimensionless LapackPositivityException("DPOTRF", "Minor " + info + " was negative. Matrix must be positive definite.");
      }
      clearLower(result);
      return result;
  }

  private static void clearLower(@Dimensionless DoubleMatrix A) {
      for (@Dimensionless int j = ((@Dimensionless int) (0)); j < A.columns; j++)
          for (@Dimensionless int i = j + ((@Dimensionless int) (1)); i < A.rows; i++)
              A.put(i, j, ((@Dimensionless double) (0.0)));
  }

  /**
   * Class to represent a QR decomposition.
   *
   * @param <T>
   */
  @Dimensionless
  public static class QRDecomposition<@Dimensionless T extends @Dimensionless Object> {
    public @Dimensionless T q;
    public @Dimensionless T r;

    QRDecomposition(Decompose.@Dimensionless QRDecomposition<T> this, @Dimensionless T q, @Dimensionless T r) {
      this.q = q;
      this.r = r;
    }

    @Override
    public @Dimensionless String toString(Decompose.@Dimensionless QRDecomposition<@Dimensionless T> this) {
      return "<Q=" + q + " R=" + r + ">";
    }
  }

  /**
   * QR decomposition.
   *
   * Decomposes (m,n) matrix A into a (m,m) matrix Q and an (m,n) matrix R such that
   * Q is orthogonal, R is upper triangular and Q * R = A
   *
   * Note that if A has more rows than columns, then the lower rows of R will contain
   * only zeros, such that the corresponding later columns of Q do not enter the computation
   * at all. For some reason, LAPACK does not properly normalize those columns.
   *
   * @param A matrix
   * @return QR decomposition
   */
  public static @Dimensionless QRDecomposition<@Dimensionless DoubleMatrix> qr(@Dimensionless DoubleMatrix A) {
    @Dimensionless
    int minmn = min(A.rows, A.columns);
    @Dimensionless
    DoubleMatrix result = A.dup();
    @Dimensionless
    DoubleMatrix tau = new @Dimensionless DoubleMatrix(minmn);
    SimpleBlas.geqrf(result, tau);
    @Dimensionless
    DoubleMatrix R = new @Dimensionless DoubleMatrix(A.rows, A.columns);
    for (@Dimensionless int i = ((@Dimensionless int) (0)); i < A.rows; i++) {
      for (@Dimensionless int j = i; j < A.columns; j++) {
        R.put(i, j, result.get(i, j));
      }
    }
    @Dimensionless
    DoubleMatrix Q = DoubleMatrix.eye(A.rows);
    SimpleBlas.ormqr('L', 'N', result, tau, Q);
    return new @Dimensionless QRDecomposition<@Dimensionless DoubleMatrix>(Q, R);
  }
  
  /**
   * QR decomposition.
   *
   * Decomposes (m,n) matrix A into a (m,m) matrix Q and an (m,n) matrix R such that
   * Q is orthogonal, R is upper triangular and Q * R = A
   *
   * @param A matrix
   * @return QR decomposition
   */
  public static @Dimensionless QRDecomposition<@Dimensionless FloatMatrix> qr(@Dimensionless FloatMatrix A) {
    @Dimensionless
    int minmn = min(A.rows, A.columns);
    @Dimensionless
    FloatMatrix result = A.dup();
    @Dimensionless
    FloatMatrix tau = new @Dimensionless FloatMatrix(minmn);
    SimpleBlas.geqrf(result, tau);
    @Dimensionless
    FloatMatrix R = new @Dimensionless FloatMatrix(A.rows, A.columns);
    for (@Dimensionless int i = ((@Dimensionless int) (0)); i < A.rows; i++) {
      for (@Dimensionless int j = i; j < A.columns; j++) {
        R.put(i, j, result.get(i, j));
      }
    }
    @Dimensionless
    FloatMatrix Q = FloatMatrix.eye(A.rows);
    SimpleBlas.ormqr('L', 'N', result, tau, Q);
    return new @Dimensionless QRDecomposition<@Dimensionless FloatMatrix>(Q, R);
  }
}
