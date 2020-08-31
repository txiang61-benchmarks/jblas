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
import org.jblas.util.Logger;

/**
 * Native BLAS and LAPACK functions.
 *
 * <p>The NativeBlas class contains the native BLAS and LAPACK functions. Each
 * Fortran function is mapped to a static method of this class. For each array argument,
 * an additional parameter is introduced which gives the offset from the beginning of
 * the passed array. In C, you would be able to pass a different pointer, but
 * in Java, you can only pass the whole array.</p>
 *
 * <p>Note that due to the way the JNI is usually implemented, the arrays are first
 * copied outside of the JVM before the function is called. This means that
 * functions whose runtime is linear in the amount of memory do usually not
 * run faster just because you are using a native implementation. This holds true
 * for most Level 1 BLAS routines (like vector addition), and unfortunately also
 * for most Level 2 BLAS routines (like matrix-vector multiplication). For these,
 * there exists a class JavaBlas which contains Java implementations.</p>
 *
 * <p>In LAPACK, there exist routines which require workspace to be allocated together
 * with a standard procedure for computing the size of these workspaces. jblas
 * automatically also generates wrappers for these routines with automatic
 * workspace allocation. These routines have the same name, but the workspace
 * arguments are removed.</p>
 *
 * <p>Finally, an example: The fortran routine<pre>
 * SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
 *     DOUBLE PRECISION DA
 *     INTEGER INCX,INCY,N
 *     DOUBLE PRECISION DX(*),DY(*)
 * </pre>
 * becomes <pre>
 * public static native void daxpy(int n, double da, double[] dx, int dxIdx,
 *    int incx, double[] dy, int dyIdx, int incy);
 * </pre>
 */
public class NativeBlas {

  static {
    NativeBlasLibraryLoader.loadLibraryAndCheckErrors();
  }

  private static @Dimensionless int @Dimensionless [] intDummy = new @Dimensionless int @Dimensionless [((@Dimensionless int) (1))];
  private static @Dimensionless double @Dimensionless [] doubleDummy = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
  private static @Dimensionless float @Dimensionless [] floatDummy = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];

  public static native void ccopy(int n, float[] cx, int cxIdx, int incx, float[] cy, int cyIdx, int incy);
  public static native void dcopy(int n, double[] dx, int dxIdx, int incx, double[] dy, int dyIdx, int incy);
  public static native void scopy(int n, float[] sx, int sxIdx, int incx, float[] sy, int syIdx, int incy);
  public static native void zcopy(int n, double[] zx, int zxIdx, int incx, double[] zy, int zyIdx, int incy);
  public static native void cswap(int n, float[] cx, int cxIdx, int incx, float[] cy, int cyIdx, int incy);
  public static native void dswap(int n, double[] dx, int dxIdx, int incx, double[] dy, int dyIdx, int incy);
  public static native void sswap(int n, float[] sx, int sxIdx, int incx, float[] sy, int syIdx, int incy);
  public static native void zswap(int n, double[] zx, int zxIdx, int incx, double[] zy, int zyIdx, int incy);
  public static native void caxpy(int n, ComplexFloat ca, float[] cx, int cxIdx, int incx, float[] cy, int cyIdx, int incy);
  public static native void daxpy(@Dimensionless int n, @Dimensionless double da, @Dimensionless double @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy);
  public static native void saxpy(@Dimensionless int n, @Dimensionless float sa, @Dimensionless float @Dimensionless [] sx, @Dimensionless int sxIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] sy, @Dimensionless int syIdx, @Dimensionless int incy);
  public static native void zaxpy(int n, ComplexDouble za, double[] zx, int zxIdx, int incx, double[] zy, int zyIdx, int incy);
  public static native void cscal(@Dimensionless int n, @Dimensionless ComplexFloat ca, @Dimensionless float @Dimensionless [] cx, @Dimensionless int cxIdx, @Dimensionless int incx);
  public static native void dscal(int n, double da, double[] dx, int dxIdx, int incx);
  public static native void sscal(int n, float sa, float[] sx, int sxIdx, int incx);
  public static native void zscal(@Dimensionless int n, @Dimensionless ComplexDouble za, @Dimensionless double @Dimensionless [] zx, @Dimensionless int zxIdx, @Dimensionless int incx);
  public static native void csscal(@Dimensionless int n, @Dimensionless float sa, @Dimensionless float @Dimensionless [] cx, @Dimensionless int cxIdx, @Dimensionless int incx);
  public static native void zdscal(@Dimensionless int n, @Dimensionless double da, @Dimensionless double @Dimensionless [] zx, @Dimensionless int zxIdx, @Dimensionless int incx);
  public static native @Dimensionless ComplexFloat cdotc(@Dimensionless int n, @Dimensionless float @Dimensionless [] cx, @Dimensionless int cxIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] cy, @Dimensionless int cyIdx, @Dimensionless int incy);
  public static native @Dimensionless ComplexFloat cdotu(@Dimensionless int n, @Dimensionless float @Dimensionless [] cx, @Dimensionless int cxIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] cy, @Dimensionless int cyIdx, @Dimensionless int incy);
  public static native @Dimensionless double ddot(@Dimensionless int n, @Dimensionless double @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy);
  public static native @Dimensionless float sdot(@Dimensionless int n, @Dimensionless float @Dimensionless [] sx, @Dimensionless int sxIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] sy, @Dimensionless int syIdx, @Dimensionless int incy);
  public static native @Dimensionless ComplexDouble zdotc(@Dimensionless int n, @Dimensionless double @Dimensionless [] zx, @Dimensionless int zxIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] zy, @Dimensionless int zyIdx, @Dimensionless int incy);
  public static native @Dimensionless ComplexDouble zdotu(@Dimensionless int n, @Dimensionless double @Dimensionless [] zx, @Dimensionless int zxIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] zy, @Dimensionless int zyIdx, @Dimensionless int incy);
  public static native @Dimensionless double dnrm2(@Dimensionless int n, @Dimensionless double @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx);
  public static native @Dimensionless double dznrm2(@Dimensionless int n, @Dimensionless double @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx);
  public static native @Dimensionless float scnrm2(@Dimensionless int n, @Dimensionless float @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx);
  public static native @Dimensionless float snrm2(@Dimensionless int n, @Dimensionless float @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx);
  public static native @Dimensionless double dasum(@Dimensionless int n, @Dimensionless double @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx);
  public static native @Dimensionless double dzasum(@Dimensionless int n, @Dimensionless double @Dimensionless [] zx, @Dimensionless int zxIdx, @Dimensionless int incx);
  public static native @Dimensionless float sasum(@Dimensionless int n, @Dimensionless float @Dimensionless [] sx, @Dimensionless int sxIdx, @Dimensionless int incx);
  public static native @Dimensionless float scasum(@Dimensionless int n, @Dimensionless float @Dimensionless [] cx, @Dimensionless int cxIdx, @Dimensionless int incx);
  public static native @Dimensionless int icamax(@Dimensionless int n, @Dimensionless float @Dimensionless [] cx, @Dimensionless int cxIdx, @Dimensionless int incx);
  public static native @Dimensionless int idamax(@Dimensionless int n, @Dimensionless double @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx);
  public static native @Dimensionless int isamax(@Dimensionless int n, @Dimensionless float @Dimensionless [] sx, @Dimensionless int sxIdx, @Dimensionless int incx);
  public static native @Dimensionless int izamax(@Dimensionless int n, @Dimensionless double @Dimensionless [] zx, @Dimensionless int zxIdx, @Dimensionless int incx);
  public static native void cgemv(@Dimensionless char trans, @Dimensionless int m, @Dimensionless int n, @Dimensionless ComplexFloat alpha, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless ComplexFloat beta, @Dimensionless float @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy);
  public static native void dgemv(@Dimensionless char trans, @Dimensionless int m, @Dimensionless int n, @Dimensionless double alpha, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless double beta, @Dimensionless double @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy);
  public static native void sgemv(@Dimensionless char trans, @Dimensionless int m, @Dimensionless int n, @Dimensionless float alpha, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless float beta, @Dimensionless float @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy);
  public static native void zgemv(@Dimensionless char trans, @Dimensionless int m, @Dimensionless int n, @Dimensionless ComplexDouble alpha, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless ComplexDouble beta, @Dimensionless double @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy);
  public static native void cgerc(@Dimensionless int m, @Dimensionless int n, @Dimensionless ComplexFloat alpha, @Dimensionless float @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda);
  public static native void cgeru(@Dimensionless int m, @Dimensionless int n, @Dimensionless ComplexFloat alpha, @Dimensionless float @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda);
  public static native void dger(@Dimensionless int m, @Dimensionless int n, @Dimensionless double alpha, @Dimensionless double @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda);
  public static native void sger(@Dimensionless int m, @Dimensionless int n, @Dimensionless float alpha, @Dimensionless float @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda);
  public static native void zgerc(@Dimensionless int m, @Dimensionless int n, @Dimensionless ComplexDouble alpha, @Dimensionless double @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda);
  public static native void zgeru(@Dimensionless int m, @Dimensionless int n, @Dimensionless ComplexDouble alpha, @Dimensionless double @Dimensionless [] x, @Dimensionless int xIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] y, @Dimensionless int yIdx, @Dimensionless int incy, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda);
  public static native void cgemm(@Dimensionless char transa, @Dimensionless char transb, @Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless ComplexFloat alpha, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless ComplexFloat beta, @Dimensionless float @Dimensionless [] c, @Dimensionless int cIdx, @Dimensionless int ldc);
  public static native void dgemm(char transa, char transb, int m, int n, int k, double alpha, double[] a, int aIdx, int lda, double[] b, int bIdx, int ldb, double beta, double[] c, int cIdx, int ldc);
  public static native void sgemm(@Dimensionless char transa, @Dimensionless char transb, @Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless float alpha, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless float beta, @Dimensionless float @Dimensionless [] c, @Dimensionless int cIdx, @Dimensionless int ldc);
  public static native void zgemm(@Dimensionless char transa, @Dimensionless char transb, @Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless ComplexDouble alpha, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless ComplexDouble beta, @Dimensionless double @Dimensionless [] c, @Dimensionless int cIdx, @Dimensionless int ldc);
  public static native @Dimensionless int dgesv(@Dimensionless int n, @Dimensionless int nrhs, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless int @Dimensionless [] ipiv, @Dimensionless int ipivIdx, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb);
  public static native @Dimensionless int sgesv(@Dimensionless int n, @Dimensionless int nrhs, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless int @Dimensionless [] ipiv, @Dimensionless int ipivIdx, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb);
  public static native @Dimensionless int dsysv(@Dimensionless char uplo, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless int @Dimensionless [] ipiv, @Dimensionless int ipivIdx, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int dsysv(@Dimensionless char uplo, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless int @Dimensionless [] ipiv, @Dimensionless int ipivIdx, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dsysv(uplo, n, nrhs, doubleDummy, ((@Dimensionless int) (0)), lda, intDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), ldb, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dsysv(uplo, n, nrhs, a, aIdx, lda, ipiv, ipivIdx, b, bIdx, ldb, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int ssysv(@Dimensionless char uplo, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless int @Dimensionless [] ipiv, @Dimensionless int ipivIdx, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int ssysv(@Dimensionless char uplo, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless int @Dimensionless [] ipiv, @Dimensionless int ipivIdx, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = ssysv(uplo, n, nrhs, floatDummy, ((@Dimensionless int) (0)), lda, intDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), ldb, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = ssysv(uplo, n, nrhs, a, aIdx, lda, ipiv, ipivIdx, b, bIdx, ldb, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int dsyev(@Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int dsyev(@Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dsyev(jobz, uplo, n, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dsyev(jobz, uplo, n, a, aIdx, lda, w, wIdx, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int ssyev(@Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int ssyev(@Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = ssyev(jobz, uplo, n, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = ssyev(jobz, uplo, n, a, aIdx, lda, w, wIdx, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int dsyevd(@Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int liwork);
  public static @Dimensionless int dsyevd(@Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    @Dimensionless
    int @Dimensionless [] iwork = new @Dimensionless int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int liwork;
    info = dsyevd(jobz, uplo, n, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), iwork, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    liwork = (@Dimensionless int) iwork[((@Dimensionless int) (0))]; iwork = new @Dimensionless int @Dimensionless [liwork];
    info = dsyevd(jobz, uplo, n, a, aIdx, lda, w, wIdx, work, ((@Dimensionless int) (0)), lwork, iwork, ((@Dimensionless int) (0)), liwork);
    return info;
  }

  public static native @Dimensionless int dsyevr(@Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double vl, @Dimensionless double vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless double abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless int @Dimensionless [] isuppz, @Dimensionless int isuppzIdx, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int liwork);
  public static @Dimensionless int dsyevr(@Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double vl, @Dimensionless double vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless double abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless int @Dimensionless [] isuppz, @Dimensionless int isuppzIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    @Dimensionless
    int @Dimensionless [] iwork = new @Dimensionless int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int liwork;
    info = dsyevr(jobz, range, uplo, n, doubleDummy, ((@Dimensionless int) (0)), lda, vl, vu, il, iu, abstol, intDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), ldz, intDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), iwork, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    liwork = (@Dimensionless int) iwork[((@Dimensionless int) (0))]; iwork = new @Dimensionless int @Dimensionless [liwork];
    info = dsyevr(jobz, range, uplo, n, a, aIdx, lda, vl, vu, il, iu, abstol, m, mIdx, w, wIdx, z, zIdx, ldz, isuppz, isuppzIdx, work, ((@Dimensionless int) (0)), lwork, iwork, ((@Dimensionless int) (0)), liwork);
    return info;
  }

  public static native @Dimensionless int dsyevx(@Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double vl, @Dimensionless double vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless double abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int @Dimensionless [] ifail, @Dimensionless int ifailIdx);
  public static @Dimensionless int dsyevx(@Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double vl, @Dimensionless double vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless double abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int @Dimensionless [] ifail, @Dimensionless int ifailIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dsyevx(jobz, range, uplo, n, doubleDummy, ((@Dimensionless int) (0)), lda, vl, vu, il, iu, abstol, intDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), ldz, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), intDummy, ((@Dimensionless int) (0)), intDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dsyevx(jobz, range, uplo, n, a, aIdx, lda, vl, vu, il, iu, abstol, m, mIdx, w, wIdx, z, zIdx, ldz, work, ((@Dimensionless int) (0)), lwork, iwork, iworkIdx, ifail, ifailIdx);
    return info;
  }

  public static native @Dimensionless int ssyevd(@Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int liwork);
  public static @Dimensionless int ssyevd(@Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    @Dimensionless
    int @Dimensionless [] iwork = new @Dimensionless int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int liwork;
    info = ssyevd(jobz, uplo, n, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), iwork, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    liwork = (@Dimensionless int) iwork[((@Dimensionless int) (0))]; iwork = new @Dimensionless int @Dimensionless [liwork];
    info = ssyevd(jobz, uplo, n, a, aIdx, lda, w, wIdx, work, ((@Dimensionless int) (0)), lwork, iwork, ((@Dimensionless int) (0)), liwork);
    return info;
  }

  public static native @Dimensionless int ssyevr(@Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float vl, @Dimensionless float vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless float abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless int @Dimensionless [] isuppz, @Dimensionless int isuppzIdx, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int liwork);
  public static @Dimensionless int ssyevr(@Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float vl, @Dimensionless float vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless float abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless int @Dimensionless [] isuppz, @Dimensionless int isuppzIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    @Dimensionless
    int @Dimensionless [] iwork = new @Dimensionless int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int liwork;
    info = ssyevr(jobz, range, uplo, n, floatDummy, ((@Dimensionless int) (0)), lda, vl, vu, il, iu, abstol, intDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), ldz, intDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), iwork, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    liwork = (@Dimensionless int) iwork[((@Dimensionless int) (0))]; iwork = new @Dimensionless int @Dimensionless [liwork];
    info = ssyevr(jobz, range, uplo, n, a, aIdx, lda, vl, vu, il, iu, abstol, m, mIdx, w, wIdx, z, zIdx, ldz, isuppz, isuppzIdx, work, ((@Dimensionless int) (0)), lwork, iwork, ((@Dimensionless int) (0)), liwork);
    return info;
  }

  public static native @Dimensionless int ssyevx(@Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float vl, @Dimensionless float vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless float abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int @Dimensionless [] ifail, @Dimensionless int ifailIdx);
  public static @Dimensionless int ssyevx(@Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float vl, @Dimensionless float vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless float abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int @Dimensionless [] ifail, @Dimensionless int ifailIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = ssyevx(jobz, range, uplo, n, floatDummy, ((@Dimensionless int) (0)), lda, vl, vu, il, iu, abstol, intDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), ldz, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), intDummy, ((@Dimensionless int) (0)), intDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = ssyevx(jobz, range, uplo, n, a, aIdx, lda, vl, vu, il, iu, abstol, m, mIdx, w, wIdx, z, zIdx, ldz, work, ((@Dimensionless int) (0)), lwork, iwork, iworkIdx, ifail, ifailIdx);
    return info;
  }

  public static native @Dimensionless int dposv(@Dimensionless char uplo, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb);
  public static native @Dimensionless int sposv(@Dimensionless char uplo, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb);
  public static native @Dimensionless int cgeev(@Dimensionless char jobvl, @Dimensionless char jobvr, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] vl, @Dimensionless int vlIdx, @Dimensionless int ldvl, @Dimensionless float @Dimensionless [] vr, @Dimensionless int vrIdx, @Dimensionless int ldvr, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless float @Dimensionless [] rwork, @Dimensionless int rworkIdx);
  public static @Dimensionless int cgeev(@Dimensionless char jobvl, @Dimensionless char jobvr, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] vl, @Dimensionless int vlIdx, @Dimensionless int ldvl, @Dimensionless float @Dimensionless [] vr, @Dimensionless int vrIdx, @Dimensionless int ldvr, @Dimensionless float @Dimensionless [] rwork, @Dimensionless int rworkIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))* ((@Dimensionless int) (2))];
    @Dimensionless
    int lwork;
    info = cgeev(jobvl, jobvr, n, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), ldvl, floatDummy, ((@Dimensionless int) (0)), ldvr, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), floatDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork* ((@Dimensionless int) (2))];
    info = cgeev(jobvl, jobvr, n, a, aIdx, lda, w, wIdx, vl, vlIdx, ldvl, vr, vrIdx, ldvr, work, ((@Dimensionless int) (0)), lwork, rwork, rworkIdx);
    return info;
  }

  public static native @Dimensionless int dgeev(@Dimensionless char jobvl, @Dimensionless char jobvr, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] wr, @Dimensionless int wrIdx, @Dimensionless double @Dimensionless [] wi, @Dimensionless int wiIdx, @Dimensionless double @Dimensionless [] vl, @Dimensionless int vlIdx, @Dimensionless int ldvl, @Dimensionless double @Dimensionless [] vr, @Dimensionless int vrIdx, @Dimensionless int ldvr, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int dgeev(@Dimensionless char jobvl, @Dimensionless char jobvr, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] wr, @Dimensionless int wrIdx, @Dimensionless double @Dimensionless [] wi, @Dimensionless int wiIdx, @Dimensionless double @Dimensionless [] vl, @Dimensionless int vlIdx, @Dimensionless int ldvl, @Dimensionless double @Dimensionless [] vr, @Dimensionless int vrIdx, @Dimensionless int ldvr) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dgeev(jobvl, jobvr, n, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), ldvl, doubleDummy, ((@Dimensionless int) (0)), ldvr, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dgeev(jobvl, jobvr, n, a, aIdx, lda, wr, wrIdx, wi, wiIdx, vl, vlIdx, ldvl, vr, vrIdx, ldvr, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int sgeev(@Dimensionless char jobvl, @Dimensionless char jobvr, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] wr, @Dimensionless int wrIdx, @Dimensionless float @Dimensionless [] wi, @Dimensionless int wiIdx, @Dimensionless float @Dimensionless [] vl, @Dimensionless int vlIdx, @Dimensionless int ldvl, @Dimensionless float @Dimensionless [] vr, @Dimensionless int vrIdx, @Dimensionless int ldvr, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int sgeev(@Dimensionless char jobvl, @Dimensionless char jobvr, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] wr, @Dimensionless int wrIdx, @Dimensionless float @Dimensionless [] wi, @Dimensionless int wiIdx, @Dimensionless float @Dimensionless [] vl, @Dimensionless int vlIdx, @Dimensionless int ldvl, @Dimensionless float @Dimensionless [] vr, @Dimensionless int vrIdx, @Dimensionless int ldvr) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = sgeev(jobvl, jobvr, n, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), ldvl, floatDummy, ((@Dimensionless int) (0)), ldvr, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = sgeev(jobvl, jobvr, n, a, aIdx, lda, wr, wrIdx, wi, wiIdx, vl, vlIdx, ldvl, vr, vrIdx, ldvr, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int zgeev(@Dimensionless char jobvl, @Dimensionless char jobvr, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] vl, @Dimensionless int vlIdx, @Dimensionless int ldvl, @Dimensionless double @Dimensionless [] vr, @Dimensionless int vrIdx, @Dimensionless int ldvr, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless double @Dimensionless [] rwork, @Dimensionless int rworkIdx);
  public static @Dimensionless int zgeev(@Dimensionless char jobvl, @Dimensionless char jobvr, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] vl, @Dimensionless int vlIdx, @Dimensionless int ldvl, @Dimensionless double @Dimensionless [] vr, @Dimensionless int vrIdx, @Dimensionless int ldvr, @Dimensionless double @Dimensionless [] rwork, @Dimensionless int rworkIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))* ((@Dimensionless int) (2))];
    @Dimensionless
    int lwork;
    info = zgeev(jobvl, jobvr, n, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), ldvl, doubleDummy, ((@Dimensionless int) (0)), ldvr, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), doubleDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork* ((@Dimensionless int) (2))];
    info = zgeev(jobvl, jobvr, n, a, aIdx, lda, w, wIdx, vl, vlIdx, ldvl, vr, vrIdx, ldvr, work, ((@Dimensionless int) (0)), lwork, rwork, rworkIdx);
    return info;
  }

  public static native @Dimensionless int dgetrf(@Dimensionless int m, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless int @Dimensionless [] ipiv, @Dimensionless int ipivIdx);
  public static native @Dimensionless int sgetrf(@Dimensionless int m, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless int @Dimensionless [] ipiv, @Dimensionless int ipivIdx);
  public static native @Dimensionless int dpotrf(@Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda);
  public static native @Dimensionless int spotrf(@Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda);
  public static native @Dimensionless int cgesvd(@Dimensionless char jobu, @Dimensionless char jobvt, @Dimensionless int m, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless float @Dimensionless [] u, @Dimensionless int uIdx, @Dimensionless int ldu, @Dimensionless float @Dimensionless [] vt, @Dimensionless int vtIdx, @Dimensionless int ldvt, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless float @Dimensionless [] rwork, @Dimensionless int rworkIdx);
  public static @Dimensionless int cgesvd(@Dimensionless char jobu, @Dimensionless char jobvt, @Dimensionless int m, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless float @Dimensionless [] u, @Dimensionless int uIdx, @Dimensionless int ldu, @Dimensionless float @Dimensionless [] vt, @Dimensionless int vtIdx, @Dimensionless int ldvt, @Dimensionless float @Dimensionless [] rwork, @Dimensionless int rworkIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))* ((@Dimensionless int) (2))];
    @Dimensionless
    int lwork;
    info = cgesvd(jobu, jobvt, m, n, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), ldu, floatDummy, ((@Dimensionless int) (0)), ldvt, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), floatDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork* ((@Dimensionless int) (2))];
    info = cgesvd(jobu, jobvt, m, n, a, aIdx, lda, s, sIdx, u, uIdx, ldu, vt, vtIdx, ldvt, work, ((@Dimensionless int) (0)), lwork, rwork, rworkIdx);
    return info;
  }

  public static native @Dimensionless int dgesvd(@Dimensionless char jobu, @Dimensionless char jobvt, @Dimensionless int m, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless double @Dimensionless [] u, @Dimensionless int uIdx, @Dimensionless int ldu, @Dimensionless double @Dimensionless [] vt, @Dimensionless int vtIdx, @Dimensionless int ldvt, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int dgesvd(@Dimensionless char jobu, @Dimensionless char jobvt, @Dimensionless int m, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless double @Dimensionless [] u, @Dimensionless int uIdx, @Dimensionless int ldu, @Dimensionless double @Dimensionless [] vt, @Dimensionless int vtIdx, @Dimensionless int ldvt) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dgesvd(jobu, jobvt, m, n, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), ldu, doubleDummy, ((@Dimensionless int) (0)), ldvt, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dgesvd(jobu, jobvt, m, n, a, aIdx, lda, s, sIdx, u, uIdx, ldu, vt, vtIdx, ldvt, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int sgesvd(@Dimensionless char jobu, @Dimensionless char jobvt, @Dimensionless int m, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless float @Dimensionless [] u, @Dimensionless int uIdx, @Dimensionless int ldu, @Dimensionless float @Dimensionless [] vt, @Dimensionless int vtIdx, @Dimensionless int ldvt, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int sgesvd(@Dimensionless char jobu, @Dimensionless char jobvt, @Dimensionless int m, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless float @Dimensionless [] u, @Dimensionless int uIdx, @Dimensionless int ldu, @Dimensionless float @Dimensionless [] vt, @Dimensionless int vtIdx, @Dimensionless int ldvt) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = sgesvd(jobu, jobvt, m, n, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), ldu, floatDummy, ((@Dimensionless int) (0)), ldvt, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = sgesvd(jobu, jobvt, m, n, a, aIdx, lda, s, sIdx, u, uIdx, ldu, vt, vtIdx, ldvt, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int zgesvd(@Dimensionless char jobu, @Dimensionless char jobvt, @Dimensionless int m, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless double @Dimensionless [] u, @Dimensionless int uIdx, @Dimensionless int ldu, @Dimensionless double @Dimensionless [] vt, @Dimensionless int vtIdx, @Dimensionless int ldvt, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless double @Dimensionless [] rwork, @Dimensionless int rworkIdx);
  public static @Dimensionless int zgesvd(@Dimensionless char jobu, @Dimensionless char jobvt, @Dimensionless int m, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless double @Dimensionless [] u, @Dimensionless int uIdx, @Dimensionless int ldu, @Dimensionless double @Dimensionless [] vt, @Dimensionless int vtIdx, @Dimensionless int ldvt, @Dimensionless double @Dimensionless [] rwork, @Dimensionless int rworkIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))* ((@Dimensionless int) (2))];
    @Dimensionless
    int lwork;
    info = zgesvd(jobu, jobvt, m, n, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), ldu, doubleDummy, ((@Dimensionless int) (0)), ldvt, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), doubleDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork* ((@Dimensionless int) (2))];
    info = zgesvd(jobu, jobvt, m, n, a, aIdx, lda, s, sIdx, u, uIdx, ldu, vt, vtIdx, ldvt, work, ((@Dimensionless int) (0)), lwork, rwork, rworkIdx);
    return info;
  }

  public static native @Dimensionless int dsygvd(@Dimensionless int itype, @Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int liwork);
  public static @Dimensionless int dsygvd(@Dimensionless int itype, @Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    @Dimensionless
    int @Dimensionless [] iwork = new @Dimensionless int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int liwork;
    info = dsygvd(itype, jobz, uplo, n, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), ldb, doubleDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), iwork, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    liwork = (@Dimensionless int) iwork[((@Dimensionless int) (0))]; iwork = new @Dimensionless int @Dimensionless [liwork];
    info = dsygvd(itype, jobz, uplo, n, a, aIdx, lda, b, bIdx, ldb, w, wIdx, work, ((@Dimensionless int) (0)), lwork, iwork, ((@Dimensionless int) (0)), liwork);
    return info;
  }

  public static native @Dimensionless int ssygvd(@Dimensionless int itype, @Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int liwork);
  public static @Dimensionless int ssygvd(@Dimensionless int itype, @Dimensionless char jobz, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    @Dimensionless
    int @Dimensionless [] iwork = new @Dimensionless int @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int liwork;
    info = ssygvd(itype, jobz, uplo, n, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), ldb, floatDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), iwork, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    liwork = (@Dimensionless int) iwork[((@Dimensionless int) (0))]; iwork = new @Dimensionless int @Dimensionless [liwork];
    info = ssygvd(itype, jobz, uplo, n, a, aIdx, lda, b, bIdx, ldb, w, wIdx, work, ((@Dimensionless int) (0)), lwork, iwork, ((@Dimensionless int) (0)), liwork);
    return info;
  }

  public static native @Dimensionless int dgelsd(@Dimensionless int m, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless double @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless double rcond, @Dimensionless int @Dimensionless [] rank, @Dimensionless int rankIdx, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx);
  public static @Dimensionless int dgelsd(@Dimensionless int m, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless double @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless double rcond, @Dimensionless int @Dimensionless [] rank, @Dimensionless int rankIdx, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dgelsd(m, n, nrhs, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), ldb, doubleDummy, ((@Dimensionless int) (0)), rcond, intDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), intDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dgelsd(m, n, nrhs, a, aIdx, lda, b, bIdx, ldb, s, sIdx, rcond, rank, rankIdx, work, ((@Dimensionless int) (0)), lwork, iwork, iworkIdx);
    return info;
  }

  public static native @Dimensionless int sgelsd(@Dimensionless int m, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless float @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless float rcond, @Dimensionless int @Dimensionless [] rank, @Dimensionless int rankIdx, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx);
  public static @Dimensionless int sgelsd(@Dimensionless int m, @Dimensionless int n, @Dimensionless int nrhs, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless float @Dimensionless [] s, @Dimensionless int sIdx, @Dimensionless float rcond, @Dimensionless int @Dimensionless [] rank, @Dimensionless int rankIdx, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = sgelsd(m, n, nrhs, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), ldb, floatDummy, ((@Dimensionless int) (0)), rcond, intDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), intDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = sgelsd(m, n, nrhs, a, aIdx, lda, b, bIdx, ldb, s, sIdx, rcond, rank, rankIdx, work, ((@Dimensionless int) (0)), lwork, iwork, iworkIdx);
    return info;
  }

  public static native @Dimensionless int ilaenv(@Dimensionless int ispec, @Dimensionless String name, @Dimensionless String opts, @Dimensionless int n1, @Dimensionless int n2, @Dimensionless int n3, @Dimensionless int n4);
  public static native @Dimensionless int dgeqrf(@Dimensionless int m, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] tau, @Dimensionless int tauIdx, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int dgeqrf(@Dimensionless int m, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] tau, @Dimensionless int tauIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dgeqrf(m, n, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dgeqrf(m, n, a, aIdx, lda, tau, tauIdx, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int sgeqrf(@Dimensionless int m, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] tau, @Dimensionless int tauIdx, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int sgeqrf(@Dimensionless int m, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] tau, @Dimensionless int tauIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = sgeqrf(m, n, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = sgeqrf(m, n, a, aIdx, lda, tau, tauIdx, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int dormqr(@Dimensionless char side, @Dimensionless char trans, @Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] tau, @Dimensionless int tauIdx, @Dimensionless double @Dimensionless [] c, @Dimensionless int cIdx, @Dimensionless int ldc, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int dormqr(@Dimensionless char side, @Dimensionless char trans, @Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] tau, @Dimensionless int tauIdx, @Dimensionless double @Dimensionless [] c, @Dimensionless int cIdx, @Dimensionless int ldc) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dormqr(side, trans, m, n, k, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), ldc, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dormqr(side, trans, m, n, k, a, aIdx, lda, tau, tauIdx, c, cIdx, ldc, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int sormqr(@Dimensionless char side, @Dimensionless char trans, @Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] tau, @Dimensionless int tauIdx, @Dimensionless float @Dimensionless [] c, @Dimensionless int cIdx, @Dimensionless int ldc, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int sormqr(@Dimensionless char side, @Dimensionless char trans, @Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] tau, @Dimensionless int tauIdx, @Dimensionless float @Dimensionless [] c, @Dimensionless int cIdx, @Dimensionless int ldc) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = sormqr(side, trans, m, n, k, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), ldc, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = sormqr(side, trans, m, n, k, a, aIdx, lda, tau, tauIdx, c, cIdx, ldc, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int dorgqr(@Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] tau, @Dimensionless int tauIdx, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int dorgqr(@Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] tau, @Dimensionless int tauIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dorgqr(m, n, k, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dorgqr(m, n, k, a, aIdx, lda, tau, tauIdx, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int sorgqr(@Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] tau, @Dimensionless int tauIdx, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork);
  public static @Dimensionless int sorgqr(@Dimensionless int m, @Dimensionless int n, @Dimensionless int k, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] tau, @Dimensionless int tauIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = sorgqr(m, n, k, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = sorgqr(m, n, k, a, aIdx, lda, tau, tauIdx, work, ((@Dimensionless int) (0)), lwork);
    return info;
  }

  public static native @Dimensionless int dsygvx(@Dimensionless int itype, @Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless double vl, @Dimensionless double vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless double abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless double @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int @Dimensionless [] ifail, @Dimensionless int ifailIdx);
  public static @Dimensionless int dsygvx(@Dimensionless int itype, @Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless double @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless double @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless double vl, @Dimensionless double vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless double abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless double @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless double @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int @Dimensionless [] ifail, @Dimensionless int ifailIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    double @Dimensionless [] work = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = dsygvx(itype, jobz, range, uplo, n, doubleDummy, ((@Dimensionless int) (0)), lda, doubleDummy, ((@Dimensionless int) (0)), ldb, vl, vu, il, iu, abstol, intDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), doubleDummy, ((@Dimensionless int) (0)), ldz, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), intDummy, ((@Dimensionless int) (0)), intDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless double @Dimensionless [lwork];
    info = dsygvx(itype, jobz, range, uplo, n, a, aIdx, lda, b, bIdx, ldb, vl, vu, il, iu, abstol, m, mIdx, w, wIdx, z, zIdx, ldz, work, ((@Dimensionless int) (0)), lwork, iwork, iworkIdx, ifail, ifailIdx);
    return info;
  }

  public static native @Dimensionless int ssygvx(@Dimensionless int itype, @Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless float vl, @Dimensionless float vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless float abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless float @Dimensionless [] work, @Dimensionless int workIdx, @Dimensionless int lwork, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int @Dimensionless [] ifail, @Dimensionless int ifailIdx);
  public static @Dimensionless int ssygvx(@Dimensionless int itype, @Dimensionless char jobz, @Dimensionless char range, @Dimensionless char uplo, @Dimensionless int n, @Dimensionless float @Dimensionless [] a, @Dimensionless int aIdx, @Dimensionless int lda, @Dimensionless float @Dimensionless [] b, @Dimensionless int bIdx, @Dimensionless int ldb, @Dimensionless float vl, @Dimensionless float vu, @Dimensionless int il, @Dimensionless int iu, @Dimensionless float abstol, @Dimensionless int @Dimensionless [] m, @Dimensionless int mIdx, @Dimensionless float @Dimensionless [] w, @Dimensionless int wIdx, @Dimensionless float @Dimensionless [] z, @Dimensionless int zIdx, @Dimensionless int ldz, @Dimensionless int @Dimensionless [] iwork, @Dimensionless int iworkIdx, @Dimensionless int @Dimensionless [] ifail, @Dimensionless int ifailIdx) {
    @Dimensionless
    int info;
    @Dimensionless
    float @Dimensionless [] work = new @Dimensionless float @Dimensionless [((@Dimensionless int) (1))];
    @Dimensionless
    int lwork;
    info = ssygvx(itype, jobz, range, uplo, n, floatDummy, ((@Dimensionless int) (0)), lda, floatDummy, ((@Dimensionless int) (0)), ldb, vl, vu, il, iu, abstol, intDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), floatDummy, ((@Dimensionless int) (0)), ldz, work, ((@Dimensionless int) (0)), ((@Dimensionless int) (-1)), intDummy, ((@Dimensionless int) (0)), intDummy, ((@Dimensionless int) (0)));
    if (info != ((@Dimensionless int) (0)))
      return info;
    lwork = (@Dimensionless int) work[((@Dimensionless int) (0))]; work = new @Dimensionless float @Dimensionless [lwork];
    info = ssygvx(itype, jobz, range, uplo, n, a, aIdx, lda, b, bIdx, ldb, vl, vu, il, iu, abstol, m, mIdx, w, wIdx, z, zIdx, ldz, work, ((@Dimensionless int) (0)), lwork, iwork, iworkIdx, ifail, ifailIdx);
    return info;
  }


}
