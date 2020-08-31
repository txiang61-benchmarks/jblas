// --- BEGIN LICENSE BLOCK ---
/*
 * Copyright (c) 2009-2011, Mikio L. Braun
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
package org.jblas.util;

import units.qual.Dimensionless;
import org.jblas.ComplexDouble;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.NativeBlas;
import org.jblas.DoubleMatrix;

/**
 * Run a few sanity checks on the installation to see whether
 * everything runs as expected.
 *
 * @author Mikio L. Braun
 */
public class SanityChecks {

    public static @Dimensionless int checksFailed;

    public static void check(@Dimensionless String message, @Dimensionless boolean condition) {
        System.out.print(message + "... ");
        if (condition) {
            System.out.println("ok");
        } else {
            System.out.println("failed");
            checksFailed++;
        }
    }

    /** Check whether vector addition works. This is pure Java code and should work. */
    public static void checkVectorAddition() {
        @Dimensionless
        DoubleMatrix x = new @Dimensionless DoubleMatrix(((@Dimensionless int) (3)), ((@Dimensionless int) (1)), ((@Dimensionless double) (1.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (3.0)));
        @Dimensionless
        DoubleMatrix y = new @Dimensionless DoubleMatrix(((@Dimensionless int) (3)), ((@Dimensionless int) (1)), ((@Dimensionless double) (4.0)), ((@Dimensionless double) (5.0)), ((@Dimensionless double) (6.0)));
        @Dimensionless
        DoubleMatrix z = new @Dimensionless DoubleMatrix(((@Dimensionless int) (3)), ((@Dimensionless int) (1)), ((@Dimensionless double) (5.0)), ((@Dimensionless double) (7.0)), ((@Dimensionless double) (9.0)));

        check("checking vector addition", x.add(y).equals(z));
    }

    /** Check matrix multiplication. This is already ATLAS/BLAS code. */
    public static void checkMatrixMultiplication() {

        @Dimensionless
        DoubleMatrix A = new @Dimensionless DoubleMatrix(new @Dimensionless double @Dimensionless [] @Dimensionless []{
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (1.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (3.0))},
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (4.0)), ((@Dimensionless double) (5.0)), ((@Dimensionless double) (6.0))},
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (7.0)), ((@Dimensionless double) (8.0)), ((@Dimensionless double) (9.0))}
                });
        @Dimensionless
        DoubleMatrix E = new @Dimensionless DoubleMatrix(new @Dimensionless double @Dimensionless [] @Dimensionless []{
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (0.0)), ((@Dimensionless double) (0.0)), ((@Dimensionless double) (1.0))},
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (0.0)), ((@Dimensionless double) (1.0)), ((@Dimensionless double) (0.0))},
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (1.0)), ((@Dimensionless double) (0.0)), ((@Dimensionless double) (0.0))}
                });
        @Dimensionless
        DoubleMatrix B = new @Dimensionless DoubleMatrix(new @Dimensionless double @Dimensionless [] @Dimensionless []{
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (3.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (1.0))},
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (6.0)), ((@Dimensionless double) (5.0)), ((@Dimensionless double) (4.0))},
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (9.0)), ((@Dimensionless double) (8.0)), ((@Dimensionless double) (7.0))}
                });

        check("checking matrix multiplication", A.mmul(E).equals(B));
    }

    /**
     * Check whether error handling works. If it works, you should see an
     * ok, otherwise, you might see the actual error message and then
     * the program exits.
     */
    public static void checkXerbla() {
        @Dimensionless
        double @Dimensionless [] x = new @Dimensionless double @Dimensionless [((@Dimensionless int) (9))];
        System.out.println("Check whether we're catching XERBLA errors. If you see something like \"** On entry to DGEMM  parameter number  4 had an illegal value\", it didn't work!");
        try {
            NativeBlas.dgemm('N', 'N', ((@Dimensionless int) (3)), ((@Dimensionless int) (-1)), ((@Dimensionless int) (3)), ((@Dimensionless double) (1.0)), x, ((@Dimensionless int) (0)), ((@Dimensionless int) (3)), x, ((@Dimensionless int) (0)), ((@Dimensionless int) (3)), ((@Dimensionless double) (0.0)), x, ((@Dimensionless int) (0)), ((@Dimensionless int) (3)));
        } catch (@Dimensionless IllegalArgumentException e) {
            check("checking XERBLA", e.getMessage().contains("XERBLA"));
            return;
        }
        assert (false); // shouldn't happen
    }

    /**
     * Compute eigenvalues. This is a routine not in ATLAS, but in the original
     * LAPACK.
     */
    public static void checkEigenvalues() {
        @Dimensionless
        DoubleMatrix A = new @Dimensionless DoubleMatrix(new @Dimensionless double @Dimensionless [] @Dimensionless []{
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (3.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (0.0))},
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (2.0)), ((@Dimensionless double) (3.0)), ((@Dimensionless double) (2.0))},
                    new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (0.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (3.0))}
                });

        @Dimensionless
        DoubleMatrix E = new @Dimensionless DoubleMatrix(((@Dimensionless int) (3)), ((@Dimensionless int) (1)));

        NativeBlas.dsyev('N', 'U', ((@Dimensionless int) (3)), A.data, ((@Dimensionless int) (0)), ((@Dimensionless int) (3)), E.data, ((@Dimensionless int) (0)));
        check("checking existence of dsyev...", true);
    }

    public static void checkSVD() {
        @Dimensionless
        double @Dimensionless [] @Dimensionless [] data = new @Dimensionless double @Dimensionless [] @Dimensionless []{
            new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (1.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (3.0))},
            new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (4.0)), ((@Dimensionless double) (5.0)), ((@Dimensionless double) (6.0))},
            new @Dimensionless double @Dimensionless [] { ((@Dimensionless double) (7.0)), ((@Dimensionless double) (8.0)), ((@Dimensionless double) (9.0))},
            new @Dimensionless double @Dimensionless [] {- ((@Dimensionless double) (1.0)), - ((@Dimensionless double) (2.0)), - ((@Dimensionless double) (3.0))}
        };

        @Dimensionless
        DoubleMatrix A = new @Dimensionless DoubleMatrix(data);

        @Dimensionless
        DoubleMatrix @Dimensionless [] USV = org.jblas.Singular.sparseSVD(A);
        System.out.println(USV[((@Dimensionless int) (0))].toString());
        System.out.println(USV[((@Dimensionless int) (1))].toString());
        System.out.println(USV[((@Dimensionless int) (2))].toString());

        System.out.println(org.jblas.Singular.SVDValues(A));

        /*ComplexDoubleMatrix[] AZB = org.jblas.Singular.sparseSVD(new ComplexDoubleMatrix(data));
        System.out.println(AZB[0].toString());
        System.out.println(AZB[1].toString());
        System.out.println(AZB[2].toString());*/
        check("checking existence of dgesvd...", true);
    }

    public static void checkGeneralizedEigenvalues() {
        @Dimensionless
        DoubleMatrix A = new @Dimensionless DoubleMatrix(((@Dimensionless int) (3)), ((@Dimensionless int) (3)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (1.0)), ((@Dimensionless double) (0.0)), ((@Dimensionless double) (1.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (1.0)), ((@Dimensionless double) (0.0)), ((@Dimensionless double) (1.0)), ((@Dimensionless double) (2.0)));
        @Dimensionless
        DoubleMatrix B = new @Dimensionless DoubleMatrix(((@Dimensionless int) (3)), ((@Dimensionless int) (3)), ((@Dimensionless double) (4.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (1.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (4.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (1.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (4.0)));

        @Dimensionless
        DoubleMatrix @Dimensionless [] LA = org.jblas.Eigen.symmetricGeneralizedEigenvectors(A, B);

        check("checkign existence of gsyevd (generalized eigenvalues)...", true);
    }

    public static void checkComplexReturnValues() {
        @Dimensionless
        double @Dimensionless [] data = new @Dimensionless double @Dimensionless [] {
            ((@Dimensionless double) (1.0)), ((@Dimensionless double) (2.0)), ((@Dimensionless double) (3.0)), ((@Dimensionless double) (4.0)), ((@Dimensionless double) (5.0)), ((@Dimensionless double) (6.0))
        };

        @Dimensionless
        ComplexDoubleMatrix A = new @Dimensionless ComplexDoubleMatrix(data);

        @Dimensionless
        ComplexDouble z = A.dotu(A);

        System.out.print("Checking complex return values... ");
        System.out.println("(z = " + z.toString() + ")");
    }

    public static void main(String[] args) {
        Logger.getLogger().setLevel(((@Dimensionless CONFIG) (Logger.CONFIG)));
        for (@Dimensionless String arg : args) {
            if (arg.equals("--debug")) {
                Logger.getLogger().setLevel(((@Dimensionless DEBUG) (Logger.DEBUG)));
            }
        }
        checkVectorAddition();
        checkMatrixMultiplication();
        checkEigenvalues();
        checkSVD();
        checkComplexReturnValues();
        checkXerbla();
        printSummary();
    }

    private static void printSummary() {
        if (checksFailed == ((@Dimensionless int) (0))) {
            System.out.println("Sanity checks passed.");
        } else {
            System.out.println("Sainty checks FAILED!");
        }
    }
}
