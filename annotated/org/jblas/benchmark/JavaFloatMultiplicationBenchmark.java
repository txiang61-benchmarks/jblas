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

package org.jblas.benchmark;

import units.qual.Dimensionless;
import static org.jblas.FloatMatrix.randn;

/**
 *
 */
@Dimensionless
class JavaFloatMultiplicationBenchmark implements Benchmark {

    public @Dimensionless String getName(@Dimensionless JavaFloatMultiplicationBenchmark this) {
        return "Java matrix multiplication, single precision";
    }

    /** Compute C = A * B */
    private void mmuli(@Dimensionless JavaFloatMultiplicationBenchmark this, @Dimensionless int n, @Dimensionless float @Dimensionless [] A, @Dimensionless float @Dimensionless [] B, @Dimensionless float @Dimensionless [] C) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n * n; i++) {
            C[i] = ((@Dimensionless int) (0));
        }

        for (@Dimensionless int j = ((@Dimensionless int) (0)); j < n; j++) {
            @Dimensionless
            int jn = j * n;
            for (@Dimensionless int k = ((@Dimensionless int) (0)); k < n; k++) {
                @Dimensionless
                int kn = k * n;
                @Dimensionless
                float bkjn = B[k + jn];
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
                    C[i + jn] += A[i + kn] * bkjn;
                }
            }
        }
    }

    public @Dimensionless BenchmarkResult run(@Dimensionless JavaFloatMultiplicationBenchmark this, @Dimensionless int size, @Dimensionless double seconds) {
        @Dimensionless
        int counter = ((@Dimensionless int) (0));
        @Dimensionless
        long ops = ((@Dimensionless int) (0));

        @Dimensionless
        float @Dimensionless [] A = randn(size, size).data;
        @Dimensionless
        float @Dimensionless [] B = randn(size, size).data;
        @Dimensionless
        float @Dimensionless [] C = randn(size, size).data;

        @Dimensionless
        Timer t = new @Dimensionless Timer();
        t.start();
        while (!t.ranFor(seconds)) {
            mmuli(size, A, B, C);
            counter++;
            ops += ((@Dimensionless long) (2L)) * size * size * size;
        }
        t.stop();

        return new @Dimensionless BenchmarkResult(ops, t.elapsedSeconds(), counter);
    }
}
