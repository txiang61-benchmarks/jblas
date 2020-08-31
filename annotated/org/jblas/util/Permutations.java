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

package org.jblas.util;

import units.qual.Dimensionless;
import java.util.Random;
import org.jblas.DoubleMatrix;
import org.jblas.FloatMatrix;

/**
 * Functions which generate random permutations.
 *
 * @author Mikio L. Braun
 */
@Dimensionless
public class Permutations {
    /**
     * Create a random permutation of the numbers 0, ..., size - 1.
     *
     * see Algorithm P, D.E. Knuth: The Art of Computer Programming, Vol. 2, p. 145
     */
    public static @Dimensionless int @Dimensionless [] randomPermutation(@Dimensionless int size) {
        @Dimensionless
        Random r = new @Dimensionless Random();
        @Dimensionless
        int @Dimensionless [] result = new @Dimensionless int @Dimensionless [size];

        for (@Dimensionless int j = ((@Dimensionless int) (0)); j < size; j++) {
            result[j] = j;
        }
        
        for (@Dimensionless int j = size - ((@Dimensionless int) (1)); j > ((@Dimensionless int) (0)); j--) {
            @Dimensionless
            int k = r.nextInt(j);
            @Dimensionless
            int temp = result[j];
            result[j] = result[k];
            result[k] = temp;
        }

        return result;
    }
    
    /**
     * Get a random sample of k out of n elements.
     *
     * See Algorithm S, D. E. Knuth, The Art of Computer Programming, Vol. 2, p.142.
     */
    public static @Dimensionless int @Dimensionless [] randomSubset(@Dimensionless int k, @Dimensionless int n) {
        assert(((@Dimensionless int) (0)) < k && k <= n);
        @Dimensionless
        Random r = new @Dimensionless Random();
        int t = ((@Dimensionless int) (0)), m = ((@Dimensionless int) (0));
        @Dimensionless
        int @Dimensionless [] result = new @Dimensionless int @Dimensionless [k];

        while (m < k) {
            @Dimensionless
            double u = r.nextDouble();
            if ( (n - t) * u < k - m ) {
                result[m] = t;
                m++;
            }
            t++;
        }
        return result;
    }

    /**
     * Create a permutation matrix from a LAPACK-style 'ipiv' vector.
     *
     * @param ipiv row i was interchanged with row ipiv[i]
     */
    public static @Dimensionless DoubleMatrix permutationDoubleMatrixFromPivotIndices(@Dimensionless int size, @Dimensionless int @Dimensionless [] ipiv) {
        @Dimensionless
        int n = ipiv.length;
        //System.out.printf("size = %d n = %d\n", size, n);
        @Dimensionless
        int indices @Dimensionless [] = new int @Dimensionless [size];
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < size; i++)
            indices[i] = i;

        //for (int i = 0; i < n; i++)
        //    System.out.printf("ipiv[%d] = %d\n", i, ipiv[i]);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
            @Dimensionless
            int j = ipiv[i] - ((@Dimensionless int) (1));
            @Dimensionless
            int t = indices[i];
            indices[i] = indices[j];
            indices[j] = t;
        }
        @Dimensionless
        DoubleMatrix result = new @Dimensionless DoubleMatrix(size, size);
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < size; i++)
            result.put(indices[i], i, ((@Dimensionless double) (1.0)));
        return result;
    }

  /**
   * Create a permutation matrix from a LAPACK-style 'ipiv' vector.
   *
   * @param ipiv row i was interchanged with row ipiv[i]
   */
  public static @Dimensionless FloatMatrix permutationFloatMatrixFromPivotIndices(@Dimensionless int size, @Dimensionless int @Dimensionless [] ipiv) {
      @Dimensionless
      int n = ipiv.length;
      //System.out.printf("size = %d n = %d\n", size, n);
      @Dimensionless
      int indices @Dimensionless [] = new int @Dimensionless [size];
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < size; i++)
          indices[i] = i;

      //for (int i = 0; i < n; i++)
      //    System.out.printf("ipiv[%d] = %d\n", i, ipiv[i]);

      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
          @Dimensionless
          int j = ipiv[i] - ((@Dimensionless int) (1));
          @Dimensionless
          int t = indices[i];
          indices[i] = indices[j];
          indices[j] = t;
      }
      @Dimensionless
      FloatMatrix result = new @Dimensionless FloatMatrix(size, size);
      for (@Dimensionless int i = ((@Dimensionless int) (0)); i < size; i++)
          result.put(indices[i], i, ((@Dimensionless float) (1.0f)));
      return result;
  }
}
