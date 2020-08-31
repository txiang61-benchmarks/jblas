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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.jblas;

import units.qual.Dimensionless;
import org.jblas.exceptions.LapackException;

/**
 * <p>Implementation of some Blas functions, mostly those which require linear runtime
 * in the number of matrix elements. Because of the copying overhead when passing
 * primitive arrays to native code, it doesn't make sense for these functions
 * to be implemented in native code. The Java code is about as fast.</p>
 * 
 * <p>The same conventions were used as in the native code, that is, for each array
 * you also pass an index pointing to the starting index.</p>
 * 
 * <p>These methods are mostly optimized for the case where the starting index is 0
 * and the increment is 1.</p>
 */
public class JavaBlas {

    /** Exchange two vectors. */
    public static void rswap(@Dimensionless int n, @Dimensionless double @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy) {
        if (incx == ((@Dimensionless int) (1)) && incy == ((@Dimensionless int) (1)) && dxIdx == ((@Dimensionless int) (0)) && dyIdx == ((@Dimensionless int) (0))) {
            @Dimensionless
            double z;
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
                z = dx[i];
                dx[i] = dy[i];
                dy[i] = z;
            }
        } else {
            @Dimensionless
            double z;
            for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; xi += incx, yi += incy, c++) {
                z = dx[xi];
                dx[xi] = dy[yi];
                dy[yi] = z;
            }
        }
    }

    /** Copy dx to dy. */
    public static void rcopy(int n, double[] dx, int dxIdx, int incx, double[] dy, int dyIdx, int incy) {
        if (dxIdx < ((@Dimensionless int) (0)) || dxIdx + (n - ((@Dimensionless int) (1))) * incx >= dx.length) {
            throw new @Dimensionless LapackException("Java.raxpy", "Parameters for x aren't valid! (n = " + n + ", dx.length = " + dx.length + ", dxIdx = " + dxIdx + ", incx = " + incx + ")");
        }
        if (dyIdx < ((@Dimensionless int) (0)) || dyIdx + (n - ((@Dimensionless int) (1))) * incy >= dy.length) {
            throw new @Dimensionless LapackException("Java.raxpy", "Parameters for y aren't valid! (n = " + n + ", dy.length = " + dy.length + ", dyIdx = " + dyIdx + ", incy = " + incy + ")");
        }
        if (incx == ((@Dimensionless int) (1)) && incy == ((@Dimensionless int) (1))) {
            System.arraycopy(dx, dxIdx, dy, dyIdx, n);
        } else {
            for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; xi += incx, yi += incy, c++) {
                dy[yi] = dx[xi];
            }
        }
    }

    /** Compute dy <- da * dx + dy. */
    public static void raxpy(@Dimensionless int n, @Dimensionless double da, @Dimensionless double @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy) {
        if (dxIdx < ((@Dimensionless int) (0)) || dxIdx + (n - ((@Dimensionless int) (1))) * incx >= dx.length) {
            throw new @Dimensionless LapackException("Java.raxpy", "Parameters for x aren't valid! (n = " + n + ", dx.length = " + dx.length + ", dxIdx = " + dxIdx + ", incx = " + incx + ")");
        }
        if (dyIdx < ((@Dimensionless int) (0)) || dyIdx + (n - ((@Dimensionless int) (1))) * incy >= dy.length) {
            throw new @Dimensionless LapackException("Java.raxpy", "Parameters for y aren't valid! (n = " + n + ", dy.length = " + dy.length + ", dyIdx = " + dyIdx + ", incy = " + incy + ")");
        }
        
        if (incx == ((@Dimensionless int) (1)) && incy == ((@Dimensionless int) (1)) && dxIdx == ((@Dimensionless int) (0)) && dyIdx == ((@Dimensionless int) (0))) {
            if (da == ((@Dimensionless double) (1.0))) {
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
                    dy[i] += dx[i];
                }
            } else {
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
                    dy[i] += da * dx[i];
                }
            }
        } else {
            if (da == ((@Dimensionless double) (1.0))) {
                for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; c++, xi += incx, yi += incy) {
                    dy[yi] += dx[xi];
                }

            } else {
                for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; c++, xi += incx, yi += incy) {
                    dy[yi] += da * dx[xi];
                }
            }
        }
    }

    /** Computes dz <- dx + dy */
    public static void rzaxpy(@Dimensionless int n, @Dimensionless double @Dimensionless [] dz, @Dimensionless int dzIdx, @Dimensionless int incz, @Dimensionless double da, @Dimensionless double @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy) {
        if (dxIdx == ((@Dimensionless int) (0)) && incx == ((@Dimensionless int) (1)) && dyIdx == ((@Dimensionless int) (0)) && incy == ((@Dimensionless int) (1)) && dzIdx == ((@Dimensionless int) (0)) && incz == ((@Dimensionless int) (1))) {
            if (da == ((@Dimensionless double) (1.0))) {
                for (@Dimensionless int c = ((@Dimensionless int) (0)); c < n; c++)
                    dz[c] = dx[c] + dy[c];
            } else {
                for (@Dimensionless int c = ((@Dimensionless int) (0)); c < n; c++)
                    dz[c] = da*dx[c] + dy[c];
            }
        } else {
            if (da == ((@Dimensionless double) (1.0))) {
                for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx, zi = dzIdx; c < n; c++, xi += incx, yi += incy, zi += incz) {
                    dz[zi] = dx[xi] + dy[yi];
                }
            } else {
                for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx, zi = dzIdx; c < n; c++, xi += incx, yi += incy, zi += incz) {
                    dz[zi] = da*dx[xi] + dy[yi];
                }
            }
        }
    }

    public static void rzgxpy(int n, double[] dz, double[] dx, double[] dy) {
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < n; c++)
            dz[c] = dx[c] + dy[c];       
    }

    /** Compute scalar product between dx and dy. */
    public static @Dimensionless double rdot(@Dimensionless int n, @Dimensionless double @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless double @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy) {
        @Dimensionless
        double s = ((@Dimensionless double) (0.0));
        if (incx == ((@Dimensionless int) (1)) && incy == ((@Dimensionless int) (1)) && dxIdx == ((@Dimensionless int) (0)) && dyIdx == ((@Dimensionless int) (0))) {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++)
                s += dx[i] * dy[i];
        }
        else {
            for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; c++, xi += incx, yi += incy) {
                s += dx[xi] * dy[yi];
            }
        }
        return s;
    }
//BEGIN
  // The code below has been automatically generated.
  // DO NOT EDIT!

    /** Exchange two vectors. */
    public static void rswap(@Dimensionless int n, @Dimensionless float @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy) {
        if (incx == ((@Dimensionless int) (1)) && incy == ((@Dimensionless int) (1)) && dxIdx == ((@Dimensionless int) (0)) && dyIdx == ((@Dimensionless int) (0))) {
            @Dimensionless
            float z;
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
                z = dx[i];
                dx[i] = dy[i];
                dy[i] = z;
            }
        } else {
            @Dimensionless
            float z;
            for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; xi += incx, yi += incy, c++) {
                z = dx[xi];
                dx[xi] = dy[yi];
                dy[yi] = z;
            }
        }
    }

    /** Copy dx to dy. */
    public static void rcopy(int n, float[] dx, int dxIdx, int incx, float[] dy, int dyIdx, int incy) {
        if (dxIdx < ((@Dimensionless int) (0)) || dxIdx + (n - ((@Dimensionless int) (1))) * incx >= dx.length) {
            throw new @Dimensionless LapackException("Java.raxpy", "Parameters for x aren't valid! (n = " + n + ", dx.length = " + dx.length + ", dxIdx = " + dxIdx + ", incx = " + incx + ")");
        }
        if (dyIdx < ((@Dimensionless int) (0)) || dyIdx + (n - ((@Dimensionless int) (1))) * incy >= dy.length) {
            throw new @Dimensionless LapackException("Java.raxpy", "Parameters for y aren't valid! (n = " + n + ", dy.length = " + dy.length + ", dyIdx = " + dyIdx + ", incy = " + incy + ")");
        }
        if (incx == ((@Dimensionless int) (1)) && incy == ((@Dimensionless int) (1))) {
            System.arraycopy(dx, dxIdx, dy, dyIdx, n);
        } else {
            for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; xi += incx, yi += incy, c++) {
                dy[yi] = dx[xi];
            }
        }
    }

    /** Compute dy <- da * dx + dy. */
    public static void raxpy(@Dimensionless int n, @Dimensionless float da, @Dimensionless float @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy) {
        if (dxIdx < ((@Dimensionless int) (0)) || dxIdx + (n - ((@Dimensionless int) (1))) * incx >= dx.length) {
            throw new @Dimensionless LapackException("Java.raxpy", "Parameters for x aren't valid! (n = " + n + ", dx.length = " + dx.length + ", dxIdx = " + dxIdx + ", incx = " + incx + ")");
        }
        if (dyIdx < ((@Dimensionless int) (0)) || dyIdx + (n - ((@Dimensionless int) (1))) * incy >= dy.length) {
            throw new @Dimensionless LapackException("Java.raxpy", "Parameters for y aren't valid! (n = " + n + ", dy.length = " + dy.length + ", dyIdx = " + dyIdx + ", incy = " + incy + ")");
        }
        
        if (incx == ((@Dimensionless int) (1)) && incy == ((@Dimensionless int) (1)) && dxIdx == ((@Dimensionless int) (0)) && dyIdx == ((@Dimensionless int) (0))) {
            if (da == ((@Dimensionless float) (1.0f))) {
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
                    dy[i] += dx[i];
                }
            } else {
                for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++) {
                    dy[i] += da * dx[i];
                }
            }
        } else {
            if (da == ((@Dimensionless float) (1.0f))) {
                for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; c++, xi += incx, yi += incy) {
                    dy[yi] += dx[xi];
                }

            } else {
                for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; c++, xi += incx, yi += incy) {
                    dy[yi] += da * dx[xi];
                }
            }
        }
    }

    /** Computes dz <- dx + dy */
    public static void rzaxpy(@Dimensionless int n, @Dimensionless float @Dimensionless [] dz, @Dimensionless int dzIdx, @Dimensionless int incz, @Dimensionless float da, @Dimensionless float @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy) {
        if (dxIdx == ((@Dimensionless int) (0)) && incx == ((@Dimensionless int) (1)) && dyIdx == ((@Dimensionless int) (0)) && incy == ((@Dimensionless int) (1)) && dzIdx == ((@Dimensionless int) (0)) && incz == ((@Dimensionless int) (1))) {
            if (da == ((@Dimensionless float) (1.0f))) {
                for (@Dimensionless int c = ((@Dimensionless int) (0)); c < n; c++)
                    dz[c] = dx[c] + dy[c];
            } else {
                for (@Dimensionless int c = ((@Dimensionless int) (0)); c < n; c++)
                    dz[c] = da*dx[c] + dy[c];
            }
        } else {
            if (da == ((@Dimensionless float) (1.0f))) {
                for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx, zi = dzIdx; c < n; c++, xi += incx, yi += incy, zi += incz) {
                    dz[zi] = dx[xi] + dy[yi];
                }
            } else {
                for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx, zi = dzIdx; c < n; c++, xi += incx, yi += incy, zi += incz) {
                    dz[zi] = da*dx[xi] + dy[yi];
                }
            }
        }
    }

    public static void rzgxpy(int n, float[] dz, float[] dx, float[] dy) {
        for (@Dimensionless int c = ((@Dimensionless int) (0)); c < n; c++)
            dz[c] = dx[c] + dy[c];       
    }

    /** Compute scalar product between dx and dy. */
    public static @Dimensionless float rdot(@Dimensionless int n, @Dimensionless float @Dimensionless [] dx, @Dimensionless int dxIdx, @Dimensionless int incx, @Dimensionless float @Dimensionless [] dy, @Dimensionless int dyIdx, @Dimensionless int incy) {
        @Dimensionless
        float s = ((@Dimensionless float) (0.0f));
        if (incx == ((@Dimensionless int) (1)) && incy == ((@Dimensionless int) (1)) && dxIdx == ((@Dimensionless int) (0)) && dyIdx == ((@Dimensionless int) (0))) {
            for (@Dimensionless int i = ((@Dimensionless int) (0)); i < n; i++)
                s += dx[i] * dy[i];
        }
        else {
            for (int c = ((@Dimensionless int) (0)), xi = dxIdx, yi = dyIdx; c < n; c++, xi += incx, yi += incy) {
                s += dx[xi] * dy[yi];
            }
        }
        return s;
    }
//END
}
