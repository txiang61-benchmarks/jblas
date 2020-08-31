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
import units.qual.rad;

/**
 * A complex value with double precision.
 *
 * @author Mikio L. Braun
 */
public class ComplexDouble {

  private double r, i;
  public static final ComplexDouble UNIT = new @Dimensionless ComplexDouble(((@Dimensionless double) (1.0)), ((@Dimensionless double) (0.0)));
  public static final @Dimensionless ComplexDouble I = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)), ((@Dimensionless double) (1.0)));
  public static final ComplexDouble NEG_UNIT = new @Dimensionless ComplexDouble(- ((@Dimensionless double) (1.0)), ((@Dimensionless double) (0.0)));
  public static final @Dimensionless ComplexDouble NEG_I = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)), - ((@Dimensionless double) (1.0)));
  public static final ComplexDouble ZERO = new @Dimensionless ComplexDouble(((@Dimensionless double) (0.0)));

  public ComplexDouble(double real, double imag) {
    r = real;
    i = imag;
  }

  public ComplexDouble(double real) {
    this(real, ((@Dimensionless double) (0.0)));
  }

  public @Dimensionless String toString(@Dimensionless ComplexDouble this) {
    if (i >= ((@Dimensionless int) (0))) {
      return r + " + " + i + "i";
    } else {
      return r + " - " + (-i) + "i";
    }
  }

  public ComplexDouble set(@Dimensionless ComplexDouble this, double real, double imag) {
    r = real;
    i = imag;
    return this;
  }

  public double real(@Dimensionless ComplexDouble this) {
    return r;
  }

  public double imag(@Dimensionless ComplexDouble this) {
    return i;
  }

  public @Dimensionless ComplexDouble dup(@Dimensionless ComplexDouble this) {
    return new @Dimensionless ComplexDouble(r, i);
  }

  public ComplexDouble copy(@Dimensionless ComplexDouble this, ComplexDouble other) {
    r = other.r;
    i = other.i;
    return this;
  }

  /**
   * Add two complex numbers in-place
   *
   * @param c      other complex number
   * @param result complex number where result is stored
   * @return same as result
   */
  public @Dimensionless ComplexDouble addi(@Dimensionless ComplexDouble this, @Dimensionless ComplexDouble c, @Dimensionless ComplexDouble result) {
    if (this == result) {
      r += c.r;
      i += c.i;
    } else {
      result.r = r + c.r;
      result.i = i + c.i;
    }
    return result;
  }

  /**
   * Add two complex numbers in-place storing the result in this.
   *
   * @param c other complex number
   * @return resulting complex number
   */
  public ComplexDouble addi(@Dimensionless ComplexDouble this, ComplexDouble c) {
    return addi(c, this);
  }

  /**
   * Add two complex numbers.
   *
   * @param c other complex number
   * @return new complex number with result
   */
  public ComplexDouble add(@Dimensionless ComplexDouble this, ComplexDouble c) {
    return dup().addi(c);
  }

  /**
   * Add a real number to a complex number in-place.
   *
   * @param a      real number to add
   * @param result complex number to hold result
   * @return same as result
   */
  public @Dimensionless ComplexDouble addi(@Dimensionless ComplexDouble this, @Dimensionless double a, @Dimensionless ComplexDouble result) {
    if (this == result) {
      r += a;
    } else {
      result.r = r + a;
      result.i = i;
    }
    return result;
  }

  /**
   * Add a real number to complex number in-place, storing the result in this.
   *
   * @param c real number to add
   * @return resulting complex number
   */
  public @Dimensionless ComplexDouble addi(@Dimensionless ComplexDouble this, @Dimensionless double c) {
    return addi(c, this);
  }

  /**
   * Add a real number to a complex number.
   *
   * @param c real number to add
   * @return new complex number with result
   */
  public @Dimensionless ComplexDouble add(@Dimensionless ComplexDouble this, @Dimensionless double c) {
    return dup().addi(c);
  }

  /**
   * Subtract two complex numbers, in-place
   *
   * @param c      complex number to subtract
   * @param result resulting complex number
   * @return same as result
   */
  public @Dimensionless ComplexDouble subi(@Dimensionless ComplexDouble this, @Dimensionless ComplexDouble c, @Dimensionless ComplexDouble result) {
    if (this == result) {
      r -= c.r;
      i -= c.i;
    } else {
      result.r = r - c.r;
      result.i = i - c.i;
    }
    return this;
  }

  public @Dimensionless ComplexDouble subi(@Dimensionless ComplexDouble this, @Dimensionless ComplexDouble c) {
    return subi(c, this);
  }

  /**
   * Subtract two complex numbers
   *
   * @param c complex number to subtract
   * @return new complex number with result
   */
  public ComplexDouble sub(@Dimensionless ComplexDouble this, ComplexDouble c) {
    return dup().subi(c);
  }

  public @Dimensionless ComplexDouble subi(@Dimensionless ComplexDouble this, @Dimensionless double a, @Dimensionless ComplexDouble result) {
    if (this == result) {
      r -= a;
    } else {
      result.r = r - a;
      result.i = i;
    }
    return result;
  }

  public @Dimensionless ComplexDouble subi(@Dimensionless ComplexDouble this, @Dimensionless double a) {
    return subi(a, this);
  }

  public @Dimensionless ComplexDouble sub(@Dimensionless ComplexDouble this, @Dimensionless double r) {
    return dup().subi(r);
  }

  /**
   * Multiply two complex numbers, in-place
   *
   * @param c      other complex number
   * @param result complex number where product is stored
   * @return same as result
   */
  public @Dimensionless ComplexDouble muli(@Dimensionless ComplexDouble this, @Dimensionless ComplexDouble c, @Dimensionless ComplexDouble result) {
    @Dimensionless
    double newR = r * c.r - i * c.i;
    @Dimensionless
    double newI = r * c.i + i * c.r;
    result.r = newR;
    result.i = newI;
    return result;
  }

  public ComplexDouble muli(@Dimensionless ComplexDouble this, ComplexDouble c) {
    return muli(c, this);
  }

  /**
   * Multiply two complex numbers
   *
   * @param c other complex number
   * @return new complex number with product of this and c
   */
  public @Dimensionless ComplexDouble mul(@Dimensionless ComplexDouble this, @Dimensionless ComplexDouble c) {
    return dup().muli(c);
  }

  public @Dimensionless ComplexDouble mul(@Dimensionless ComplexDouble this, @Dimensionless double v) {
    return dup().muli(v);
  }

  public @Dimensionless ComplexDouble muli(@Dimensionless ComplexDouble this, @Dimensionless double v, @Dimensionless ComplexDouble result) {
    if (this == result) {
      r *= v;
      i *= v;
    } else {
      result.r = r * v;
      result.i = i * v;
    }
    return this;
  }

  public @Dimensionless ComplexDouble muli(@Dimensionless ComplexDouble this, @Dimensionless double v) {
    return muli(v, this);
  }

  /**
   * Divide two complex numbers
   *
   * @param c complex number to divide this by
   * @return new complex number with quotient of this and c
   */
  public @Dimensionless ComplexDouble div(@Dimensionless ComplexDouble this, @Dimensionless ComplexDouble c) {
    return dup().divi(c);
  }

  /**
   * Divide two complex numbers, in-place
   *
   * @param c      complex number to divide this by
   * @param result complex number to hold result
   * @return same as result
   */
  public @Dimensionless ComplexDouble divi(@Dimensionless ComplexDouble this, @Dimensionless ComplexDouble c, @Dimensionless ComplexDouble result) {
    @Dimensionless
    double d = c.r * c.r + c.i * c.i;
    @Dimensionless
    double newR = (r * c.r + i * c.i) / d;
    @Dimensionless
    double newI = (i * c.r - r * c.i) / d;
    result.r = newR;
    result.i = newI;
    return result;
  }

  public ComplexDouble divi(@Dimensionless ComplexDouble this, ComplexDouble c) {
    return divi(c, this);
  }

  public @Dimensionless ComplexDouble divi(@Dimensionless ComplexDouble this, @Dimensionless double v, @Dimensionless ComplexDouble result) {
    if (this == result) {
      r /= v;
      i /= v;
    } else {
      result.r = r / v;
      result.i = i / v;
    }
    return this;
  }

  public @Dimensionless ComplexDouble divi(@Dimensionless ComplexDouble this, @Dimensionless double v) {
    return divi(v, this);
  }

  public ComplexDouble div(@Dimensionless ComplexDouble this, double v) {
    return dup().divi(v);
  }

  /**
   * Return the absolute value
   *
   * @return the result (length of the vector in 2d plane)
   */
  public double abs(@Dimensionless ComplexDouble this) {
    return (@Dimensionless double) Math.sqrt(r * r + i * i);
  }

  /**
   * Returns the argument of a complex number.
   *
   * @return the result (angle in radians of the vector in 2d plane)
   */
  public @rad double arg(@Dimensionless ComplexDouble this) {
    return (@rad double) Math.atan2(i, r);
  }

  public @Dimensionless ComplexDouble invi(@Dimensionless ComplexDouble this) {
    @Dimensionless
    double d = r * r + i * i;
    r = r / d;
    i = -i / d;
    return this;
  }

  public @Dimensionless ComplexDouble inv(@Dimensionless ComplexDouble this) {
    return dup().invi();
  }

  public @Dimensionless ComplexDouble neg(@Dimensionless ComplexDouble this) {
    return dup().negi();
  }

  public ComplexDouble negi(@Dimensionless ComplexDouble this) {
    r = -r;
    i = -i;
    return this;
  }

  public ComplexDouble conji(@Dimensionless ComplexDouble this) {
    i = -i;
    return this;
  }

  public @Dimensionless ComplexDouble conj(@Dimensionless ComplexDouble this) {
    return dup().conji();
  }

  public @Dimensionless ComplexDouble sqrt(@Dimensionless ComplexDouble this) {
    @Dimensionless
    double a = abs();
    @Dimensionless
    double s2 = (@Dimensionless double) Math.sqrt(((@Dimensionless int) (2)));
    @Dimensionless
    double p = (@Dimensionless double) Math.sqrt(a + r) / s2;
    @Dimensionless
    double sgn = Math.signum(i);
    if (sgn == ((@Dimensionless double) (0.0))) {
      sgn = ((@Dimensionless double) (1.0));
    }
    @Dimensionless
    double q = (@Dimensionless double) Math.sqrt(a - r) / s2 * Math.signum(sgn);
    return new @Dimensionless ComplexDouble(p, q);
  }

  /**
   * Comparing two ComplexDouble values.
   *
   * @param o object to compare this against
   * @return true if both numbers have the same value
   */
  public @Dimensionless boolean equals(@Dimensionless ComplexDouble this, @Dimensionless Object o) {
    if (!(o instanceof ComplexDouble)) {
      return false;
    }
    @Dimensionless
    ComplexDouble c = (@Dimensionless ComplexDouble) o;

    return r == c.r && i == c.i;
  }

  public @Dimensionless int hashCode(@Dimensionless ComplexDouble this) {
    return Double.valueOf(r).hashCode() ^ Double.valueOf(i).hashCode();
  }

  public boolean eq(@Dimensionless ComplexDouble this, ComplexDouble c) {
    return Math.abs(r - c.r) + Math.abs(i - c.i) < (@Dimensionless double) ((@Dimensionless double) (1e-6));
  }

  public @Dimensionless boolean ne(@Dimensionless ComplexDouble this, @Dimensionless ComplexDouble c) {
    return !eq(c);
  }

  public boolean isZero(@Dimensionless ComplexDouble this) {
    return r == ((@Dimensionless double) (0.0)) && i == ((@Dimensionless double) (0.0));
  }

  public @Dimensionless boolean isReal(@Dimensionless ComplexDouble this) {
    return i == ((@Dimensionless double) (0.0));
  }

  public @Dimensionless boolean isImag(@Dimensionless ComplexDouble this) {
    return r == ((@Dimensionless double) (0.0));
  }
}
