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
 * A complex value with float precision.
 *
 * @author Mikio L. Braun
 */
public class ComplexFloat {

  private float r, i;
  public static final ComplexFloat UNIT = new @Dimensionless ComplexFloat(((@Dimensionless float) (1.0f)), ((@Dimensionless float) (0.0f)));
  public static final @Dimensionless ComplexFloat I = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)), ((@Dimensionless float) (1.0f)));
  public static final ComplexFloat NEG_UNIT = new @Dimensionless ComplexFloat(- ((@Dimensionless float) (1.0f)), ((@Dimensionless float) (0.0f)));
  public static final @Dimensionless ComplexFloat NEG_I = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)), - ((@Dimensionless float) (1.0f)));
  public static final ComplexFloat ZERO = new @Dimensionless ComplexFloat(((@Dimensionless float) (0.0f)));

  public ComplexFloat(float real, float imag) {
    r = real;
    i = imag;
  }

  public ComplexFloat(float real) {
    this(real, ((@Dimensionless float) (0.0f)));
  }

  public @Dimensionless String toString(@Dimensionless ComplexFloat this) {
    if (i >= ((@Dimensionless int) (0))) {
      return r + " + " + i + "i";
    } else {
      return r + " - " + (-i) + "i";
    }
  }

  public ComplexFloat set(@Dimensionless ComplexFloat this, float real, float imag) {
    r = real;
    i = imag;
    return this;
  }

  public float real(@Dimensionless ComplexFloat this) {
    return r;
  }

  public float imag(@Dimensionless ComplexFloat this) {
    return i;
  }

  public @Dimensionless ComplexFloat dup(@Dimensionless ComplexFloat this) {
    return new @Dimensionless ComplexFloat(r, i);
  }

  public ComplexFloat copy(@Dimensionless ComplexFloat this, ComplexFloat other) {
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
  public @Dimensionless ComplexFloat addi(@Dimensionless ComplexFloat this, @Dimensionless ComplexFloat c, @Dimensionless ComplexFloat result) {
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
  public ComplexFloat addi(@Dimensionless ComplexFloat this, ComplexFloat c) {
    return addi(c, this);
  }

  /**
   * Add two complex numbers.
   *
   * @param c other complex number
   * @return new complex number with result
   */
  public ComplexFloat add(@Dimensionless ComplexFloat this, ComplexFloat c) {
    return dup().addi(c);
  }

  /**
   * Add a real number to a complex number in-place.
   *
   * @param a      real number to add
   * @param result complex number to hold result
   * @return same as result
   */
  public @Dimensionless ComplexFloat addi(@Dimensionless ComplexFloat this, @Dimensionless float a, @Dimensionless ComplexFloat result) {
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
  public @Dimensionless ComplexFloat addi(@Dimensionless ComplexFloat this, @Dimensionless float c) {
    return addi(c, this);
  }

  /**
   * Add a real number to a complex number.
   *
   * @param c real number to add
   * @return new complex number with result
   */
  public @Dimensionless ComplexFloat add(@Dimensionless ComplexFloat this, @Dimensionless float c) {
    return dup().addi(c);
  }

  /**
   * Subtract two complex numbers, in-place
   *
   * @param c      complex number to subtract
   * @param result resulting complex number
   * @return same as result
   */
  public @Dimensionless ComplexFloat subi(@Dimensionless ComplexFloat this, @Dimensionless ComplexFloat c, @Dimensionless ComplexFloat result) {
    if (this == result) {
      r -= c.r;
      i -= c.i;
    } else {
      result.r = r - c.r;
      result.i = i - c.i;
    }
    return this;
  }

  public @Dimensionless ComplexFloat subi(@Dimensionless ComplexFloat this, @Dimensionless ComplexFloat c) {
    return subi(c, this);
  }

  /**
   * Subtract two complex numbers
   *
   * @param c complex number to subtract
   * @return new complex number with result
   */
  public ComplexFloat sub(@Dimensionless ComplexFloat this, ComplexFloat c) {
    return dup().subi(c);
  }

  public @Dimensionless ComplexFloat subi(@Dimensionless ComplexFloat this, @Dimensionless float a, @Dimensionless ComplexFloat result) {
    if (this == result) {
      r -= a;
    } else {
      result.r = r - a;
      result.i = i;
    }
    return result;
  }

  public @Dimensionless ComplexFloat subi(@Dimensionless ComplexFloat this, @Dimensionless float a) {
    return subi(a, this);
  }

  public @Dimensionless ComplexFloat sub(@Dimensionless ComplexFloat this, @Dimensionless float r) {
    return dup().subi(r);
  }

  /**
   * Multiply two complex numbers, in-place
   *
   * @param c      other complex number
   * @param result complex number where product is stored
   * @return same as result
   */
  public @Dimensionless ComplexFloat muli(@Dimensionless ComplexFloat this, @Dimensionless ComplexFloat c, @Dimensionless ComplexFloat result) {
    @Dimensionless
    float newR = r * c.r - i * c.i;
    @Dimensionless
    float newI = r * c.i + i * c.r;
    result.r = newR;
    result.i = newI;
    return result;
  }

  public ComplexFloat muli(@Dimensionless ComplexFloat this, ComplexFloat c) {
    return muli(c, this);
  }

  /**
   * Multiply two complex numbers
   *
   * @param c other complex number
   * @return new complex number with product of this and c
   */
  public @Dimensionless ComplexFloat mul(@Dimensionless ComplexFloat this, @Dimensionless ComplexFloat c) {
    return dup().muli(c);
  }

  public @Dimensionless ComplexFloat mul(@Dimensionless ComplexFloat this, @Dimensionless float v) {
    return dup().muli(v);
  }

  public @Dimensionless ComplexFloat muli(@Dimensionless ComplexFloat this, @Dimensionless float v, @Dimensionless ComplexFloat result) {
    if (this == result) {
      r *= v;
      i *= v;
    } else {
      result.r = r * v;
      result.i = i * v;
    }
    return this;
  }

  public @Dimensionless ComplexFloat muli(@Dimensionless ComplexFloat this, @Dimensionless float v) {
    return muli(v, this);
  }

  /**
   * Divide two complex numbers
   *
   * @param c complex number to divide this by
   * @return new complex number with quotient of this and c
   */
  public @Dimensionless ComplexFloat div(@Dimensionless ComplexFloat this, @Dimensionless ComplexFloat c) {
    return dup().divi(c);
  }

  /**
   * Divide two complex numbers, in-place
   *
   * @param c      complex number to divide this by
   * @param result complex number to hold result
   * @return same as result
   */
  public @Dimensionless ComplexFloat divi(@Dimensionless ComplexFloat this, @Dimensionless ComplexFloat c, @Dimensionless ComplexFloat result) {
    @Dimensionless
    float d = c.r * c.r + c.i * c.i;
    @Dimensionless
    float newR = (r * c.r + i * c.i) / d;
    @Dimensionless
    float newI = (i * c.r - r * c.i) / d;
    result.r = newR;
    result.i = newI;
    return result;
  }

  public ComplexFloat divi(@Dimensionless ComplexFloat this, ComplexFloat c) {
    return divi(c, this);
  }

  public @Dimensionless ComplexFloat divi(@Dimensionless ComplexFloat this, @Dimensionless float v, @Dimensionless ComplexFloat result) {
    if (this == result) {
      r /= v;
      i /= v;
    } else {
      result.r = r / v;
      result.i = i / v;
    }
    return this;
  }

  public @Dimensionless ComplexFloat divi(@Dimensionless ComplexFloat this, @Dimensionless float v) {
    return divi(v, this);
  }

  public ComplexFloat div(@Dimensionless ComplexFloat this, float v) {
    return dup().divi(v);
  }

  /**
   * Return the absolute value
   *
   * @return the result (length of the vector in 2d plane)
   */
  public float abs(@Dimensionless ComplexFloat this) {
    return (@Dimensionless float) Math.sqrt(r * r + i * i);
  }

  /**
   * Returns the argument of a complex number.
   *
   * @return the result (angle in radians of the vector in 2d plane)
   */
  public @rad float arg(@Dimensionless ComplexFloat this) {
    return (@rad float) Math.atan2(i, r);
  }

  public @Dimensionless ComplexFloat invi(@Dimensionless ComplexFloat this) {
    @Dimensionless
    float d = r * r + i * i;
    r = r / d;
    i = -i / d;
    return this;
  }

  public @Dimensionless ComplexFloat inv(@Dimensionless ComplexFloat this) {
    return dup().invi();
  }

  public @Dimensionless ComplexFloat neg(@Dimensionless ComplexFloat this) {
    return dup().negi();
  }

  public ComplexFloat negi(@Dimensionless ComplexFloat this) {
    r = -r;
    i = -i;
    return this;
  }

  public ComplexFloat conji(@Dimensionless ComplexFloat this) {
    i = -i;
    return this;
  }

  public @Dimensionless ComplexFloat conj(@Dimensionless ComplexFloat this) {
    return dup().conji();
  }

  public @Dimensionless ComplexFloat sqrt(@Dimensionless ComplexFloat this) {
    @Dimensionless
    float a = abs();
    @Dimensionless
    float s2 = (@Dimensionless float) Math.sqrt(((@Dimensionless int) (2)));
    @Dimensionless
    float p = (@Dimensionless float) Math.sqrt(a + r) / s2;
    @Dimensionless
    float sgn = Math.signum(i);
    if (sgn == ((@Dimensionless float) (0.0f))) {
      sgn = ((@Dimensionless float) (1.0f));
    }
    @Dimensionless
    float q = (@Dimensionless float) Math.sqrt(a - r) / s2 * Math.signum(sgn);
    return new @Dimensionless ComplexFloat(p, q);
  }

  /**
   * Comparing two ComplexFloat values.
   *
   * @param o object to compare this against
   * @return true if both numbers have the same value
   */
  public @Dimensionless boolean equals(@Dimensionless ComplexFloat this, @Dimensionless Object o) {
    if (!(o instanceof ComplexFloat)) {
      return false;
    }
    @Dimensionless
    ComplexFloat c = (@Dimensionless ComplexFloat) o;

    return r == c.r && i == c.i;
  }

  public @Dimensionless int hashCode(@Dimensionless ComplexFloat this) {
    return Float.valueOf(r).hashCode() ^ Float.valueOf(i).hashCode();
  }

  public boolean eq(@Dimensionless ComplexFloat this, ComplexFloat c) {
    return Math.abs(r - c.r) + Math.abs(i - c.i) < (@Dimensionless float) ((@Dimensionless double) (1e-6));
  }

  public @Dimensionless boolean ne(@Dimensionless ComplexFloat this, @Dimensionless ComplexFloat c) {
    return !eq(c);
  }

  public boolean isZero(@Dimensionless ComplexFloat this) {
    return r == ((@Dimensionless float) (0.0f)) && i == ((@Dimensionless float) (0.0f));
  }

  public @Dimensionless boolean isReal(@Dimensionless ComplexFloat this) {
    return i == ((@Dimensionless float) (0.0f));
  }

  public @Dimensionless boolean isImag(@Dimensionless ComplexFloat this) {
    return r == ((@Dimensionless float) (0.0f));
  }
}