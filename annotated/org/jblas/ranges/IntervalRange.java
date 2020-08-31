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
package org.jblas.ranges;
import units.qual.Dimensionless;

/**
 * Range which varies from a given interval. Endpoint is exclusive!
 * <p/>
 * "new IntervalRange(0, 3)" enumerates 0, 1, 2.
 */
public class IntervalRange implements Range {
  private @Dimensionless int start;
  private @Dimensionless int end;
  private @Dimensionless int value;

  /**
   * Construct new interval range. Endpoints are inclusive.
   */
  public IntervalRange(int a, int b) {
    start = a;
    end = b;
  }

  public void init(@Dimensionless IntervalRange this, @Dimensionless int lower, @Dimensionless int upper) {
    value = start;
    if (start < lower || end > upper + ((@Dimensionless int) (1))) {
      throw new @Dimensionless IllegalArgumentException("Bounds " + lower + " to " + upper + " are beyond range interval " + start + " to " + end + ".");
    }
  }

  public @Dimensionless int length(@Dimensionless IntervalRange this) {
    return end - start;
  }

  public void next(@Dimensionless IntervalRange this) {
    value++;
  }

  public @Dimensionless int index(@Dimensionless IntervalRange this) {
    return value - start;
  }

  public @Dimensionless int value(@Dimensionless IntervalRange this) {
    return value;
  }

  public @Dimensionless boolean hasMore(@Dimensionless IntervalRange this) {
    return value < end;
  }

  @Override
  public @Dimensionless String toString(@Dimensionless IntervalRange this) {
    return String.format("<Interval Range from %d to %d, length %d index=%d value=%d>", start, end, length(), index(), value());
  }
}
