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
 * A range over all available indices. Can be used to address whole columns or rows. Like
 * the ":" index in matlab. Don't forget to call init() before using this range.
 */
@Dimensionless
public class AllRange implements Range {
  private @Dimensionless int lower;
  private @Dimensionless int upper;
  private @Dimensionless int value;
  private @Dimensionless int counter;

  public AllRange() {
  }

  public void init(@Dimensionless AllRange this, @Dimensionless int l, @Dimensionless int u) {
    lower = l;
    upper = u;
    value = l;
    counter = ((@Dimensionless int) (0));
  }

  public @Dimensionless int length(@Dimensionless AllRange this) {
    return upper - lower;
  }

  public @Dimensionless int value(@Dimensionless AllRange this) {
    return value;
  }

  public @Dimensionless int index(@Dimensionless AllRange this) {
    return counter;
  }

  public void next(@Dimensionless AllRange this) {
    counter++;
    value++;
  }

  public @Dimensionless boolean hasMore(@Dimensionless AllRange this) {
    return value < upper;
  }

  @Override
  public @Dimensionless String toString(@Dimensionless AllRange this) {
    return String.format("<AllRange from %d to %d, with length %d, index=%d, value=%d>", lower, upper, length(), index(), value());
  }
}
