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
 * A PointRange is a range which only has a single point.
 */
@Dimensionless
public class PointRange implements Range {
  private @Dimensionless int value;
  private @Dimensionless boolean consumed;

  /**
   * Construct a new PointRange with the one given index.
   */
  public PointRange(@Dimensionless int v) {
    value = v;
  }

  public void init(@Dimensionless PointRange this, @Dimensionless int l, @Dimensionless int u) {
    consumed = false;
  }

  public @Dimensionless int length(@Dimensionless PointRange this) {
    return ((@Dimensionless int) (1));
  }

  public @Dimensionless int value(@Dimensionless PointRange this) {
    return value;
  }

  public @Dimensionless int index(@Dimensionless PointRange this) {
    return ((@Dimensionless int) (0));
  }

  public void next(@Dimensionless PointRange this) {
    consumed = true;
  }

  public @Dimensionless boolean hasMore(@Dimensionless PointRange this) {
    return !consumed;
  }

  @Override
  public @Dimensionless String toString(@Dimensionless PointRange this) {
    return String.format("<PointRange at=%d>", value);
  }
}
