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

/**
 *
 */
@Dimensionless
public class Logger {
    public static final @Dimensionless int ERROR = ((@Dimensionless int) (5));
    public static final @Dimensionless int WARNING = ((@Dimensionless int) (4));
    public static final @Dimensionless int INFO = ((@Dimensionless int) (3));
    public static final @Dimensionless int CONFIG = ((@Dimensionless int) (2));
    public static final @Dimensionless int DEBUG = ((@Dimensionless int) (1));

    public static final @Dimensionless String levelNames @Dimensionless [] = new @Dimensionless String @Dimensionless [] {
        "DEBUG", "CONFIG", "INFO", "WARNING", "ERROR"
    };

    private static @Dimensionless Logger theLogger = new @Dimensionless Logger();
    private @Dimensionless int level;

    private Logger() {
        level = ((@Dimensionless int) (INFO));
    }

    public static @Dimensionless Logger getLogger() {
        return theLogger;
    }

    public void log(@Dimensionless Logger this, @Dimensionless int messageLevel, @Dimensionless String msg) {
        if (level <= messageLevel) {
            System.err.println("-- org.jblas " + levelNames[messageLevel - ((@Dimensionless int) (1))] + " "+ msg);
        }
    }

    public void debug(@Dimensionless Logger this, @Dimensionless String msg) {
        log(((@Dimensionless int) (DEBUG)), msg);
    }

    public void config(@Dimensionless Logger this, @Dimensionless String msg) {
        log(((@Dimensionless int) (CONFIG)), msg);
    }

    public void info(@Dimensionless Logger this, @Dimensionless String msg) {
        log(((@Dimensionless int) (INFO)), msg);
    }

    public void warning(@Dimensionless Logger this, @Dimensionless String msg) {
        log(((@Dimensionless int) (WARNING)), msg);
    }

    public void error(@Dimensionless Logger this, @Dimensionless String msg) {
        log(((@Dimensionless int) (ERROR)), msg);
    }

    public void setLevel(@Dimensionless Logger this, @Dimensionless int level) {
        this.level = level;
    }
}
