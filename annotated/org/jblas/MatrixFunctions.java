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
import units.qual.*;
import units.UnitsTools;

/**
 * This class provides the functions from java.lang.Math for matrices. The
 * functions are applied to each element of the matrix.
 * 
 * @author Mikio Braun
 */
public class MatrixFunctions {

	/*#
	def mapfct(f); <<-EOS
	   for (int i = 0; i < x.length; i++)
	      x.put(i, (double) #{f}(x.get(i)));
	   return x;
	   EOS
  	end
  	
  	def cmapfct(f); <<-EOS
	   for (int i = 0; i < x.length; i++)
	      x.put(i, x.get(i).#{f}());
	   return x;
	   EOS
  	end
	#*/

	/**
	 * Sets all elements in this matrix to their absolute values. Note
	 * that this operation is in-place.
	 * @see MatrixFunctions#abs(DoubleMatrix)
	 * @return this matrix
	 */
	public static DoubleMatrix absi(DoubleMatrix x) { 
		/*# mapfct('Math.abs') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.abs(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	public static @Dimensionless ComplexDoubleMatrix absi(@Dimensionless ComplexDoubleMatrix x) {
		/*# cmapfct('abs') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, x.get(i).abs());
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the trigonometric <i>arccosine</i> function element wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#acos(DoubleMatrix)
	 * @return this matrix
	 */
	public static @Dimensionless DoubleMatrix acosi(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.acos') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@rad double) Math.acos(x.get(i)) / UnitsTools.rad);
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the trigonometric <i>arcsine</i> function element wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#asin(DoubleMatrix)
	 * @return this matrix
	 */	
	public static @Dimensionless DoubleMatrix asini(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.asin') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@rad double) Math.asin(x.get(i)) / UnitsTools.rad);
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the trigonometric <i>arctangend</i> function element wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#atan(DoubleMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless DoubleMatrix atani(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.atan') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@rad double) Math.atan(x.get(i)) / UnitsTools.rad);
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>cube root</i> function element wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#cbrt(DoubleMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless DoubleMatrix cbrti(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.cbrt') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.cbrt(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Element-wise round up by applying the <i>ceil</i> function on each 
	 * element. Note that this is an in-place operation.
	 * @see MatrixFunctions#ceil(DoubleMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless DoubleMatrix ceili(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.ceil') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.ceil(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>cosine</i> function element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#cos(DoubleMatrix)
	 * @return this matrix
	 */
	public static @Dimensionless DoubleMatrix cosi(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.cos') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.cos(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>hyperbolic cosine</i> function element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#cosh(DoubleMatrix)
	 * @return this matrix
	 */	
	public static @Dimensionless DoubleMatrix coshi(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.cosh') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.cosh(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>exponential</i> function element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#exp(DoubleMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless DoubleMatrix expi(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.exp') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.exp(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Element-wise round down by applying the <i>floor</i> function on each 
	 * element. Note that this is an in-place operation.
	 * @see MatrixFunctions#floor(DoubleMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless DoubleMatrix floori(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.floor') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.floor(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>natural logarithm</i> function element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#log(DoubleMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless DoubleMatrix logi(@Dimensionless DoubleMatrix x) {
		/*# mapfct('Math.log') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.log(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>logarithm with basis to 10</i> element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#log10(DoubleMatrix)
	 * @return this matrix
	 */
	public static @Dimensionless DoubleMatrix log10i(@Dimensionless DoubleMatrix x) {
		/*# mapfct('Math.log10') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.log10(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Element-wise power function. Replaces each element with its
	 * power of <tt>d</tt>.Note that this is an in-place operation.
	 * @param d the exponent
	 * @see MatrixFunctions#pow(DoubleMatrix,double)
	 * @return this matrix
	 */	
	public static @Dimensionless DoubleMatrix powi(@Dimensionless DoubleMatrix x, @Dimensionless double d) {
		if (d == ((@Dimensionless double) (2.0)))
			return x.muli(x);
		else {
			for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
				x.put(i, (@Dimensionless double) Math.pow(x.get(i), d));
			return x;
		}
	}

    public static @Dimensionless DoubleMatrix powi(@Dimensionless double base, @Dimensionless DoubleMatrix x) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
            x.put(i, (@Dimensionless double) Math.pow(base, x.get(i)));
        return x;
    }

    public static @Dimensionless DoubleMatrix powi(@Dimensionless DoubleMatrix x, @Dimensionless DoubleMatrix e) {
        x.checkLength(e.length);
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
            x.put(i, (@Dimensionless double) Math.pow(x.get(i), e.get(i)));
        return x;
    }

    public static @Dimensionless DoubleMatrix signumi(@Dimensionless DoubleMatrix x) {
		/*# mapfct('Math.signum') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.signum(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	public static @Dimensionless DoubleMatrix sini(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.sin') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.sin(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}

	public static @Dimensionless DoubleMatrix sinhi(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.sinh') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.sinh(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	public static @Dimensionless DoubleMatrix sqrti(@Dimensionless DoubleMatrix x) { 
		/*# mapfct('Math.sqrt') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.sqrt(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	public static @Dimensionless DoubleMatrix tani(@Dimensionless DoubleMatrix x) {
		/*# mapfct('Math.tan') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.tan(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	public static @Dimensionless DoubleMatrix tanhi(@Dimensionless DoubleMatrix x) {
		/*# mapfct('Math.tanh') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless double) Math.tanh(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}

	/**
	 * Returns a copy of this matrix where all elements are set to their
	 * absolute values. 
	 * @see MatrixFunctions#absi(DoubleMatrix)
	 * @return copy of this matrix
	 */
	public static @Dimensionless DoubleMatrix abs(@Dimensionless DoubleMatrix x) { return absi(x.dup()); }
	
	/**
	 * Returns a copy of this matrix where the trigonometric <i>acos</i> function is applied
	 * element wise.
	 * @see MatrixFunctions#acosi(DoubleMatrix)
	 * @return copy of this matrix
	 */
	public static @Dimensionless DoubleMatrix acos(@Dimensionless DoubleMatrix x)   { return acosi(x.dup()); }
	public static @Dimensionless DoubleMatrix asin(@Dimensionless DoubleMatrix x)   { return asini(x.dup()); }
	public static @Dimensionless DoubleMatrix atan(@Dimensionless DoubleMatrix x)   { return atani(x.dup()); }
	public static @Dimensionless DoubleMatrix cbrt(@Dimensionless DoubleMatrix x)   { return cbrti(x.dup()); }
    public static @Dimensionless DoubleMatrix ceil(@Dimensionless DoubleMatrix x)   { return ceili(x.dup()); }
    public static @Dimensionless DoubleMatrix cos(@Dimensionless DoubleMatrix x)    { return cosi(x.dup()); }
    public static @Dimensionless DoubleMatrix cosh(@Dimensionless DoubleMatrix x)   { return coshi(x.dup()); }
    public static @Dimensionless DoubleMatrix exp(@Dimensionless DoubleMatrix x)    { return expi(x.dup()); }
    public static @Dimensionless DoubleMatrix floor(@Dimensionless DoubleMatrix x)  { return floori(x.dup()); }
    public static @Dimensionless DoubleMatrix log(@Dimensionless DoubleMatrix x)    { return logi(x.dup()); }
    public static @Dimensionless DoubleMatrix log10(@Dimensionless DoubleMatrix x)  { return log10i(x.dup()); }
    public static @Dimensionless double pow(@Dimensionless double x, @Dimensionless double y) { return (@Dimensionless double)Math.pow(x, y); }
    public static @Dimensionless DoubleMatrix pow(@Dimensionless DoubleMatrix x, @Dimensionless double e) { return powi(x.dup(), e); }
    public static @Dimensionless DoubleMatrix pow(@Dimensionless double b, @Dimensionless DoubleMatrix x) { return powi(b, x.dup()); }
    public static @Dimensionless DoubleMatrix pow(@Dimensionless DoubleMatrix x, @Dimensionless DoubleMatrix e) { return powi(x.dup(), e); }
    public static @Dimensionless DoubleMatrix signum(@Dimensionless DoubleMatrix x) { return signumi(x.dup()); }
    public static @Dimensionless DoubleMatrix sin(@Dimensionless DoubleMatrix x)    { return sini(x.dup()); }
    public static @Dimensionless DoubleMatrix sinh(@Dimensionless DoubleMatrix x)   { return sinhi(x.dup()); }
    public static @Dimensionless DoubleMatrix sqrt(@Dimensionless DoubleMatrix x)   { return sqrti(x.dup()); }
    public static @Dimensionless DoubleMatrix tan(@Dimensionless DoubleMatrix x)    { return tani(x.dup()); }
    public static @Dimensionless DoubleMatrix tanh(@Dimensionless DoubleMatrix x)   { return tanhi(x.dup()); }

    /*# %w{abs acos asin atan cbrt ceil cos cosh exp floor log log10 signum sin sinh sqrt tan tanh}.map do |fct| <<-EOS
    public static double #{fct}(double x) { return (double)Math.#{fct}(x); }
    EOS
        end   
     #*/
//RJPP-BEGIN------------------------------------------------------------
    public static @Dimensionless double abs(@Dimensionless double x) { return (@Dimensionless double)Math.abs(x); }
    public static @rad double acos(@Dimensionless double x) { return (@rad double)Math.acos(x); }
    public static @rad double asin(@Dimensionless double x) { return (@rad double)Math.asin(x); }
    public static @rad double atan(@Dimensionless double x) { return (@rad double)Math.atan(x); }
    public static @Dimensionless double cbrt(@Dimensionless double x) { return (@Dimensionless double)Math.cbrt(x); }
    public static @Dimensionless double ceil(@Dimensionless double x) { return (@Dimensionless double)Math.ceil(x); }
    public static @Dimensionless double cos(@rad double x) { return (@Dimensionless double)Math.cos(x); }
    public static @Dimensionless double cosh(@rad double x) { return (@Dimensionless double)Math.cosh(x); }
    public static @Dimensionless double exp(@Dimensionless double x) { return (@Dimensionless double)Math.exp(x); }
    public static @Dimensionless double floor(@Dimensionless double x) { return (@Dimensionless double)Math.floor(x); }
    public static @Dimensionless double log(@Dimensionless double x) { return (@Dimensionless double)Math.log(x); }
    public static @Dimensionless double log10(@Dimensionless double x) { return (@Dimensionless double)Math.log10(x); }
    public static @Dimensionless double signum(@Dimensionless double x) { return (@Dimensionless double)Math.signum(x); }
    public static @Dimensionless double sin(@rad double x) { return (@Dimensionless double)Math.sin(x); }
    public static @Dimensionless double sinh(@rad double x) { return (@Dimensionless double)Math.sinh(x); }
    public static @Dimensionless double sqrt(@Dimensionless double x) { return (@Dimensionless double)Math.sqrt(x); }
    public static @Dimensionless double tan(@rad double x) { return (@Dimensionless double)Math.tan(x); }
    public static @Dimensionless double tanh(@rad double x) { return (@Dimensionless double)Math.tanh(x); }
//RJPP-END--------------------------------------------------------------

    /**
     * Calculate matrix exponential of a square matrix.
     *
     * A scaled Pade approximation algorithm is used.
     * The algorithm has been directly translated from Golub & Van Loan "Matrix Computations",
     * algorithm 11.3.1. Special Horner techniques from 11.2 are also used to minimize the number
     * of matrix multiplications.
     *
     * @param A square matrix
     * @return matrix exponential of A
     */
    public static @Dimensionless DoubleMatrix expm(@Dimensionless DoubleMatrix A) {
        // constants for pade approximation
        final @Dimensionless double c0 = ((@Dimensionless double) (1.0));
        final @Dimensionless double c1 = ((@Dimensionless double) (0.5));
        final @Dimensionless double c2 = ((@Dimensionless double) (0.12));
        final @Dimensionless double c3 = ((@Dimensionless double) (0.01833333333333333));
        final @Dimensionless double c4 = ((@Dimensionless double) (0.0019927536231884053));
        final @Dimensionless double c5 = ((@Dimensionless double) (1.630434782608695E-4));
        final @Dimensionless double c6 = ((@Dimensionless double) (1.0351966873706E-5));
        final @Dimensionless double c7 = ((@Dimensionless double) (5.175983436853E-7));
        final @Dimensionless double c8 = ((@Dimensionless double) (2.0431513566525E-8));
        final @Dimensionless double c9 = ((@Dimensionless double) (6.306022705717593E-10));
        final @Dimensionless double c10 = ((@Dimensionless double) (1.4837700484041396E-11));
        final @Dimensionless double c11 = ((@Dimensionless double) (2.5291534915979653E-13));
        final @Dimensionless double c12 = ((@Dimensionless double) (2.8101705462199615E-15));
        final @Dimensionless double c13 = ((@Dimensionless double) (1.5440497506703084E-17));

        @Dimensionless
        int j = Math.max(((@Dimensionless int) (0)), ((@Dimensionless int) (1)) + (@Dimensionless int) Math.floor(Math.log(A.normmax()) / Math.log(((@Dimensionless int) (2)))));
        @Dimensionless
        DoubleMatrix As = A.div((@Dimensionless double) Math.pow(((@Dimensionless int) (2)), j)); // scaled version of A
        @Dimensionless
        int n = A.getRows();

        // calculate D and N using special Horner techniques
        @Dimensionless
        DoubleMatrix As_2 = As.mmul(As);
        @Dimensionless
        DoubleMatrix As_4 = As_2.mmul(As_2);
        @Dimensionless
        DoubleMatrix As_6 = As_4.mmul(As_2);
        // U = c0*I + c2*A^2 + c4*A^4 + (c6*I + c8*A^2 + c10*A^4 + c12*A^6)*A^6
        @Dimensionless
        DoubleMatrix U = DoubleMatrix.eye(n).muli(c0).addi(As_2.mul(c2)).addi(As_4.mul(c4)).addi(
                DoubleMatrix.eye(n).muli(c6).addi(As_2.mul(c8)).addi(As_4.mul(c10)).addi(As_6.mul(c12)).mmuli(As_6));
        // V = c1*I + c3*A^2 + c5*A^4 + (c7*I + c9*A^2 + c11*A^4 + c13*A^6)*A^6
        @Dimensionless
        DoubleMatrix V = DoubleMatrix.eye(n).muli(c1).addi(As_2.mul(c3)).addi(As_4.mul(c5)).addi(
                DoubleMatrix.eye(n).muli(c7).addi(As_2.mul(c9)).addi(As_4.mul(c11)).addi(As_6.mul(c13)).mmuli(As_6));

        @Dimensionless
        DoubleMatrix AV = As.mmuli(V);
        @Dimensionless
        DoubleMatrix N = U.add(AV);
        @Dimensionless
        DoubleMatrix D = U.subi(AV);

        // solve DF = N for F
        @Dimensionless
        DoubleMatrix F = Solve.solve(D, N);

        // now square j times
        for (@Dimensionless int k = ((@Dimensionless int) (0)); k < j; k++) {
            F.mmuli(F);
        }

        return F;
    }


//STOP
    public static @Dimensionless DoubleMatrix floatToDouble(@Dimensionless FloatMatrix fm) {
    	@Dimensionless
    	DoubleMatrix dm = new @Dimensionless DoubleMatrix(fm.rows, fm.columns);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < fm.length; i++)
            dm.put(i, (@Dimensionless double) fm.get(i));

        return dm;
    }

    public static @Dimensionless FloatMatrix doubleToFloat(@Dimensionless DoubleMatrix dm) {
        @Dimensionless
        FloatMatrix fm = new @Dimensionless FloatMatrix(dm.rows, dm.columns);

        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < dm.length; i++)
            fm.put(i, (@Dimensionless float) dm.get(i));

        return fm;
    }
//START
    
//BEGIN
  // The code below has been automatically generated.
  // DO NOT EDIT!

	/*#
	def mapfct(f); <<-EOS
	   for (int i = 0; i < x.length; i++)
	      x.put(i, (float) #{f}(x.get(i)));
	   return x;
	   EOS
  	end
  	
  	def cmapfct(f); <<-EOS
	   for (int i = 0; i < x.length; i++)
	      x.put(i, x.get(i).#{f}());
	   return x;
	   EOS
  	end
	#*/

	/**
	 * Sets all elements in this matrix to their absolute values. Note
	 * that this operation is in-place.
	 * @see MatrixFunctions#abs(FloatMatrix)
	 * @return this matrix
	 */
	public static FloatMatrix absi(FloatMatrix x) { 
		/*# mapfct('Math.abs') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.abs(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	public static @Dimensionless ComplexFloatMatrix absi(@Dimensionless ComplexFloatMatrix x) {
		/*# cmapfct('abs') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, x.get(i).abs());
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the trigonometric <i>arccosine</i> function element wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#acos(FloatMatrix)
	 * @return this matrix
	 */
	public static @Dimensionless FloatMatrix acosi(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.acos') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@rad float) Math.acos(x.get(i)) / UnitsTools.rad);
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the trigonometric <i>arcsine</i> function element wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#asin(FloatMatrix)
	 * @return this matrix
	 */	
	public static @Dimensionless FloatMatrix asini(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.asin') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@rad float) Math.asin(x.get(i)) / UnitsTools.rad);
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the trigonometric <i>arctangend</i> function element wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#atan(FloatMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless FloatMatrix atani(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.atan') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@rad float) Math.atan(x.get(i)) / UnitsTools.rad);
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>cube root</i> function element wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#cbrt(FloatMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless FloatMatrix cbrti(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.cbrt') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.cbrt(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Element-wise round up by applying the <i>ceil</i> function on each 
	 * element. Note that this is an in-place operation.
	 * @see MatrixFunctions#ceil(FloatMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless FloatMatrix ceili(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.ceil') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.ceil(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>cosine</i> function element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#cos(FloatMatrix)
	 * @return this matrix
	 */
	public static @Dimensionless FloatMatrix cosi(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.cos') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.cos(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>hyperbolic cosine</i> function element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#cosh(FloatMatrix)
	 * @return this matrix
	 */	
	public static @Dimensionless FloatMatrix coshi(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.cosh') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.cosh(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>exponential</i> function element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#exp(FloatMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless FloatMatrix expi(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.exp') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.exp(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Element-wise round down by applying the <i>floor</i> function on each 
	 * element. Note that this is an in-place operation.
	 * @see MatrixFunctions#floor(FloatMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless FloatMatrix floori(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.floor') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.floor(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>natural logarithm</i> function element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#log(FloatMatrix)
	 * @return this matrix
	 */		
	public static @Dimensionless FloatMatrix logi(@Dimensionless FloatMatrix x) {
		/*# mapfct('Math.log') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.log(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Applies the <i>logarithm with basis to 10</i> element-wise on this
	 * matrix. Note that this is an in-place operation.
	 * @see MatrixFunctions#log10(FloatMatrix)
	 * @return this matrix
	 */
	public static @Dimensionless FloatMatrix log10i(@Dimensionless FloatMatrix x) {
		/*# mapfct('Math.log10') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.log10(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	/**
	 * Element-wise power function. Replaces each element with its
	 * power of <tt>d</tt>.Note that this is an in-place operation.
	 * @param d the exponent
	 * @see MatrixFunctions#pow(FloatMatrix,float)
	 * @return this matrix
	 */	
	public static @Dimensionless FloatMatrix powi(@Dimensionless FloatMatrix x, @Dimensionless float d) {
		if (d == ((@Dimensionless float) (2.0f)))
			return x.muli(x);
		else {
			for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
				x.put(i, (@Dimensionless float) Math.pow(x.get(i), d));
			return x;
		}
	}

    public static @Dimensionless FloatMatrix powi(@Dimensionless float base, @Dimensionless FloatMatrix x) {
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
            x.put(i, (@Dimensionless float) Math.pow(base, x.get(i)));
        return x;
    }

    public static @Dimensionless FloatMatrix powi(@Dimensionless FloatMatrix x, @Dimensionless FloatMatrix e) {
        x.checkLength(e.length);
        for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
            x.put(i, (@Dimensionless float) Math.pow(x.get(i), e.get(i)));
        return x;
    }

    public static @Dimensionless FloatMatrix signumi(@Dimensionless FloatMatrix x) {
		/*# mapfct('Math.signum') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.signum(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	
	public static @Dimensionless FloatMatrix sini(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.sin') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.sin(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}

	public static @Dimensionless FloatMatrix sinhi(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.sinh') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.sinh(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	public static @Dimensionless FloatMatrix sqrti(@Dimensionless FloatMatrix x) { 
		/*# mapfct('Math.sqrt') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.sqrt(x.get(i)));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	public static @Dimensionless FloatMatrix tani(@Dimensionless FloatMatrix x) {
		/*# mapfct('Math.tan') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.tan(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}
	public static @Dimensionless FloatMatrix tanhi(@Dimensionless FloatMatrix x) {
		/*# mapfct('Math.tanh') #*/
//RJPP-BEGIN------------------------------------------------------------
	   for (@Dimensionless int i = ((@Dimensionless int) (0)); i < x.length; i++)
	      x.put(i, (@Dimensionless float) Math.tanh(x.get(i) * UnitsTools.rad));
	   return x;
//RJPP-END--------------------------------------------------------------
	}

	/**
	 * Returns a copy of this matrix where all elements are set to their
	 * absolute values. 
	 * @see MatrixFunctions#absi(FloatMatrix)
	 * @return copy of this matrix
	 */
	public static @Dimensionless FloatMatrix abs(@Dimensionless FloatMatrix x) { return absi(x.dup()); }
	
	/**
	 * Returns a copy of this matrix where the trigonometric <i>acos</i> function is applied
	 * element wise.
	 * @see MatrixFunctions#acosi(FloatMatrix)
	 * @return copy of this matrix
	 */
	public static @Dimensionless FloatMatrix acos(@Dimensionless FloatMatrix x)   { return acosi(x.dup()); }
	public static @Dimensionless FloatMatrix asin(@Dimensionless FloatMatrix x)   { return asini(x.dup()); }
	public static @Dimensionless FloatMatrix atan(@Dimensionless FloatMatrix x)   { return atani(x.dup()); }
	public static @Dimensionless FloatMatrix cbrt(@Dimensionless FloatMatrix x)   { return cbrti(x.dup()); }
    public static @Dimensionless FloatMatrix ceil(@Dimensionless FloatMatrix x)   { return ceili(x.dup()); }
    public static @Dimensionless FloatMatrix cos(@Dimensionless FloatMatrix x)    { return cosi(x.dup()); }
    public static @Dimensionless FloatMatrix cosh(@Dimensionless FloatMatrix x)   { return coshi(x.dup()); }
    public static @Dimensionless FloatMatrix exp(@Dimensionless FloatMatrix x)    { return expi(x.dup()); }
    public static @Dimensionless FloatMatrix floor(@Dimensionless FloatMatrix x)  { return floori(x.dup()); }
    public static @Dimensionless FloatMatrix log(@Dimensionless FloatMatrix x)    { return logi(x.dup()); }
    public static @Dimensionless FloatMatrix log10(@Dimensionless FloatMatrix x)  { return log10i(x.dup()); }
    public static @Dimensionless float pow(@Dimensionless float x, @Dimensionless float y) { return (@Dimensionless float)Math.pow(x, y); }
    public static @Dimensionless FloatMatrix pow(@Dimensionless FloatMatrix x, @Dimensionless float e) { return powi(x.dup(), e); }
    public static @Dimensionless FloatMatrix pow(@Dimensionless float b, @Dimensionless FloatMatrix x) { return powi(b, x.dup()); }
    public static @Dimensionless FloatMatrix pow(@Dimensionless FloatMatrix x, @Dimensionless FloatMatrix e) { return powi(x.dup(), e); }
    public static @Dimensionless FloatMatrix signum(@Dimensionless FloatMatrix x) { return signumi(x.dup()); }
    public static @Dimensionless FloatMatrix sin(@Dimensionless FloatMatrix x)    { return sini(x.dup()); }
    public static @Dimensionless FloatMatrix sinh(@Dimensionless FloatMatrix x)   { return sinhi(x.dup()); }
    public static @Dimensionless FloatMatrix sqrt(@Dimensionless FloatMatrix x)   { return sqrti(x.dup()); }
    public static @Dimensionless FloatMatrix tan(@Dimensionless FloatMatrix x)    { return tani(x.dup()); }
    public static @Dimensionless FloatMatrix tanh(@Dimensionless FloatMatrix x)   { return tanhi(x.dup()); }

    /*# %w{abs acos asin atan cbrt ceil cos cosh exp floor log log10 signum sin sinh sqrt tan tanh}.map do |fct| <<-EOS
    public static float #{fct}(float x) { return (float)Math.#{fct}(x); }
    EOS
        end   
     #*/
//RJPP-BEGIN------------------------------------------------------------
    public static @Dimensionless float abs(@Dimensionless float x) { return (@Dimensionless float)Math.abs(x); }
    public static @rad float acos(@Dimensionless float x) { return (@rad float)Math.acos(x); }
    public static @rad float asin(@Dimensionless float x) { return (@rad float)Math.asin(x); }
    public static @rad float atan(@Dimensionless float x) { return (@rad float)Math.atan(x); }
    public static @Dimensionless float cbrt(@Dimensionless float x) { return (@Dimensionless float)Math.cbrt(x); }
    public static @Dimensionless float ceil(@Dimensionless float x) { return (@Dimensionless float)Math.ceil(x); }
    public static @Dimensionless float cos(@rad float x) { return (@Dimensionless float)Math.cos(x); }
    public static @Dimensionless float cosh(@rad float x) { return (@Dimensionless float)Math.cosh(x); }
    public static @Dimensionless float exp(@Dimensionless float x) { return (@Dimensionless float)Math.exp(x); }
    public static @Dimensionless float floor(@Dimensionless float x) { return (@Dimensionless float)Math.floor(x); }
    public static @Dimensionless float log(@Dimensionless float x) { return (@Dimensionless float)Math.log(x); }
    public static @Dimensionless float log10(@Dimensionless float x) { return (@Dimensionless float)Math.log10(x); }
    public static @Dimensionless float signum(@Dimensionless float x) { return (@Dimensionless float)Math.signum(x); }
    public static @Dimensionless float sin(@rad float x) { return (@Dimensionless float)Math.sin(x); }
    public static @Dimensionless float sinh(@rad float x) { return (@Dimensionless float)Math.sinh(x); }
    public static @Dimensionless float sqrt(@Dimensionless float x) { return (@Dimensionless float)Math.sqrt(x); }
    public static @Dimensionless float tan(@rad float x) { return (@Dimensionless float)Math.tan(x); }
    public static @Dimensionless float tanh(@rad float x) { return (@Dimensionless float)Math.tanh(x); }
//RJPP-END--------------------------------------------------------------

    /**
     * Calculate matrix exponential of a square matrix.
     *
     * A scaled Pade approximation algorithm is used.
     * The algorithm has been directly translated from Golub & Van Loan "Matrix Computations",
     * algorithm 11.3f.1. Special Horner techniques from 11.2f are also used to minimize the number
     * of matrix multiplications.
     *
     * @param A square matrix
     * @return matrix exponential of A
     */
    public static @Dimensionless FloatMatrix expm(@Dimensionless FloatMatrix A) {
        // constants for pade approximation
        final @Dimensionless float c0 = ((@Dimensionless float) (1.0f));
        final @Dimensionless float c1 = ((@Dimensionless float) (0.5f));
        final @Dimensionless float c2 = ((@Dimensionless float) (0.12f));
        final @Dimensionless float c3 = ((@Dimensionless float) (0.01833333333333333f));
        final @Dimensionless float c4 = ((@Dimensionless float) (0.0019927536231884053f));
        final @Dimensionless float c5 = ((@Dimensionless float) (1.630434782608695E-4f));
        final @Dimensionless float c6 = ((@Dimensionless float) (1.0351966873706E-5f));
        final @Dimensionless float c7 = ((@Dimensionless float) (5.175983436853E-7f));
        final @Dimensionless float c8 = ((@Dimensionless float) (2.0431513566525E-8f));
        final @Dimensionless float c9 = ((@Dimensionless float) (6.306022705717593E-10f));
        final @Dimensionless float c10 = ((@Dimensionless float) (1.4837700484041396E-11f));
        final @Dimensionless float c11 = ((@Dimensionless float) (2.5291534915979653E-13f));
        final @Dimensionless float c12 = ((@Dimensionless float) (2.8101705462199615E-15f));
        final @Dimensionless float c13 = ((@Dimensionless float) (1.5440497506703084E-17f));

        @Dimensionless
        int j = Math.max(((@Dimensionless int) (0)), ((@Dimensionless int) (1)) + (@Dimensionless int) Math.floor(Math.log(A.normmax()) / Math.log(((@Dimensionless int) (2)))));
        @Dimensionless
        FloatMatrix As = A.div((@Dimensionless float) Math.pow(((@Dimensionless int) (2)), j)); // scaled version of A
        @Dimensionless
        int n = A.getRows();

        // calculate D and N using special Horner techniques
        @Dimensionless
        FloatMatrix As_2 = As.mmul(As);
        @Dimensionless
        FloatMatrix As_4 = As_2.mmul(As_2);
        @Dimensionless
        FloatMatrix As_6 = As_4.mmul(As_2);
        // U = c0*I + c2*A^2 + c4*A^4 + (c6*I + c8*A^2 + c10*A^4 + c12*A^6)*A^6
        @Dimensionless
        FloatMatrix U = FloatMatrix.eye(n).muli(c0).addi(As_2.mul(c2)).addi(As_4.mul(c4)).addi(
                FloatMatrix.eye(n).muli(c6).addi(As_2.mul(c8)).addi(As_4.mul(c10)).addi(As_6.mul(c12)).mmuli(As_6));
        // V = c1*I + c3*A^2 + c5*A^4 + (c7*I + c9*A^2 + c11*A^4 + c13*A^6)*A^6
        @Dimensionless
        FloatMatrix V = FloatMatrix.eye(n).muli(c1).addi(As_2.mul(c3)).addi(As_4.mul(c5)).addi(
                FloatMatrix.eye(n).muli(c7).addi(As_2.mul(c9)).addi(As_4.mul(c11)).addi(As_6.mul(c13)).mmuli(As_6));

        @Dimensionless
        FloatMatrix AV = As.mmuli(V);
        @Dimensionless
        FloatMatrix N = U.add(AV);
        @Dimensionless
        FloatMatrix D = U.subi(AV);

        // solve DF = N for F
        @Dimensionless
        FloatMatrix F = Solve.solve(D, N);

        // now square j times
        for (@Dimensionless int k = ((@Dimensionless int) (0)); k < j; k++) {
            F.mmuli(F);
        }

        return F;
    }


    
//END
}
