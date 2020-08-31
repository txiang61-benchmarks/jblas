/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.jblas.exceptions;
import units.qual.Dimensionless;

/**
 *
 */
@Dimensionless
public class LapackPositivityException extends LapackException {
    public LapackPositivityException(@Dimensionless String fct, @Dimensionless String msg) {
        super(fct, msg);
    }
}
