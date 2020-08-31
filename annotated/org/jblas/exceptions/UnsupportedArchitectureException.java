package org.jblas.exceptions;
import units.qual.Dimensionless;

/**
 * <one line description>
 * <p/>
 * <longer description>
 * <p/>
 * User: mikio
 * Date: 2/13/13
 * Time: 12:28 PM
 */
public class UnsupportedArchitectureException extends RuntimeException {
  public UnsupportedArchitectureException(@Dimensionless String message) {
    super(message);
  }
}
