package org.jblas.util;
import units.qual.Dimensionless;

/**
 * Created by IntelliJ IDEA.
 * User: mikio
 * Date: 6/24/11
 * Time: 10:45 AM
 * To change this template use File | Settings | File Templates.
 */

public class Random {
    private static java.util.@Dimensionless Random r = new java.util.@Dimensionless Random();

    public static void seed(@Dimensionless long newSeed) {
        r = new java.util.@Dimensionless Random(newSeed);
    }

    public static double nextDouble() {
        return r.nextDouble();
    }

    public static float nextFloat() {
        return r.nextFloat();
    }

    public static @Dimensionless int nextInt(@Dimensionless int max) {
        return r.nextInt(max);
    }

    public static double nextGaussian() {
        return r.nextGaussian();
    }
}
