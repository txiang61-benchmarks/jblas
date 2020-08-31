package org.jblas;

import units.qual.Dimensionless;
import org.jblas.exceptions.UnsupportedArchitectureException;
import org.jblas.util.LibraryLoader;
import org.jblas.util.Logger;

/**
 * Help class for loading libraries needed for NativeBlas
 *
 * The only use of this class is to have NativeBlas inherit from this class.
 *
 * User: Mikio L. Braun
 * Date: Oct 24, 2012
 */
@Dimensionless
class NativeBlasLibraryLoader {
  static void loadLibraryAndCheckErrors() {
    try {
      try {
        // Try to load it first, probably it's in the path
        System.loadLibrary("jblas");
      } catch (@Dimensionless UnsatisfiedLinkError e) {
        // Nope, ok, so let's copy it.
        Logger.getLogger().config(
            "BLAS native library not found in path. Copying native library "
                + "from the archive. Consider installing the library somewhere "
                + "in the path (for Windows: PATH, for Linux: LD_LIBRARY_PATH).");

        // potentially load dependet libraries (mostly Cygwin libs for Windows)
        loadDependentLibraries();

        // Ok, and now load it!
        new @Dimensionless LibraryLoader().loadLibrary("jblas", true);
      }

      // Let's do some quick tests to see whether we trigger some errors
      // when dependent libraries cannot be found
      @Dimensionless
      double @Dimensionless [] a = new @Dimensionless double @Dimensionless [((@Dimensionless int) (1))];
      NativeBlas.dgemm('N', 'N', ((@Dimensionless int) (1)), ((@Dimensionless int) (1)), ((@Dimensionless int) (1)), ((@Dimensionless double) (1.0)), a, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), a, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)), ((@Dimensionless double) (1.0)), a, ((@Dimensionless int) (0)), ((@Dimensionless int) (1)));
    } catch (@Dimensionless UnsatisfiedLinkError e) {
      @Dimensionless
      String arch = System.getProperty("os.arch");
      @Dimensionless
      String name = System.getProperty("os.name");

      if (name.startsWith("Windows") && e.getMessage().contains("Can't find dependent libraries")) {
        System.err.println("On Windows, you need some additional support libraries.\n" +
          "For example, you can install the two packages in cygwin:\n" +
          "\n" +
          "   mingw64-x86_64-gcc-core   mingw64-x86_64-gfortran\n" +
          "\n" +
          "and add the directory <cygwin-home>\\usr\\x86_64-w64-mingw32\\sys-root\\mingw\\bin to your path.\n\n" +
          "For more information, see http://github.com/mikiobraun/jblas/wiki/Missing-Libraries");
      } else if (name.equals("Linux") && arch.equals("amd64")) {
        System.err.println("On Linux 64bit, you need additional support libraries.\n" +
          "You need to install libgfortran3.\n\n" +
          "For example for debian or Ubuntu, type \"sudo apt-get install libgfortran3\"\n\n" +
          "For more information, see https://github.com/mikiobraun/jblas/wiki/Missing-Libraries");
      }
    } catch (@Dimensionless UnsupportedArchitectureException e) {
      System.err.println(e.getMessage());
    }
  }

  public static void loadDependentLibraries() {
    @Dimensionless
    String arch = System.getProperty("os.arch");
    @Dimensionless
    String name = System.getProperty("os.name");

    @Dimensionless
    LibraryLoader loader = new @Dimensionless LibraryLoader();

    if (name.startsWith("Windows") && arch.equals("amd64")) {
      loader.loadLibrary("libgcc_s_sjlj-1", false);
      loader.loadLibrary("libgfortran-3", false);
    } else if (name.startsWith("Windows") && arch.equals("x86")) {
      loader.loadLibrary("libgcc_s_dw2-1", false);
      loader.loadLibrary("libgfortran-3", false);
    }
  }
}
