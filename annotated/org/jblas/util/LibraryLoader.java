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
import units.qual.ms;
import org.jblas.exceptions.UnsupportedArchitectureException;

import java.io.*;

/**
 * Class which allows to load a dynamic file as resource (for example, from a
 * jar-file)
 */
public class LibraryLoader {

  private @Dimensionless Logger logger;
  private @Dimensionless String libpath;

  private static @Dimensionless File tempDir;

  static {
    final @Dimensionless Logger logger = Logger.getLogger();

    try {
      tempDir = File.createTempFile("jblas", "");

      if (!tempDir.delete() || !tempDir.mkdir()) {
        throw new @Dimensionless IOException(String.format("Couldn't create directory \"%s\"", tempDir.getAbsolutePath()));
      }

      /*
       * Different cleanup strategies for Windows and Linux.
       *
       * For *NIX operating systems: A shutdown hook to clean up the files created. Under
       * Windows this won't work because 
       */
      if (getUnifiedOSName() != "Windows") {
        Runtime.getRuntime().addShutdownHook(new @Dimensionless Thread() @Dimensionless {
          @Override
          public void run() {
            for (@Dimensionless File f : tempDir.listFiles()) {
              logger.info("Deleting " + f.getAbsolutePath());
              if (!f.delete()) {
                logger.warning(String.format("Couldn't delete temporary file \"%s\"", f.getAbsolutePath()));
              }
            }
            logger.info("Deleting " + tempDir.getAbsolutePath());
            if (!tempDir.delete()) {
              logger.warning(String.format("Couldn't delete temporary directory \"%s\"", tempDir.getAbsolutePath()));
            }
          }
        });
      } else {
        new @Dimensionless Thread() @Dimensionless {
          @Override
          public void run() {
            try {
              Thread.sleep(((@ms int) (1000)));

              logger.info("Starting temp DLL cleanup task.");

              @Dimensionless
              int deletedFiles = ((@Dimensionless int) (0));

              @Dimensionless
              File jblasTempDir = new @Dimensionless File(System.getProperty("java.io.tmpdir"));
              for (@Dimensionless File jblasDir : jblasTempDir.listFiles()) {
                assert (jblasDir != null);
                if (jblasDir != tempDir && jblasDir.isDirectory() && jblasDir.getName().startsWith("jblas")) {
                  for (@Dimensionless File oldJblasFile : jblasDir.listFiles()) {
                    if (!oldJblasFile.delete()) {
                      logger.debug("Couldn't delete " + oldJblasFile.getAbsolutePath());
                    } else {
                      logger.debug("Deleted " + oldJblasFile.getAbsolutePath());
                      deletedFiles++;
                    }
                  }
                }
              }

              if (deletedFiles > ((@Dimensionless int) (0))) {
                logger.info(String.format("Deleted %d unused temp DLL libraries from %s", deletedFiles, jblasTempDir.getAbsolutePath()));
              }
            } catch (@Dimensionless InterruptedException ex) {
              //
            }
          }
        }.start();
      }
    } catch (@Dimensionless IOException ex) {
      logger.error("Couldn't create temporary directory: " + ex.getMessage());
    }
  }

  public LibraryLoader() {
    logger = Logger.getLogger();
    libpath = null;
  }

  /**
   * <p>Find the library <tt>libname</tt> as a resource, copy it to a tempfile
   * and load it using System.load(). The name of the library has to be the
   * base name, it is mapped to the corresponding system name using
   * System.mapLibraryName(). For example, the library "foo" is called "libfoo.so"
   * under Linux and "foo.dll" under Windows, but you just have to pass "foo"
   * the loadLibrary().</p>
   * <p/>
   * <p>I'm not quite sure if this doesn't open all kinds of security holes. Any ideas?</p>
   * <p/>
   * <p>This function reports some more information to the "org.jblas" logger at
   * the FINE level.</p>
   *
   * @param libname basename of the library
   * @throws UnsatisfiedLinkError if library cannot be founds
   */
  public void loadLibrary(@Dimensionless LibraryLoader this, String libname, boolean withFlavor) {
    // preload flavor libraries
    @Dimensionless
    String flavor = null;
    if (withFlavor) {
      logger.debug("Preloading ArchFlavor library.");
      flavor = ArchFlavor.archFlavor();
      if (flavor != null && flavor.equals("sse2")) {
        throw new @Dimensionless UnsupportedArchitectureException("Support for SSE2 processors stopped with version 1.2.2. Sorry.");
      }
    }
    logger.debug("Found flavor = '" + flavor + "'");

    libname = System.mapLibraryName(libname);

    /*
     * JDK 7 changed the ending for Mac OS from "jnilib" to "dylib".
     *
     * If that is the case, remap the filename.
     */
    @Dimensionless
    String loadLibname = libname;
    if (libname.endsWith("dylib")) {
      loadLibname = libname.replace(".dylib", ".jnilib");
      logger.config("Replaced .dylib with .jnilib");
    }

    logger.debug("Attempting to load \"" + loadLibname + "\".");

    @Dimensionless
    String @Dimensionless [] paths = new @Dimensionless String @Dimensionless [] {
        fatJarLibraryPath("static", flavor),
        fatJarLibraryPath("dynamic", flavor),
    };

    @Dimensionless
    InputStream is = findLibrary(paths, loadLibname);

    // Haven't found the lib anywhere? Throw a reception.
    if (is == null) {
      throw new @Dimensionless UnsatisfiedLinkError("Couldn't find the resource " + loadLibname + ".");
    }

    logger.config("Loading " + loadLibname + " from " + libpath + ", copying to " + libname + ".");
    loadLibraryFromStream(libname, is);
  }

  private @Dimensionless InputStream findLibrary(@Dimensionless LibraryLoader this, @Dimensionless String @Dimensionless [] paths, @Dimensionless String libname) {
    @Dimensionless
    InputStream is = null;
    for (@Dimensionless String path : paths) {
      is = tryPath(path + libname);
      if (is != null) {
        logger.debug("Found " + libname + " in " + path);
        libpath = path;
        break;
      }
    }
    return is;
  }

  /**
   * Translate all those Windows to "Windows". ("Windows XP", "Windows Vista", "Windows 7", etc.)
   */
  private static @Dimensionless String unifyOSName(@Dimensionless String osname) {
    if (osname.startsWith("Windows")) {
      return "Windows";
    }
    return osname;
  }

  private static @Dimensionless String getUnifiedOSName() {
    return unifyOSName(System.getProperty("os.name"));
  }

  /**
   * Compute the path to the library. The path is basically
   * "/" + os.name + "/" + os.arch + "/" + libname.
   */
  private @Dimensionless String fatJarLibraryPath(@Dimensionless LibraryLoader this, @Dimensionless String linkage, @Dimensionless String flavor) {
    @Dimensionless
    String sep = "/"; //System.getProperty("file.separator");
    @Dimensionless
    String os_name = getUnifiedOSName();
    @Dimensionless
    String os_arch = System.getProperty("os.arch");
    @Dimensionless
    String path = sep + "lib" + sep + linkage + sep + os_name + sep + os_arch + sep;
    if (null != flavor)
      path += flavor + sep;
    return path;
  }

  /**
   * Try to open a file at the given position.
   */
  private @Dimensionless InputStream tryPath(@Dimensionless LibraryLoader this, @Dimensionless String path) {
    Logger.getLogger().debug("Trying path \"" + path + "\".");
    return getClass().getResourceAsStream(path);
  }

  private @Dimensionless File createTempFile(@Dimensionless LibraryLoader this, @Dimensionless String name) throws IOException {
    return new @Dimensionless File(tempDir + File.separator + name);
  }

  /**
   * Load a system library from a stream. Copies the library to a temp file
   * and loads from there.
   *
   * @param libname name of the library (just used in constructing the library name)
   * @param is      InputStream pointing to the library
   */
  private void loadLibraryFromStream(@Dimensionless LibraryLoader this, @Dimensionless String libname, @Dimensionless InputStream is) {
    try {
      @Dimensionless
      File tempfile = createTempFile(libname);
      @Dimensionless
      OutputStream os = new @Dimensionless FileOutputStream(tempfile);

      logger.debug("tempfile.getPath() = " + tempfile.getPath());

      @ms
      long savedTime = System.currentTimeMillis();

      // Leo says 8k block size is STANDARD ;)
      @Dimensionless
      byte buf @Dimensionless [] = new @Dimensionless byte @Dimensionless [((@Dimensionless int) (8192))];
      @Dimensionless
      int len;
      while ((len = is.read(buf)) > ((@Dimensionless int) (0))) {
        os.write(buf, ((@Dimensionless int) (0)), len);
      }

      os.flush();
      @Dimensionless
      InputStream lock = new @Dimensionless FileInputStream(tempfile);
      os.close();

      @Dimensionless
      double seconds = (@ms double) (System.currentTimeMillis() - savedTime) / ((@ms double) (1e3));
      logger.debug("Copying took " + seconds + " seconds.");

      logger.debug("Loading library from " + tempfile.getPath() + ".");
      System.load(tempfile.getPath());

      lock.close();
    } catch (@Dimensionless IOException io) {
      logger.error("Could not create the temp file: " + io.toString() + ".\n");
    } catch (@Dimensionless UnsatisfiedLinkError ule) {
      logger.error("Couldn't load copied link file: " + ule.toString() + ".\n");
      throw ule;
    }
  }
}
