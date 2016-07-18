Triangle
========

From http://www.cs.cmu.edu/~quake/triangle.html:
> Triangle generates exact Delaunay triangulations, constrained Delaunay triangulations, conforming Delaunay triangulations, Voronoi diagrams, and high-quality triangular meshes. The latter can be generated with no small or large angles, and are thus suitable for finite element analysis.

The original Triangle code is intended to be compiled as a standalone application. Though it is possible to compile the code as a library by using the `TRILIBRARY` symbol, there are a couple of problems with the approach (for example printing error messages to console and the use of `exit(1)`).

The main goal of this project is to turn Triangle into a re-usable library and the introduction of a simplified C API.

## Contents ##

The Triangle repository contains the following directory structure:

    src/Triangle                 Static library (original Triangle and aCute code)
    src/examples/libtriangle     Dynamic library code (new Triangle C API)
    src/examples/triangle-cli    Command-line interface
    src/examples/triangle-test   Simple test program

The static library code is based on the following sources:

- Triangle (version 1.6), released 07/28/2005  
  Copyright 1993, 1995, 1997, 1998, 2002, 2005 Jonathan Richard Shewchuk  
  http://www.cs.cmu.edu/~quake/triangle.html
- aCute (version 1.0), released 06/15/2009  
  Copyright Hale Erten, Alper Üngör  
  https://www.cise.ufl.edu/~ungor/aCute/download.html

A Visual Studio solution (*Triangle.sln*) can be found in the `src` directory (you can use [Microsoft Visual Studio Community](https://www.visualstudio.com/en-us/products/visual-studio-community-vs.aspx) edition to compile all projects). No platform specific code is used, so all projects should as well compile on Linux or Mac.

If you don't want to use the aCute extension, add `NO_ACUTE` to preprocessor definitions.

##Changes

Changes to Triangle:

 - Remove non-ANSI function declarations (`ANSI_DECLARATORS` symbol no longer used)
 - Remove all non-library code (`TRILIBRARY` symbol no longer used)
 - Remove main `triangulate` function (`NO_TIMER` symbol no longer used)
 - Move structure definitions to `triangle.h` header
 - Move configuration (`#define` constants) to `triangle_config.h` header
 - Create `triangle_internal.h` header containing function prototypes
 - Move robust predicates to separate source file
 - Move file I/O routines to separate source file
 - Remove most ```exit(1)``` calls and return error codes instead
 - Remove unused members from `mesh` and `behavior` structs
 - Add *experimental* support for x64 compilation
 - Include aCute for quality mesh generation
 
Changes to aCute:

 - Introduction of memory pool
 - Cleanup and minor fixes

Please refer to the commit history if you need a complete changelog.

##License

Please note that although both Triangle and aCute are freely available to researchers, they may not be sold or included in commercial products without a license. Make sure to take a look at the original [README](https://github.com/wo80/Triangle/tree/master/src/Triangle) included in the Triangle source dir.
