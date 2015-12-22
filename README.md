Triangle
========

From http://www.cs.cmu.edu/~quake/triangle.html:
> Triangle generates exact Delaunay triangulations, constrained Delaunay triangulations, conforming Delaunay triangulations, Voronoi diagrams, and high-quality triangular meshes. The latter can be generated with no small or large angles, and are thus suitable for finite element analysis.

The original code already allows for building a DLL by using the TRILIBRARY symbol. There is however a problem, since error handling is done by printing a message to console and then calling ```exit(1)```. The main goal of this project is to introduce error codes and return them to the calling code, so using the library from a GUI should be safe.

The project also includes an extension written by Hale Erten and Alper Üngör (see https://www.cise.ufl.edu/~ungor/aCute/index.html).

##Instructions.

The Visual Studio project contains the patched files ready to build.

- Triangle, http://www.cs.cmu.edu/~quake/triangle.html, version 1.6, released 07/28/2005.
- aCute, https://www.cise.ufl.edu/~ungor/aCute/download.html, version 1.0, released 06/15/2009.


**Remarks.**
 - Patch files are no longer included. If you don't want to use the aCute extension, you can make a ```diff``` with the original ```triangle.c``` and remove the aCute specific changes.
 - To compile the project, you can use [Microsoft Visual Studio Community](https://www.visualstudio.com/en-us/products/visual-studio-community-vs.aspx) edition.

##Changes.

triangle.h:
 - adds an integer "errorcode" field to triangulateio struct
 - adds __declspec(dllexport) to exported methods
 
triangle.c:
 - activates some preprocessor definitions (ANSI_DECLARATORS, NO_TIMER, TRILIBRARY, REDUCED)
 - adds error codes to critical functions (so no exit(1) will be called by these functions)
 - removes statistic functions (using REDUCED symbol)
 - adds support for x64 compilation
 - fixes an issue with vertex attributes interpolation
 - applies all changes defined in aCute's "instructions" file
 - adds acute memory pool (acute.h)
 
newSPLocation.h
 - adds acute memory pool (acute.h)
 - corrects some conditionals in doSmoothing
 - removes unused variables, initializes some pointers to NULL
 - removes statistic functions (using TRILIBRARY symbol)

##License.

Please make sure to take a look at the original [README](https://github.com/wo80/Triangle/tree/master/src/Triangle) included in the Triangle source dir.
