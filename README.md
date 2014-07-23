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
 - If you don't want to use the aCute extension, you can use *triangle.c.patch-x* file, which will apply all patches that are not aCute specific..
 - To compile the project, you can use [Microsoft Visual Studio Express 2013 for Windows Desktop](http://www.visualstudio.com/downloads/download-visual-studio-vs#d-express-windows-desktop).

##Patches.

triangle.h.patch:
 - adds an integer "errorcode" field to triangulateio struct
 - adds __declspec(dllexport) to exported methods
 
triangle.c.patch:
 - applies all changes defined in aCute's "instructions" file
 - activates some preprocessor definitions (ANSI_DECLARATORS, NO_TIMER, TRILIBRARY, REDUCED)
 - adds error codes to critical functions (so no exit(1) will be called by these functions)
 - removes statistic functions (using REDUCED symbol)
 - adds acute memory pool (acute.h)
 
newSPLocation.h
 - adds acute memory pool (acute.h)
 - corrects some conditionals in doSmoothing
 - removes unused variables, initializes some pointers to NULL
 - removes statistic functions (using TRILIBRARY symbol)

##License.

The patch files are released to the public domain without any special license. Note, however, that the original code and produced binaries will stay under the license/copyright the original authors intended. Particularly, make sure to have a look at the README included in the Triangle source dir.
