Triangle
========

From http://www.cs.cmu.edu/~quake/triangle.html:
> Triangle generates exact Delaunay triangulations, constrained Delaunay triangulations, conforming Delaunay triangulations, Voronoi diagrams, and high-quality triangular meshes. The latter can be generated with no small or large angles, and are thus suitable for finite element analysis.

The original code even allows for building DLLs by using the TRILIBRARY symbol. There is however a problem, since error handling is done by printing a message to console and then calling ```exit(1)```. The main goal of this project is to introduce error codes and return them to the calling code, so using the library from a GUI should be safe.

##Instructions.

1. Download the current Triangle release form http://www.cs.cmu.edu/~quake/triangle.html (version 1.6, released 07/28/2005)
2. Replace the *triangle.h* and *triangle.c* files with the ones from the Triangle archive.
3. Apply patches *triangle.h.patch* and *triangle.c.patch*.
4. Download the current aCute release form https://www.cise.ufl.edu/~ungor/aCute/download.html (version 1.0, released 06/15/2009)
5. Replace *newSPLocation.h* with the one from the aCute archive.
6. Apply patch *newSPLocation.h.patch*.
7. Build the DLL using the Visual Studio project.

**Remarks.**
 - If you don't want to use the aCute extension, use the *triangle.c.patch-x* file in step 3 and ignore steps 4-6.
 - If you need a tool for applying the patch files, try http://wo80.bplaced.net/projects/patch/
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

The patch files are released to the public domain without any special license. Note, however, that the original code and produced binaries will stay under the license/copyright the orignal authors intended. Particularly, make sure to have a look at the README included in the Triangle archive.
