Triangle
========

Visual Studio project and patch files for Jonathan Shewchuk's Triangle mesh generator.

##Instructions.

1. Download the current Triangle release form http://www.cs.cmu.edu/~quake/triangle.html (version 1.6, released 07/28/2005)
2. Replace the triangle.h and triangle.c files with the ones from the Triangle archive.
3. Apply patches triangle.h.patch and triangle.c.patch.
4. Download the current aCute release form https://www.cise.ufl.edu/~ungor/aCute/download.html (version 1.0, released 06/15/2009)
5. Replace newSPLocation.h with the one from the aCute archive.
6. Apply patch newSPLocation.h.patch.
7. Build the DLL using the Visual Studio project.

**Remark** If you don't want to use the aCute extension, use the triangle.c.patch-x file in step 3 and ignore steps 4-6.

##Patches.

triangle.h.patch:
 - adds the "errorcode" field to triangulateio struct
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
