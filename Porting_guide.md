# Porting guide from legacy triangle

If you're coming from the original triangle, the new triangle_api requires some changes.
This is a first attempt to give a rough porting guide.
I'm not the author of triangle_api and don't know its in's and out's, so YMMV.

## Code snippets

### Data types
Instead of `struct triangulateio` there is now a type `triangleio`.
Additionally, the global state has been deprecated in favor of `context` and `behavior` objects.

````C
// legacy code
struct triangulateio in, out;
````
becomes:
````C
context *ctx;
triangleio in, out;
````

### Initialization and cleanup

You still need to initialize (and free) the `triangleio` data structure just like `struct triangulateio` before.

Additionally, you need:
````C
// init
ctx = triangle_context_create();
````

````C
// cleanup
triangle_context_destroy(ctx);
````

### Triangulation

This part has seen drastic changes: the `triangulate` function has been removed, which means you have to call several functions instead of one.

For this example, I leave out some details:
 - most file io was omitted.
   In particular, the following parts of the `triangleio` struct `in` have been filled already:
   - pointlist
   - pointmarkerlist
   - segmentlist
   - segmentmarkerlist
   - holelist
   - regionlist
 - no mesh refinement is done

A more detailed example can be found by studying the code of `examples/triangle-cli/main.c`.

````C
// legacy code
char cmdline[] = "...";
int tristatus;

tristatus = triangulate(cmdline, &in, &out, (struct triangulateio *) nullptr, TriMessageFunction);
if (tristatus != 0) handle_error();

// process results
...
````

becomes:
````C
char cmdline[] = "...";
int tristatus;

tristatus = triangle_context_options(ctx, cmdline);
if (tristatus != TRI_OK) handle_error();

// Triangulate the polygon
tristatus = triangle_mesh_create(ctx, &in);
if (tristatus != TRI_OK) handle_error();


// process results
...

````

----
Copyright (c) 2017 Johannes Zarl-Zierl.

SPDX-License-Identifier: GFDL-1.3
License-Filename: LICENSES/fdl-1.3.txt

Permission is granted to copy, distribute and/or modify this document under the terms of the GNU Free Documentation License, Version 1.3 or any later version published by the Free Software Foundation; with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.