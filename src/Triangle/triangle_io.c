
#include "triangle.h"
#include "triangle_internal.h"

// Defined in triangle.c
extern int plus1mod3[3];
extern int minus1mod3[3];

/********* File writing routines begin here                          *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  writenodes()   Number the vertices and write them to a .node file.       */
/*                                                                           */
/*  To save memory, the vertex numbers are written over the boundary markers */
/*  after the vertices are written to a file.                                */
/*                                                                           */
/*****************************************************************************/

int file_writenodes(mesh *m, behavior *b, FILE *nodefile)
{
	vertex vertexloop;
	long outvertices;
	int vertexnumber;
	int i;

	if (b->jettison) {
		outvertices = m->vertices.items - m->undeads;
	} else {
		outvertices = m->vertices.items;
	}

	if (nodefile == (FILE *) NULL) {
		return -1;
	}
	/* Number of vertices, number of dimensions, number of vertex attributes, */
	/*   and number of boundary markers (zero or one).                        */
	fprintf(nodefile, "%ld  %d  %d  %d\n", outvertices, m->mesh_dim,
		m->nextras, 1 - b->nobound);

	traversalinit(&m->vertices);
	vertexnumber = b->firstnumber;
	vertexloop = vertextraverse(m);
	while (vertexloop != (vertex) NULL) {
		if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) {
			/* Vertex number, x and y coordinates. */
			fprintf(nodefile, "%4d    %.17g  %.17g", vertexnumber, vertexloop[0],
				vertexloop[1]);
			for (i = 0; i < m->nextras; i++) {
				/* Write an attribute. */
				fprintf(nodefile, "  %.17g", vertexloop[i + 2]);
			}
			if (b->nobound) {
				fprintf(nodefile, "\n");
			} else {
				/* Write the boundary marker. */
				fprintf(nodefile, "    %d\n", vertexmark(vertexloop));
			}

			setvertexmark(vertexloop, vertexnumber);
			vertexnumber++;
		}
		vertexloop = vertextraverse(m);
	}
	return 0;
}

/*****************************************************************************/
/*                                                                           */
/*  writeelements()   Write the triangles to an .ele file.                   */
/*                                                                           */
/*****************************************************************************/

int file_writeelements(mesh *m, behavior *b, FILE *elefile)
{
	struct otri triangleloop;
	vertex p1, p2, p3;
	vertex mid1, mid2, mid3;
	long elementnumber;
	int i;

	if (elefile == (FILE *) NULL) {
		return -1;
	}
	/* Number of triangles, vertices per triangle, attributes per triangle. */
	fprintf(elefile, "%ld  %d  %d\n", m->triangles.items,
		(b->order + 1) * (b->order + 2) / 2, m->eextras);

	traversalinit(&m->triangles);
	triangleloop.tri = triangletraverse(m);
	triangleloop.orient = 0;
	elementnumber = b->firstnumber;
	while (triangleloop.tri != (triangle *) NULL) {
		org(triangleloop, p1);
		dest(triangleloop, p2);
		apex(triangleloop, p3);
		if (b->order == 1) {
			/* Triangle number, indices for three vertices. */
			fprintf(elefile, "%4ld    %4d  %4d  %4d", elementnumber,
				vertexmark(p1), vertexmark(p2), vertexmark(p3));
		} else {
			mid1 = (vertex) triangleloop.tri[m->highorderindex + 1];
			mid2 = (vertex) triangleloop.tri[m->highorderindex + 2];
			mid3 = (vertex) triangleloop.tri[m->highorderindex];
			/* Triangle number, indices for six vertices. */
			fprintf(elefile, "%4ld    %4d  %4d  %4d  %4d  %4d  %4d", elementnumber,
				vertexmark(p1), vertexmark(p2), vertexmark(p3), vertexmark(mid1),
				vertexmark(mid2), vertexmark(mid3));
		}

		for (i = 0; i < m->eextras; i++) {
			fprintf(elefile, "  %.17g", elemattribute(triangleloop, i));
		}
		fprintf(elefile, "\n");

		triangleloop.tri = triangletraverse(m);
		elementnumber++;
	}
	return 0;
}

/*****************************************************************************/
/*                                                                           */
/*  writepoly()   Write the segments and holes to a .poly file.              */
/*                                                                           */
/*****************************************************************************/

int file_writepoly(mesh *m, behavior *b, FILE *polyfile,
				   REAL *holelist, int holes, REAL *regionlist, int regions)
{
	long holenumber, regionnumber;
	struct osub subsegloop;
	vertex endpoint1, endpoint2;
	long subsegnumber;

	if (polyfile == (FILE *) NULL) {
		return -1;
	}
	/* The zero indicates that the vertices are in a separate .node file. */
	/*   Followed by number of dimensions, number of vertex attributes,   */
	/*   and number of boundary markers (zero or one).                    */
	fprintf(polyfile, "%d  %d  %d  %d\n", 0, m->mesh_dim, m->nextras,
		1 - b->nobound);
	/* Number of segments, number of boundary markers (zero or one). */
	fprintf(polyfile, "%ld  %d\n", m->subsegs.items, 1 - b->nobound);

	traversalinit(&m->subsegs);
	subsegloop.ss = subsegtraverse(m);
	subsegloop.ssorient = 0;
	subsegnumber = b->firstnumber;
	while (subsegloop.ss != (subseg *) NULL) {
		sorg(subsegloop, endpoint1);
		sdest(subsegloop, endpoint2);
		/* Segment number, indices of its two endpoints, and possibly a marker. */
		if (b->nobound) {
			fprintf(polyfile, "%4ld    %4d  %4d\n", subsegnumber,
				vertexmark(endpoint1), vertexmark(endpoint2));
		} else {
			fprintf(polyfile, "%4ld    %4d  %4d    %4d\n", subsegnumber,
				vertexmark(endpoint1), vertexmark(endpoint2), mark(subsegloop));
		}

		subsegloop.ss = subsegtraverse(m);
		subsegnumber++;
	}

#ifndef CDT_ONLY
	fprintf(polyfile, "%d\n", holes);
	if (holes > 0) {
		for (holenumber = 0; holenumber < holes; holenumber++) {
			/* Hole number, x and y coordinates. */
			fprintf(polyfile, "%4ld   %.17g  %.17g\n", b->firstnumber + holenumber,
				holelist[2 * holenumber], holelist[2 * holenumber + 1]);
		}
	}
	if (regions > 0) {
		fprintf(polyfile, "%d\n", regions);
		for (regionnumber = 0; regionnumber < regions; regionnumber++) {
			/* Region number, x and y coordinates, attribute, maximum area. */
			fprintf(polyfile, "%4ld   %.17g  %.17g  %.17g  %.17g\n",
				b->firstnumber + regionnumber,
				regionlist[4 * regionnumber], regionlist[4 * regionnumber + 1],
				regionlist[4 * regionnumber + 2],
				regionlist[4 * regionnumber + 3]);
		}
	}
#endif /* not CDT_ONLY */
	return 0;
}

/*****************************************************************************/
/*                                                                           */
/*  writeedges()   Write the edges to an .edge file.                         */
/*                                                                           */
/*****************************************************************************/

int file_writeedges(mesh *m, behavior *b, FILE *edgefile)
{
	struct otri triangleloop, trisym;
	struct osub checkmark;
	vertex p1, p2;
	long edgenumber;
	triangle ptr;                         /* Temporary variable used by sym(). */
	subseg sptr;                      /* Temporary variable used by tspivot(). */

	if (edgefile == (FILE *) NULL) {
		return -1;
	}
	/* Number of edges, number of boundary markers (zero or one). */
	fprintf(edgefile, "%ld  %d\n", m->edges, 1 - b->nobound);

	traversalinit(&m->triangles);
	triangleloop.tri = triangletraverse(m);
	edgenumber = b->firstnumber;
	/* To loop over the set of edges, loop over all triangles, and look at   */
	/*   the three edges of each triangle.  If there isn't another triangle  */
	/*   adjacent to the edge, operate on the edge.  If there is another     */
	/*   adjacent triangle, operate on the edge only if the current triangle */
	/*   has a smaller pointer than its neighbor.  This way, each edge is    */
	/*   considered only once.                                               */
	while (triangleloop.tri != (triangle *) NULL) {
		for (triangleloop.orient = 0; triangleloop.orient < 3;
			triangleloop.orient++) {
				sym(triangleloop, trisym);
				if ((triangleloop.tri < trisym.tri) || (trisym.tri == m->dummytri)) {
					org(triangleloop, p1);
					dest(triangleloop, p2);
					if (b->nobound) {
						/* Edge number, indices of two endpoints. */
						fprintf(edgefile, "%4ld   %d  %d\n", edgenumber,
							vertexmark(p1), vertexmark(p2));
					} else {
						/* Edge number, indices of two endpoints, and a boundary marker. */
						/*   If there's no subsegment, the boundary marker is zero.      */
						if (b->usesegments) {
							tspivot(triangleloop, checkmark);
							if (checkmark.ss == m->dummysub) {
								fprintf(edgefile, "%4ld   %d  %d  %d\n", edgenumber,
									vertexmark(p1), vertexmark(p2), 0);
							} else {
								fprintf(edgefile, "%4ld   %d  %d  %d\n", edgenumber,
									vertexmark(p1), vertexmark(p2), mark(checkmark));
							}
						} else {
							fprintf(edgefile, "%4ld   %d  %d  %d\n", edgenumber,
								vertexmark(p1), vertexmark(p2), trisym.tri == m->dummytri);
						}
					}
					edgenumber++;
				}
		}
		triangleloop.tri = triangletraverse(m);
	}
	return 0;
}

int file_writeneighbors(mesh *m, behavior *b, FILE *neighborfile)
{
	struct otri triangleloop, trisym;
	long elementnumber;
	int neighbor1, neighbor2, neighbor3;
	triangle ptr;                         /* Temporary variable used by sym(). */

	if (neighborfile == (FILE *) NULL) {
		return -1;
	}
	/* Number of triangles, three neighbors per triangle. */
	fprintf(neighborfile, "%ld  %d\n", m->triangles.items, 3);

	traversalinit(&m->triangles);
	triangleloop.tri = triangletraverse(m);
	triangleloop.orient = 0;
	elementnumber = b->firstnumber;
	while (triangleloop.tri != (triangle *) NULL) {
		* (int *) (triangleloop.tri + 6) = (int) elementnumber;
		triangleloop.tri = triangletraverse(m);
		elementnumber++;
	}
	* (int *) (m->dummytri + 6) = -1;

	traversalinit(&m->triangles);
	triangleloop.tri = triangletraverse(m);
	elementnumber = b->firstnumber;
	while (triangleloop.tri != (triangle *) NULL) {
		triangleloop.orient = 1;
		sym(triangleloop, trisym);
		neighbor1 = * (int *) (trisym.tri + 6);
		triangleloop.orient = 2;
		sym(triangleloop, trisym);
		neighbor2 = * (int *) (trisym.tri + 6);
		triangleloop.orient = 0;
		sym(triangleloop, trisym);
		neighbor3 = * (int *) (trisym.tri + 6);
		/* Triangle number, neighboring triangle numbers. */
		fprintf(neighborfile, "%4ld    %d  %d  %d\n", elementnumber,
			neighbor1, neighbor2, neighbor3);

		triangleloop.tri = triangletraverse(m);
		elementnumber++;
	}
	return 0;
}

/**                                                                         **/
/**                                                                         **/
/********* File writing routines end here                            *********/


/********* File reading routines begin here                          *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  readline()   Read a nonempty line from a file.                           */
/*                                                                           */
/*  A line is considered "nonempty" if it contains something that looks like */
/*  a number.  Comments (prefaced by `#') are ignored.                       */
/*                                                                           */
/*****************************************************************************/

char *readline(char *string, FILE *infile)
{
	char *result;

	/* Search for something that looks like a number. */
	do {
		result = fgets(string, INPUTLINESIZE, infile);
		if (result == (char *) NULL) {
			return result; /* NULL must be checked by caller! */
		}
		/* Skip anything that doesn't look like a number, a comment, */
		/*   or the end of a line.                                   */
		while ((*result != '\0') && (*result != '#')
			&& (*result != '.') && (*result != '+') && (*result != '-')
			&& ((*result < '0') || (*result > '9'))) {
				result++;
		}
		/* If it's a comment or end of line, read another line and try again. */
	} while ((*result == '#') || (*result == '\0'));
	return result;
}

/*****************************************************************************/
/*                                                                           */
/*  findfield()   Find the next field of a string.                           */
/*                                                                           */
/*  Jumps past the current field by searching for whitespace, then jumps     */
/*  past the whitespace to find the next field.                              */
/*                                                                           */
/*****************************************************************************/

char *findfield(char *string)
{
	char *result;

	result = string;
	/* Skip the current field.  Stop upon reaching whitespace. */
	while ((*result != '\0') && (*result != '#')
		&& (*result != ' ') && (*result != '\t')) {
			result++;
	}
	/* Now skip the whitespace and anything else that doesn't look like a */
	/*   number, a comment, or the end of a line.                         */
	while ((*result != '\0') && (*result != '#')
		&& (*result != '.') && (*result != '+') && (*result != '-')
		&& ((*result < '0') || (*result > '9'))) {
			result++;
	}
	/* Check for a comment (prefixed with `#'). */
	if (*result == '#') {
		*result = '\0';
	}
	return result;
}

/*****************************************************************************/
/*                                                                           */
/*  file_readnodes_internal()   Read the vertices from a file, which may be  */
/*                              a .node or .poly file.                       */
/*                                                                           */
/*****************************************************************************/

int file_readnodes_internal(FILE *file, triangleio *io, int* numvertices, int* firstnode)
{
	char inputline[INPUTLINESIZE];
	char *stringptr;
	REAL x, y;
	int nodemarkers;
	int invertices;
	int mesh_dim;
	int nextras;
	int i, j;

	if (file == (FILE *) NULL) {
		return -1; // TODO: read error: cannot access file.
	}

	/* Read number of vertices, number of dimensions, number of vertex */
	/*   attributes, and number of boundary markers.                   */
	stringptr = readline(inputline, file);
	if (stringptr != (char *) NULL) {
		return -1; // TODO: read error
	}
	invertices = (int) strtol(stringptr, &stringptr, 0);
	stringptr = findfield(stringptr);
	if (*stringptr == '\0') {
		mesh_dim = 2;
	} else {
		mesh_dim = (int) strtol(stringptr, &stringptr, 0);
	}
	stringptr = findfield(stringptr);
	if (*stringptr == '\0') {
		nextras = 0;
	} else {
		nextras = (int) strtol(stringptr, &stringptr, 0);
	}
	stringptr = findfield(stringptr);
	if (*stringptr == '\0') {
		nodemarkers = 0;
	} else {
		nodemarkers = (int) strtol(stringptr, &stringptr, 0);
	}

	if (invertices < 3) {
		return -1; // TODO: error: Input must have at least three input vertices.
	}
	if (mesh_dim != 2) {
		return -1; // TODO: error: Triangle only works with two-dimensional meshes.
	}
	if (nextras == 0) {
		//b->weighted = 0;
	}

	io->numberofpoints = invertices;
	io->numberofpointattributes = nextras;
	io->pointlist = (REAL *)trimalloc(2 * invertices * sizeof(REAL));
	io->pointattributelist = (REAL *)trimalloc(nextras * invertices * sizeof(REAL));
	if (nodemarkers) {
		io->pointmarkerlist = (int *)trimalloc(invertices * sizeof(int));
	}

	/* Read the vertices. */
	for (i = 0; i < invertices; i++) {
		stringptr = readline(inputline, file);
		if (stringptr != (char *) NULL) {
			return -1; // TODO: read error
		}
		if (i == 0) {
			*firstnode = (int) strtol(stringptr, &stringptr, 0);
		}
		stringptr = findfield(stringptr);
		if (*stringptr == '\0') {
			return -1; // TODO: error: Vertex (b->firstnumber + i) has no x coordinate.
		}
		x = (REAL) strtod(stringptr, &stringptr);
		stringptr = findfield(stringptr);
		if (*stringptr == '\0') {
			return -1; // TODO: error: Vertex (b->firstnumber + i) has no y coordinate.
		}
		y = (REAL) strtod(stringptr, &stringptr);
		io->pointlist[2 * i + 0] = x;
		io->pointlist[2 * i + 1] = y;
		/* Read the vertex attributes. */
		for (j = 0; j < nextras; j++) {
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				io->pointattributelist[nextras * i + j] = 0.0;
			} else {
				io->pointattributelist[nextras * i + j] = (REAL) strtod(stringptr, &stringptr);
			}
		}
		if (nodemarkers) {
			/* Read a vertex marker. */
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				io->pointmarkerlist[i] = 0;
			} else {
				io->pointmarkerlist[i] = (int) strtol(stringptr, &stringptr, 0);
			}
		}
	}

	return 0;
}

int file_readsegments_internal(FILE *file, triangleio *io, int invertices)
{
	char inputline[INPUTLINESIZE];
	char *stringptr;

	int segmentmarkers;
	int boundmarker;
	int insegments;
	int i;

	/* Read the segments from a .poly file. */
	/* Read number of segments and number of boundary markers. */
	stringptr = readline(inputline, file);
	if (stringptr != (char *) NULL) {
		return -1; // TODO: read error
	}
	insegments = (int) strtol(stringptr, &stringptr, 0);
	stringptr = findfield(stringptr);
	if (*stringptr == '\0') {
		segmentmarkers = 0;
	} else {
		segmentmarkers = (int) strtol(stringptr, &stringptr, 0);
	}

	if (invertices == 0) {
		return - 1;
	}

	io->numberofsegments = insegments;
	io->segmentlist = (int *)trimalloc(2 * insegments * sizeof(int));
	if (segmentmarkers) {
		io->segmentmarkerlist = (int *)trimalloc(insegments * sizeof(int));
	}

	boundmarker = 0;
	/* Read the segments. */
	for (i = 0; i < insegments; i++) {
		stringptr = readline(inputline, file);
		if (stringptr != (char *) NULL) {
			return -1; // TODO: read error
		}
		stringptr = findfield(stringptr);
		if (*stringptr == '\0') {
			return -1; // TODO: error: Segment (firstnumber + i) has no endpoints.
		} else {
			io->segmentlist[2 * i] = (int) strtol(stringptr, &stringptr, 0);
		}
		stringptr = findfield(stringptr);
		if (*stringptr == '\0') {
			return -1; // TODO: error: Segment (firstnumber + i) is missing its second endpoint.
		} else {
			io->segmentlist[2 * i + 1] = (int) strtol(stringptr, &stringptr, 0);
		}
		if (segmentmarkers) {
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				io->segmentmarkerlist[i] = 0;
			} else {
				io->segmentmarkerlist[i] = (int) strtol(stringptr, &stringptr, 0);
			}
		}
	}
	return 0;
}

/*****************************************************************************/
/*                                                                           */
/*  file_readholes_internal()   Read the holes, and possibly regional attributes  */
/*                              and area  constraints, from a .poly file.    */
/*                                                                           */
/*****************************************************************************/

int file_readholes_internal(FILE *file, triangleio *io)
{
	char inputline[INPUTLINESIZE];
	char *stringptr;
	int index;
	int holes;
	int regions;
	int i;

	/* Read the holes. */
	stringptr = readline(inputline, file);
	if (stringptr != (char *) NULL) {
		return -1; // TODO: read error
	}
	holes = (int) strtol(stringptr, &stringptr, 0);
	if (holes > 0) {
		io->numberofholes = holes;
		io->holelist = (REAL *) trimalloc(2 * holes * (int) sizeof(REAL));
		for (i = 0; i < 2 * holes; i += 2) {
			stringptr = readline(inputline, file);
			if (stringptr != (char *) NULL) {
				return -1; // TODO: read error
			}
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				return -1; // TODO: error: Hole (b->firstnumber + (i >> 1)) has no x coordinate.
			} else {
				io->holelist[i] = (REAL) strtod(stringptr, &stringptr);
			}
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				return -1; // TODO: error: Hole (b->firstnumber + (i >> 1)) has no y coordinate.
			} else {
				io->holelist[i + 1] = (REAL) strtod(stringptr, &stringptr);
			}
		}
	}

#ifndef CDT_ONLY
	/* Read the area constraints. */
	stringptr = readline(inputline, file);
	if (stringptr != (char *) NULL) {
		return -1; // TODO: read error
	}
	regions = (int) strtol(stringptr, &stringptr, 0);
	if (regions > 0) {
		io->numberofregions = regions;
		io->regionlist = (REAL *) trimalloc(4 * regions * (int) sizeof(REAL));
		index = 0;
		for (i = 0; i < regions; i++) {
			stringptr = readline(inputline, file);
			if (stringptr != (char *) NULL) {
				return -1; // TODO: read error
			}
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				return -1; // TODO: error: Region (b->firstnumber + i) has no x coordinate.
			} else {
				io->regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
			}
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				return -1; // TODO: error: Region (b->firstnumber + i) has no y coordinate.
			} else {
				io->regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
			}
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				return -1; // TODO: error: Region (b->firstnumber + i) has no region attribute or area constraint.
			} else {
				io->regionlist[index++] = (REAL) strtod(stringptr, &stringptr);
			}
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				io->regionlist[index] = io->regionlist[index - 1];
			} else {
				io->regionlist[index] = (REAL) strtod(stringptr, &stringptr);
			}
			index++;
		}
	}
#endif /* not CDT_ONLY */

	return 0;
}

/*****************************************************************************/
/*                                                                           */
/*  readnodes()   Read the vertices from a .node file.                       */
/*                                                                           */
/*****************************************************************************/

int file_readnodes(FILE *nodefile, triangleio *io, int *firstnode)
{
	int numvertices;

	return file_readnodes_internal(nodefile, io, &numvertices, firstnode);
}

/*****************************************************************************/
/*                                                                           */
/*  file_readpoly()   Read vertices and segments from a .poly file.          */
/*                                                                           */
/*****************************************************************************/

int file_readpoly(FILE *polyfile, triangleio *io, int *firstnode)
{
	int numvertices;

	file_readnodes_internal(polyfile, io, &numvertices, firstnode);
	file_readsegments_internal(polyfile, io, numvertices);
	file_readholes_internal(polyfile, io);

	return 0;
}

/*****************************************************************************/
/*                                                                           */
/*  file_readelements()  Read mesh elements from an .ele file.               */
/*                                                                           */
/*****************************************************************************/

int file_readelements(FILE *file, triangleio *io)
{
	char inputline[INPUTLINESIZE];
	char *stringptr;
	int inelements;
	int eextras;

	int incorners;
	int i, j;

	/* Read the triangles from an .ele file. */
	if (file == (FILE *) NULL) {
		return -1;
	}
	/* Read number of triangles, number of vertices per triangle, and */
	/*   number of triangle attributes from .ele file.                */
	stringptr = readline(inputline, file);
	if (stringptr != (char *) NULL) {
		return -1; // TODO: read error
	}
	inelements = (int) strtol(stringptr, &stringptr, 0);
	stringptr = findfield(stringptr);
	if (*stringptr == '\0') {
		incorners = 3;
	} else {
		incorners = (int) strtol(stringptr, &stringptr, 0);
		if (incorners < 3) {
			return -1; // TODO: error: Triangles must have at least 3 vertices.
		}
	}
	stringptr = findfield(stringptr);
	if (*stringptr == '\0') {
		eextras = 0;
	} else {
		eextras = (int) strtol(stringptr, &stringptr, 0);
	}

	io->numberoftriangles = inelements;
	io->trianglelist = (int *)trimalloc(3 * inelements * sizeof(int));
	if (eextras > 0) {
		io->numberoftriangleattributes = eextras;
		io->triangleattributelist = (REAL *)trimalloc(eextras * inelements * sizeof(REAL));
	}

	/* Read the triangles from the .ele file. */
	for (i = 0; i < inelements; i++) {
		/* Read triangle number and the triangle's three corners. */
		stringptr = readline(inputline, file);
		if (stringptr != (char *) NULL) {
			return -1; // TODO: read error
		}
		for (j = 0; j < 3; j++) {
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				return -1; // TODO: error: Triangle (elementnumber) is missing vertex (j + 1).
			} else {
				io->trianglelist[3 * i + j] = (int) strtol(stringptr, &stringptr, 0);
			}
		}

		/* Read the triangle's attributes. */
		for (j = 0; j < eextras; j++) {
			stringptr = findfield(stringptr);
			if (*stringptr == '\0') {
				io->triangleattributelist[eextras * i + j] = 0;
			} else {
				io->triangleattributelist[eextras * i + j] = (REAL) strtod(stringptr, &stringptr);
			}
		}
	}

	return 0;
}

/*****************************************************************************/
/*                                                                           */
/*  file_readelementsarea()  Read elements area from an .area file.          */
/*                                                                           */
/*****************************************************************************/

int file_readelementsarea(FILE *file, triangleio *io, int numelements)
{
	char inputline[INPUTLINESIZE];
	char *stringptr;
	int areaelements;
	int i;

	/* Read the triangles area from an .area file. */
	if (file == (FILE *) NULL) {
		return -1;
	}

	/* Check for consistency with the .ele file. */
	stringptr = readline(inputline, file);
	if (stringptr != (char *) NULL) {
		return -1; // TODO: read error
	}
	areaelements = (int) strtol(stringptr, &stringptr, 0);
	if (areaelements != numelements) {
		return -1; // TODO: error: area file disagrees on number of triangles.
	}

	io->trianglearealist = (REAL *)trimalloc(numelements * sizeof(REAL));

	/* Read the triangles from the .ele file. */
	for (i = 0; i < numelements; i++) {
		/* Read an area constraint from the .area file. */
		stringptr = readline(inputline, file);
		if (stringptr != (char *) NULL) {
			return -1; // TODO: read error
		}
		stringptr = findfield(stringptr);
		if (*stringptr == '\0') {
			io->trianglearealist[i] = -1.0; /* No constraint on this triangle. */
		} else {
			io->trianglearealist[i] = (REAL) strtod(stringptr, &stringptr);
		}
	}

	return 0;
}

/**                                                                         **/
/**                                                                         **/
/********* File reading routines end here                            *********/