
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

void triangleio_reset(triangleio *io)
{
	io->pointlist = (REAL *) NULL;
	io->pointattributelist = (REAL *) NULL;
	io->pointmarkerlist = (int *) NULL;
	io->numberofpoints = 0;
	io->numberofpointattributes = 0;

	io->trianglelist = (int *) NULL;
	io->triangleattributelist = (REAL *) NULL;
	io->trianglearealist = (REAL *) NULL;
	io->neighborlist = (int *) NULL;
	io->numberoftriangles = 0;
	io->numberofcorners = 0;
	io->numberoftriangleattributes = 0;

	io->segmentlist = (int *) NULL;
	io->segmentmarkerlist = (int *) NULL;
	io->numberofsegments = 0;

	io->holelist = (REAL *) NULL;
	io->numberofholes = 0;
	io->regionlist = (REAL *) NULL;
	io->numberofregions = 0;

	io->edgelist = (int *) NULL;
	io->edgemarkerlist = (int *) NULL;
	io->numberofedges = 0;
}

void triangleio_free(triangleio *io)
{
	free(io->pointlist);
	free(io->pointattributelist);
	free(io->pointmarkerlist);

	free(io->trianglelist);
	free(io->triangleattributelist);
	free(io->trianglearealist);
	free(io->neighborlist);

	free(io->segmentlist);
	free(io->segmentmarkerlist);

	free(io->holelist);
	free(io->regionlist);

	free(io->edgelist);
	free(io->edgemarkerlist);
}

void create_rectangle(triangleio *in, REAL left, REAL top, REAL right, REAL bottom)
{
	triangleio_reset(in);

	/* Define input points. */
	in->numberofpoints = 4;
	in->numberofpointattributes = 1;

	in->pointlist = (REAL *) malloc(in->numberofpoints * 2 * sizeof(REAL));
	in->pointlist[0] = left;
	in->pointlist[1] = bottom;
	in->pointlist[2] = right;
	in->pointlist[3] = bottom;
	in->pointlist[4] = right;
	in->pointlist[5] = top;
	in->pointlist[6] = left;
	in->pointlist[7] = top;

	in->pointattributelist = (REAL *) malloc(in->numberofpoints *
		in->numberofpointattributes *
		sizeof(REAL));
	in->pointattributelist[0] = 0.0;
	in->pointattributelist[1] = 1.0;
	in->pointattributelist[2] = 11.0;
	in->pointattributelist[3] = 10.0;

	in->pointmarkerlist = (int *) malloc(in->numberofpoints * sizeof(int));
	in->pointmarkerlist[0] = 0;
	in->pointmarkerlist[1] = 2;
	in->pointmarkerlist[2] = 0;
	in->pointmarkerlist[3] = 0;

	in->numberofsegments = 0;
	in->numberofholes = 0;
	in->numberofregions = 1;

	in->regionlist = (REAL *) malloc(in->numberofregions * 4 * sizeof(REAL));
	in->regionlist[0] = 0.5;
	in->regionlist[1] = 5.0;
	in->regionlist[2] = 7.0;            /* Regional attribute (for whole mesh). */
	in->regionlist[3] = 0.1;          /* Area constraint that will not be used. */
}

void create_rectangle_mesh(triangleio *in, REAL left, REAL top, REAL right, REAL bottom)
{
	triangleio_reset(in);

	create_rectangle(in, left, top, right, bottom);
	
	/* Define segments. */
	in->numberofsegments = 4;
	in->segmentlist = (int *) malloc(in->numberofsegments * 2 * sizeof(int));
	in->segmentlist[0] = 1;
	in->segmentlist[1] = 0;
	in->segmentlist[2] = 2;
	in->segmentlist[3] = 1;
	in->segmentlist[4] = 3;
	in->segmentlist[5] = 2;
	in->segmentlist[6] = 0;
	in->segmentlist[7] = 3;

	in->segmentmarkerlist = (int *) malloc(in->numberofsegments *
		sizeof(int));
	in->segmentmarkerlist[0] = 1;
	in->segmentmarkerlist[1] = 1;
	in->segmentmarkerlist[2] = 1;
	in->segmentmarkerlist[3] = 1;

	/* Define triangles. */
	in->numberofcorners = 3;
	in->numberoftriangles = 2;
	in->numberoftriangleattributes = 1;

	in->trianglelist = (int *) malloc(in->numberoftriangles *
		in->numberofcorners *
		sizeof(int));
	in->trianglelist[0] = 3;
	in->trianglelist[1] = 0;
	in->trianglelist[2] = 1;
	in->trianglelist[3] = 1;
	in->trianglelist[4] = 2;
	in->trianglelist[5] = 3;

	in->triangleattributelist = (REAL *) malloc(in->numberoftriangles *
		sizeof(REAL));
	in->triangleattributelist[0] = 7.0;
	in->triangleattributelist[1] = 7.0;

	/* Triangle area constraints. */
	in->trianglearealist = (REAL *) malloc(in->numberoftriangles * sizeof(REAL));
	in->trianglearealist[0] = 3.0;
	in->trianglearealist[1] = 1.0;
}