
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

#include "triangle_config.h"

/*****************************************************************************/
/*                                                                           */
/*  syntax()   Print list of command line switches.                          */
/*                                                                           */
/*****************************************************************************/

void syntax()
{
#ifdef CDT_ONLY
#ifdef REDUCED
	printf("triangle [-pAcjenBPNEIOXzo_lQVh] input_file\n");
#else /* not REDUCED */
	printf("triangle [-pAcjenBPNEIOXzo_iFlCQVh] input_file\n");
#endif /* not REDUCED */
#else /* not CDT_ONLY */
#ifdef REDUCED
	printf("triangle [-prq__a__uAcDjenBPNEIOXzo_YS__lQVh] input_file\n");
#else /* not REDUCED */
	printf("triangle [-prq__a__uAcDjenBPNEIOXzo_YS__iFlsCQVh] input_file\n");
#endif /* not REDUCED */
#endif /* not CDT_ONLY */

	printf("    -p  Triangulates a Planar Straight Line Graph (.poly file).\n");
#ifndef CDT_ONLY
	printf("    -r  Refines a previously generated mesh.\n");
	printf("    -q  Quality mesh generation.  A minimum angle may be specified.\n");
	printf("    -a  Applies a maximum triangle area constraint.\n");
	printf("    -u  Applies a user-defined triangle constraint.\n");
#endif /* not CDT_ONLY */
	printf("    -A  Applies attributes to identify triangles in certain regions.\n");
	printf("    -c  Encloses the convex hull with segments.\n");
#ifndef CDT_ONLY
	printf("    -D  Conforming Delaunay:  all triangles are truly Delaunay.\n");
#endif /* not CDT_ONLY */
	printf("    -w  Weighted Delaunay triangulation.\n");
	printf("    -W  Regular triangulation (lower hull of a height field).\n");
	printf("    -j  Jettison unused vertices from output .node file.\n");
	printf("    -e  Generates an edge list.\n");
	printf("    -n  Generates a list of triangle neighbors.\n");
	printf("    -B  Suppresses output of boundary information.\n");
	printf("    -P  Suppresses output of .poly file.\n");
	printf("    -N  Suppresses output of .node file.\n");
	printf("    -E  Suppresses output of .ele file.\n");
	printf("    -I  Suppresses mesh iteration numbers.\n");
	printf("    -O  Ignores holes in .poly file.\n");
	printf("    -X  Suppresses use of exact arithmetic.\n");
	printf("    -z  Numbers all items starting from zero (rather than one).\n");
	printf("    -o2 Generates second-order subparametric elements.\n");
#ifndef CDT_ONLY
	printf("    -Y  Suppresses boundary segment splitting.\n");
	printf("    -S  Specifies maximum number of added Steiner points.\n");
#endif /* not CDT_ONLY */
#ifndef REDUCED
	printf("    -i  Uses incremental method, rather than divide-and-conquer.\n");
	printf("    -F  Uses Fortune's sweepline algorithm, rather than d-and-c.\n");
#endif /* not REDUCED */
	printf("    -l  Uses vertical cuts only, rather than alternating cuts.\n");
#ifndef REDUCED
#ifndef CDT_ONLY
	printf("    -s  Force segments into mesh by splitting (instead of using CDT).\n");
#endif /* not CDT_ONLY */
	printf("    -C  Check consistency of final mesh.\n");
#endif /* not REDUCED */
	printf("    -Q  Quiet:  No terminal output except errors.\n");
	printf("    -V  Verbose:  Detailed information on what I'm doing.\n");
	printf("    -h  Help:  Detailed instructions for Triangle.\n");
	exit(0);
}

/*****************************************************************************/
/*                                                                           */
/*  info()   Print out complete instructions.                                */
/*                                                                           */
/*****************************************************************************/

void info()
{
	printf("Triangle (version 1.6)\n");
	printf("A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.\n");
	printf("\n");
	printf("Copyright 1993, 1995, 1997, 1998, 2002, 2005 Jonathan Richard Shewchuk\n");
	printf("\n");
	printf("Detailed instructions can be found in the README file or online at\n");
	printf("http://www.cs.cmu.edu/~quake/triangle.html.\n");
	/*
#ifdef SINGLE
	printf("\nThis executable is compiled for single precision arithmetic.");
#else
	printf("\nThis executable is compiled for double precision arithmetic.");
#endif
	*/
	exit(0);
}

/*****************************************************************************/
/*                                                                           */
/*  parsecommandline()   Read the command line, identify switches, and set   */
/*                       up options and file names.                          */
/*                                                                           */
/*****************************************************************************/

void parsecommandline_legacy(int argc, char **argv, behavior_legacy *b)
{
#define STARTINDEX 1
	int increment;
	int meshnumber;
	int i, j;
	char workstring[FILENAMESIZE];

	b->firstnumber = 1;
	b->edgesout = b->epsout = b->neighbors = 0;
	b->nobound = b->nopolywritten = b->nonodewritten = b->noelewritten = 0;
	b->noiterationnum = 0;
	b->noholes = 0;
	b->docheck = 0;
	b->quiet = b->verbose = 0;
	b->innodefilename[0] = '\0';

	for (i = STARTINDEX; i < argc; i++) {
		if (argv[i][0] == '-') {
			for (j = STARTINDEX; argv[i][j] != '\0'; j++) {
				if (argv[i][j] == 'z') {
					b->firstnumber = 0;
				}
				if (argv[i][j] == 'e') {
					b->edgesout = 1;
				}
				if (argv[i][j] == 'm') {
					b->epsout = 1;
				}
				if (argv[i][j] == 'n') {
					b->neighbors = 1;
				}
				if (argv[i][j] == 'B') {
					b->nobound = 1;
				}
				if (argv[i][j] == 'P') {
					b->nopolywritten = 1;
				}
				if (argv[i][j] == 'N') {
					b->nonodewritten = 1;
				}
				if (argv[i][j] == 'E') {
					b->noelewritten = 1;
				}
				if (argv[i][j] == 'I') {
					b->noiterationnum = 1;
				}
				if (argv[i][j] == 'O') {
					b->noholes = 1;
				}
				if (argv[i][j] == 'C') {
					b->docheck = 1;
				}
				if (argv[i][j] == 'Q') {
					b->quiet = 1;
				}
				if (argv[i][j] == 'V') {
					b->verbose++;
				}
				if ((argv[i][j] == 'h') || (argv[i][j] == 'H') ||
					(argv[i][j] == '?')) {
						info();
				}
				if (argv[i][j] == 'v') {
					printf("WARNING: Voronoi option [v] not supported.\n");
				}
				if (argv[i][j] == 'g') {
					printf("WARNING: Geomview option [g] not supported.\n");
				}
			}
		} else {
			strncpy(b->innodefilename, argv[i], FILENAMESIZE - 1);
			b->innodefilename[FILENAMESIZE - 1] = '\0';
		}
	}
	if (b->innodefilename[0] == '\0') {
		syntax();
	}
	if (!strcmp(&b->innodefilename[strlen(b->innodefilename) - 5], ".node")) {
		b->innodefilename[strlen(b->innodefilename) - 5] = '\0';
	}
	if (!strcmp(&b->innodefilename[strlen(b->innodefilename) - 5], ".poly")) {
		b->innodefilename[strlen(b->innodefilename) - 5] = '\0';
	}
#ifndef CDT_ONLY
	if (!strcmp(&b->innodefilename[strlen(b->innodefilename) - 4], ".ele")) {
		b->innodefilename[strlen(b->innodefilename) - 4] = '\0';
	}
	if (!strcmp(&b->innodefilename[strlen(b->innodefilename) - 5], ".area")) {
		b->innodefilename[strlen(b->innodefilename) - 5] = '\0';
	}
#endif /* not CDT_ONLY */

	strcpy(b->inpolyfilename, b->innodefilename);
	strcpy(b->inelefilename, b->innodefilename);
	strcpy(b->areafilename, b->innodefilename);
	increment = 0;
	strcpy(workstring, b->innodefilename);
	j = 1;
	while (workstring[j] != '\0') {
		if ((workstring[j] == '.') && (workstring[j + 1] != '\0')) {
			increment = j + 1;
		}
		j++;
	}
	meshnumber = 0;
	if (increment > 0) {
		j = increment;
		do {
			if ((workstring[j] >= '0') && (workstring[j] <= '9')) {
				meshnumber = meshnumber * 10 + (int) (workstring[j] - '0');
			} else {
				increment = 0;
			}
			j++;
		} while (workstring[j] != '\0');
	}
	if (b->noiterationnum) {
		strcpy(b->outnodefilename, b->innodefilename);
		strcpy(b->outelefilename, b->innodefilename);
		strcpy(b->edgefilename, b->innodefilename);
		strcpy(b->epsfilename, b->innodefilename);
		strcpy(b->neighborfilename, b->innodefilename);
		strcat(b->outnodefilename, ".node");
		strcat(b->outelefilename, ".ele");
		strcat(b->edgefilename, ".edge");
		strcat(b->epsfilename, ".eps");
		strcat(b->neighborfilename, ".neigh");
	} else if (increment == 0) {
		strcpy(b->outnodefilename, b->innodefilename);
		strcpy(b->outpolyfilename, b->innodefilename);
		strcpy(b->outelefilename, b->innodefilename);
		strcpy(b->edgefilename, b->innodefilename);
		strcpy(b->epsfilename, b->innodefilename);
		strcpy(b->neighborfilename, b->innodefilename);
		strcat(b->outnodefilename, ".1.node");
		strcat(b->outpolyfilename, ".1.poly");
		strcat(b->outelefilename, ".1.ele");
		strcat(b->edgefilename, ".1.edge");
		strcat(b->epsfilename, ".1.eps");
		strcat(b->neighborfilename, ".1.neigh");
	} else {
		workstring[increment] = '%';
		workstring[increment + 1] = 'd';
		workstring[increment + 2] = '\0';
		sprintf(b->outnodefilename, workstring, meshnumber + 1);
		strcpy(b->outpolyfilename, b->outnodefilename);
		strcpy(b->outelefilename, b->outnodefilename);
		strcpy(b->edgefilename, b->outnodefilename);
		strcpy(b->epsfilename, b->outnodefilename);
		strcpy(b->neighborfilename, b->outnodefilename);
		strcat(b->outnodefilename, ".node");
		strcat(b->outpolyfilename, ".poly");
		strcat(b->outelefilename, ".ele");
		strcat(b->edgefilename, ".edge");
		strcat(b->epsfilename, ".eps");
		strcat(b->neighborfilename, ".neigh");
	}
	strcat(b->innodefilename, ".node");
	strcat(b->inpolyfilename, ".poly");
	strcat(b->inelefilename, ".ele");
	strcat(b->areafilename, ".area");
}

void check_behavior(behavior *b, behavior_legacy *lb)
{
	if (b->fixedarea && b->maxarea <= 0.0) {
		printf("ERROR: Maximum area must be greater than zero.\n");
		exit(TRI_OPTIONS);
	}

	if (b->refine && lb->noiterationnum) {
		printf("ERROR: You cannot use the -I switch when refining a triangulation.\n");
		exit(TRI_OPTIONS);
	}

	/* Regular/weighted triangulations are incompatible with PSLGs and meshing. */
	if (b->weighted && !lb->quiet && (b->poly || b->quality)) {
		printf("WARNING:  weighted triangulations (-w, -W) are incompatible\n");
		printf("  with PSLGs (-p) and meshing (-q, -a, -u).  Weights ignored.\n");
	}

	if (b->jettison && lb->nonodewritten && !lb->quiet) {
		printf("WARNING:  -j and -N switches are somewhat incompatible.\n");
		printf("  If any vertices are jettisoned, you will need the output\n");
		printf("  .node file to reconstruct the new node indices.");
	}
}

/*****************************************************************************/
/*                                                                           */
/*  statistics()   Print all sorts of cool facts.                            */
/*                                                                           */
/*****************************************************************************/

void print_statistics(context *ctx, int verbose)
{
	int i;
	REAL ratiotable[16];
	quality q;

	mesh *m = ctx->m;
	behavior *b = ctx->b;

	ratiotable[0]  =      1.5;      ratiotable[1]  =     2.0;
	ratiotable[2]  =      2.5;      ratiotable[3]  =     3.0;
	ratiotable[4]  =      4.0;      ratiotable[5]  =     6.0;
	ratiotable[6]  =     10.0;      ratiotable[7]  =    15.0;
	ratiotable[8]  =     25.0;      ratiotable[9]  =    50.0;
	ratiotable[10] =    100.0;      ratiotable[11] =   300.0;
	ratiotable[12] =   1000.0;      ratiotable[13] = 10000.0;

	printf("\nStatistics:\n\n");
	printf("  Input vertices: %d\n", m->invertices);
	if (b->refine) {
		printf("  Input triangles: %d\n", m->inelements);
	}
	if (b->poly) {
		printf("  Input segments: %d\n", m->insegments);
		if (!b->refine) {
			printf("  Input holes: %d\n", m->holes);
		}
	}

	printf("\n  Mesh vertices: %ld\n", m->vertices.items - m->undeads);
	printf("  Mesh triangles: %ld\n", m->triangles.items);
	printf("  Mesh edges: %ld\n", m->edges);
	printf("  Mesh exterior boundary edges: %ld\n", m->hullsize);
	if (b->poly || b->refine) {
		printf("  Mesh interior boundary edges: %ld\n",
			m->subsegs.items - m->hullsize);
		printf("  Mesh subsegments (constrained edges): %ld\n",
			m->subsegs.items);
	}
	printf("\n");

	if (verbose) {
		triangle_mesh_quality(ctx, &q);

		printf("Mesh quality statistics:\n\n");

		printf("  Smallest area: %16.5g   |  Largest area: %16.5g\n",
			q.smallestarea, q.biggestarea);
		printf("  Shortest edge: %16.5g   |  Longest edge: %16.5g\n",
			q.shortest, q.longest);
		printf("  Shortest altitude: %12.5g   |  Largest aspect ratio: %8.5g\n\n",
			q.minaltitude, q.worstaspect);

		printf("  Triangle aspect ratio histogram:\n");
		printf("  1.1547 - %-6.6g    :  %8d    | %6.6g - %-6.6g     :  %8d\n",
			ratiotable[0], q.aspecttable[0], ratiotable[7], ratiotable[8],
			q.aspecttable[8]);
		for (i = 1; i < 7; i++) {
			printf("  %6.6g - %-6.6g    :  %8d    | %6.6g - %-6.6g     :  %8d\n",
				ratiotable[i - 1], ratiotable[i], q.aspecttable[i],
				ratiotable[i + 7], ratiotable[i + 8], q.aspecttable[i + 8]);
		}
		printf("  %6.6g - %-6.6g    :  %8d    | %6.6g -            :  %8d\n",
			ratiotable[6], ratiotable[7], q.aspecttable[7], ratiotable[14],
			q.aspecttable[15]);
		printf("  (Aspect ratio is longest edge divided by shortest altitude)\n\n");

		printf("  Smallest angle: %15.5g   |  Largest angle: %15.5g\n\n",
			q.smallestangle, q.biggestangle);

		printf("  Angle histogram:\n");
		for (i = 0; i < 9; i++) {
			printf("    %3d - %3d degrees:  %8d    |    %3d - %3d degrees:  %8d\n",
				i * 10, i * 10 + 10, q.angletable[i],
				i * 10 + 90, i * 10 + 100, q.angletable[i + 9]);
		}
		printf("\n");

		printf("Memory allocation statistics:\n\n");
		printf("  Maximum number of vertices: %ld\n", m->vertices.maxitems);
		printf("  Maximum number of triangles: %ld\n", m->triangles.maxitems);
		if (m->subsegs.maxitems > 0) {
			printf("  Maximum number of subsegments: %ld\n", m->subsegs.maxitems);
		}
		if (m->viri.maxitems > 0) {
			printf("  Maximum number of viri: %ld\n", m->viri.maxitems);
		}
		if (m->badsubsegs.maxitems > 0) {
			printf("  Maximum number of encroached subsegments: %ld\n",
				m->badsubsegs.maxitems);
		}
		if (m->badtriangles.maxitems > 0) {
			printf("  Maximum number of bad triangles: %ld\n",
				m->badtriangles.maxitems);
		}
		if (m->flipstackers.maxitems > 0) {
			printf("  Maximum number of stacked triangle flips: %ld\n",
				m->flipstackers.maxitems);
		}
		if (m->splaynodes.maxitems > 0) {
			printf("  Maximum number of splay tree nodes: %ld\n",
				m->splaynodes.maxitems);
		}
		printf("  Approximate heap memory use (bytes): %ld\n\n",
			m->vertices.maxitems * m->vertices.itembytes +
			m->triangles.maxitems * m->triangles.itembytes +
			m->subsegs.maxitems * m->subsegs.itembytes +
			m->viri.maxitems * m->viri.itembytes +
			m->badsubsegs.maxitems * m->badsubsegs.itembytes +
			m->badtriangles.maxitems * m->badtriangles.itembytes +
			m->flipstackers.maxitems * m->flipstackers.itembytes +
			m->splaynodes.maxitems * m->splaynodes.itembytes);

		printf("Algorithmic statistics:\n\n");
		if (!b->weighted) {
			printf("  Number of incircle tests: %ld\n", m->incirclecount);
		} else {
			printf("  Number of 3D orientation tests: %ld\n", m->orient3dcount);
		}
		printf("  Number of 2D orientation tests: %ld\n", m->counterclockcount);
		if (m->hyperbolacount > 0) {
			printf("  Number of right-of-hyperbola tests: %ld\n",
				m->hyperbolacount);
		}
		if (m->circletopcount > 0) {
			printf("  Number of circle top computations: %ld\n",
				m->circletopcount);
		}
		if (m->circumcentercount > 0) {
			printf("  Number of triangle circumcenter computations: %ld\n",
				m->circumcentercount);
		}
		printf("\n");
	}
}

void reset_triangleio(triangleio *io)
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

void free_triangleio(triangleio *io)
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

void check_mesh(context *ctx)
{
	if (triangle_check_mesh(ctx) > 0) {
		printf("WARNING: mesh is topologically inconsistent.\n");
	}
	if (triangle_check_delaunay(ctx) > 0) {
		printf("WARNING: mesh is not (constraint) Delaunay.\n");
	}
}

void check(int code, const char *site)
{
	if (code < 0) {
		printf("ERROR: method 'triangle_%s' failed.\n", site);
		exit(code);
	}
}

void finishfile(FILE *file, int argc, char **argv)
{
	int i;

	fprintf(file, "# Generated by");
	for (i = 0; i < argc; i++) {
		fprintf(file, " ");
		fputs(argv[i], file);
	}
	fprintf(file, "\n");
	fclose(file);
}

int write_nodes(context *ctx, char *filename, int argc, char **argv) {
	FILE *file;
	int status = 0;

	file = fopen(filename, "w");
	if (file == (FILE *) NULL) {
		return -1;
	}

	status = triangle_write_nodes(ctx, file);

	finishfile(file, argc, argv);

	return status;
}

int write_elements(context *ctx, char *filename, int argc, char **argv) {
	FILE *file;
	int status = 0;

	file = fopen(filename, "w");
	if (file == (FILE *) NULL) {
		return -1;
	}

	status = triangle_write_elements(ctx, file);

	finishfile(file, argc, argv);

	return status;
}

int write_poly(context *ctx, char *filename, triangleio *in,
			   int argc, char **argv) {
	FILE *file;
	int status = 0;

	file = fopen(filename, "w");
	if (file == (FILE *) NULL) {
		return -1;
	}

	status = triangle_write_poly(ctx, file, in->holelist, in->numberofholes,
		in->regionlist, in->numberofregions);

	finishfile(file, argc, argv);

	return status;
}

int write_edges(context *ctx, char *filename, int argc, char **argv) {
	FILE *file;
	int status = 0;

	file = fopen(filename, "w");
	if (file == (FILE *) NULL) {
		return -1;
	}

	status = triangle_write_edges(ctx, file);

	finishfile(file, argc, argv);

	return status;
}

int write_neighbors(context *ctx, char *filename, int argc, char **argv) {
	FILE *file;
	int status = 0;

	file = fopen(filename, "w");
	if (file == (FILE *) NULL) {
		return -1;
	}

	status = triangle_write_neighbors(ctx, file);

	finishfile(file, argc, argv);

	return status;
}

int write_postscript(context *ctx, char *filename, int argc, char **argv) {
	FILE *file;
	int status = 0;

	file = fopen(filename, "w");
	if (file == (FILE *) NULL) {
		return -1;
	}

	status = triangle_write_eps(ctx, file);

	fclose(file);

	return status;
}