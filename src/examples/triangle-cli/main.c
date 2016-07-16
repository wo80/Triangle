
#include <stdio.h>
#include <stdlib.h>
#include <triangle_api.h>

#include "util.h"

int main(int argc, char **argv)
{
	context *ctx;
	behavior *b;

	behavior_legacy lb;
	triangleio in;

	int firstnode;

	parsecommandline_legacy(argc, argv, &lb);

	ctx = triangle_context_create();
	b = ctx->b;

	reset_triangleio(&in);

	/* At this point argc is guaranteed to be > 1. */
	if (argv[1][0] == '-') {
		triangle_behavior_parse(b, argv[1] + 1);
	}

	check_behavior(b, &lb);

	if (b->poly) {
		check(triangle_read_poly(lb.inpolyfilename, &in, &firstnode), "read_poly");
	}
	
	if (in.numberofpoints == 0) {
		/* The poly file had no vertices. */
		check(triangle_read_nodes(lb.innodefilename, &in, &firstnode), "read_nodes");
	}

#ifdef CDT_ONLY
	check(triangle_mesh_create(ctx, &in), "mesh_create"); /* Triangulate the vertices. */
#else /* not CDT_ONLY */
	if (b->refine) {
		/* Read the mesh. */
		check(triangle_read_elements(lb.inelefilename, &in), "read_elements");
		if (b->vararea) {
			check(triangle_read_area(lb.areafilename, &in), "read_area");
		}
		/* Reconstruct the mesh. */
		check(triangle_mesh_load(ctx, &in), "mesh_load");
		check(triangle_mesh_refine(ctx), "mesh_refine");
	} else {
		/* Triangulate the polygon. */
		check(triangle_mesh_create(ctx, &in), "mesh_create");
	}
#endif /* not CDT_ONLY */

	/* If not using iteration numbers, don't write a .node file if one was */
	/*   read, because the original one would be overwritten!              */
	if (lb.nonodewritten || (lb.noiterationnum && !b->poly)) {
		if (!lb.quiet) {
			printf("NOT writing a .node file.\n");
		}
		//numbernodes(&m, &b);         /* We must remember to number the vertices. */
	} else {
		/* writenodes() numbers the vertices too. */
		write_nodes(ctx, lb.outnodefilename, argc, argv);
	}
	if (lb.noelewritten) {
		if (!lb.quiet) {
			printf("NOT writing an .ele file.\n");
		}
	} else {
		write_elements(ctx, lb.outelefilename, argc, argv);
	}
	/* The -c switch (convex switch) causes a PSLG to be written */
	/*   even if none was read.                                  */
	if (b->poly || b->convex) {
		/* If not using iteration numbers, don't overwrite the .poly file. */
		if (lb.nopolywritten || lb.noiterationnum) {
			if (!lb.quiet) {
				printf("NOT writing a .poly file.\n");
			}
		} else {
			write_poly(ctx, lb.outpolyfilename, &in, argc, argv);
		}
	}

	if (lb.edgesout) {
		write_edges(ctx, lb.edgefilename, argc, argv);
	}

	if (lb.epsout) {
		write_postscript(ctx, lb.epsfilename, argc, argv);
	}

	if (lb.neighbors) {
		write_neighbors(ctx, lb.neighborfilename, argc, argv);
	}

	if (!lb.quiet) {
		print_statistics(ctx, lb.verbose);
	}

#ifndef REDUCED
	if (lb.docheck) {
		check_mesh(ctx);
	}
#endif /* not REDUCED */

	triangle_context_destroy(ctx);

	return 0;
}
