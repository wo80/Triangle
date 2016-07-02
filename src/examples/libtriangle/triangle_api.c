
#include "triangle_api.h"
#include <triangle_internal.h>

context* triangle_context_create()
{
	mesh *m = malloc(sizeof *m);
	behavior *b = malloc(sizeof *b);

	context *ctx = malloc(sizeof *ctx);

	ctx->m = m;
	ctx->b = b;

	triangleinit(ctx->m);

	return ctx;
}

VOID triangle_context_destory(context* ctx)
{
	triangledeinit(ctx->m, ctx->b);

	free(ctx->b);
	free(ctx->m);
	free(ctx);
}

int triangle_quality_statistics(context* ctx, statistics *stats)
{
	return quality_statistics(ctx->m, ctx->b, stats);
}

int triangle_options(context* ctx, char *options)
{
	int result = 0;

	parsecommandline(1, &options, ctx->b, &result);

	ctx->m->steinerleft = ctx->b->steiner;

	return result;
}

int triangle_mesh(context* ctx, triangleio *in)
{
	mesh *m = ctx->m;
	behavior *b = ctx->b;

	int result = 0;

	transfernodes(m, b, in->pointlist, in->pointattributelist,
		in->pointmarkerlist, in->numberofpoints,
		in->numberofpointattributes);

#ifdef CDT_ONLY
	m->hullsize = delaunay(m, b);                /* Triangulate the vertices. */
#else /* not CDT_ONLY */
	if (b->refine) {
		/* Read and reconstruct a mesh. */
		m->hullsize = reconstruct(m, b, in->trianglelist,
			in->triangleattributelist, in->trianglearealist,
			in->numberoftriangles, in->numberofcorners,
			in->numberoftriangleattributes,
			in->segmentlist, in->segmentmarkerlist,
			in->numberofsegments);
	} else {
		m->hullsize = delaunay(m, b);              /* Triangulate the vertices. */
	}
#endif /* not CDT_ONLY */

	/* Ensure that no vertex can be mistaken for a triangular bounding */
	/*   box vertex in insertvertex().                                 */
	m->infvertex1 = (vertex) NULL;
	m->infvertex2 = (vertex) NULL;
	m->infvertex3 = (vertex) NULL;

	if (b->usesegments) {
		m->checksegments = 1;                /* Segments will be introduced next. */
		if (!b->refine) {
			/* Insert PSLG segments and/or convex hull segments. */
			formskeleton(m, b, in->segmentlist,
				in->segmentmarkerlist, in->numberofsegments, &result);
			if (result > 0) {
				triangledeinit(m, b); /* TODO: triangledeinit ok? */
				return result;
			}
		}
	}

	if (b->poly && (m->triangles.items > 0)) {
		m->holes = in->numberofholes;
		m->regions = in->numberofregions;

		if (!b->refine) {
			/* Carve out holes and concavities. */
			carveholes(m, b, in->holelist, m->holes, in->regionlist, m->regions);
		}
	} else {
		/* Without a PSLG, there can be no holes or regional attributes   */
		/*   or area constraints.  The following are set to zero to avoid */
		/*   an accidental free() later.                                  */
		m->holes = 0;
		m->regions = 0;
	}

#ifndef CDT_ONLY
	if (b->quality && (m->triangles.items > 0)) {
		/* Enforce angle and area constraints. */
		enforcequality(m, b, &result);           
		if (result > 0) {
			triangledeinit(m, b); /* TODO: triangledeinit ok? */
			return result;
		}
	}
#endif

	/* Calculate the number of edges. */
	m->edges = (3l * m->triangles.items + m->hullsize) / 2l;

	return result;
}

int triangle_refine(context* ctx)
{
	mesh *m = ctx->m;
	behavior *b = ctx->b;

	int result = 0;

#ifndef CDT_ONLY
	if (b->quality && (m->triangles.items > 0)) {
		/* Enforce angle and area constraints. */
		enforcequality(m, b, &result);           
		if (result > 0) {
			triangledeinit(m, b); /* TODO: triangledeinit ok? */
			return result;
		}
	}

	/* Calculate the number of edges. */
	m->edges = (3l * m->triangles.items + m->hullsize) / 2l;
#endif

	return result;
}

int triangle_output(context* ctx, triangleio *out)
{
	mesh *m = ctx->m;
	behavior *b = ctx->b;

	int result = 0;

	if (b->jettison) {
		out->numberofpoints = m->vertices.items - m->undeads;
	} else {
		out->numberofpoints = m->vertices.items;
	}
	out->numberofpointattributes = m->nextras;
	out->numberoftriangles = m->triangles.items;
	out->numberofcorners = (b->order + 1) * (b->order + 2) / 2;
	out->numberoftriangleattributes = m->eextras;
	out->numberofedges = m->edges;
	if (b->usesegments) {
		out->numberofsegments = m->subsegs.items;
	} else {
		out->numberofsegments = m->hullsize;
	}

	/* If not using iteration numbers, don't write a .node file if one was */
	/*   read, because the original one would be overwritten!              */
	if (b->nonodewritten || (b->noiterationnum && m->readnodefile)) {
		numbernodes(m, b);         /* We must remember to number the vertices. */
	} else {
		/* writenodes() numbers the vertices too. */
		writenodes(m, b, &out->pointlist, &out->pointattributelist,
			&out->pointmarkerlist);
	}

	if (!b->noelewritten) {
		writeelements(m, b, &out->trianglelist, &out->triangleattributelist);
	}

	/* The -c switch (convex switch) causes a PSLG to be written */
	/*   even if none was read.                                  */
	if (b->poly || b->convex) {
		/* If not using iteration numbers, don't overwrite the .poly file. */
		if (b->nopolywritten || b->noiterationnum) {
		} else {
			writepoly(m, b, &out->segmentlist, &out->segmentmarkerlist);
			out->numberofholes = m->holes;
			out->numberofregions = m->regions;
			if (b->poly) {
				//out->holelist = in->holelist;
				//out->regionlist = in->regionlist;
			} else {
				out->holelist = (REAL *) NULL;
				out->regionlist = (REAL *) NULL;
			}
		}
	}

	if (b->edgesout) {
		writeedges(m, b, &out->edgelist, &out->edgemarkerlist);
	}

	if (b->neighbors) {
		writeneighbors(m, b, &out->neighborlist);
	}

	return result;
}

int triangle_write_nodes(context *ctx, FILE *nodefile)
{
	return file_writenodes(ctx->m, ctx->b, nodefile);
}

int triangle_write_elements(context *ctx, FILE *elefile)
{
	return file_writeelements(ctx->m, ctx->b, elefile);
}

int triangle_write_poly(context *ctx, FILE *polyfile,
		REAL *holelist, int holes, REAL *regionlist, int regions)
{
	return file_writepoly(ctx->m, ctx->b, polyfile,
				   holelist, holes, regionlist, regions);
}

int triangle_write_edges(context *ctx, FILE *edgefile)
{
	return file_writeedges(ctx->m, ctx->b, edgefile);
}

int triangle_write_neighbors(context *ctx, FILE *neighborfile)
{
	return file_writeneighbors(ctx->m, ctx->b, neighborfile);
}