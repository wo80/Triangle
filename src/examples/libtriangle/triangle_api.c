
#include "triangle_api.h"
#include <triangle_internal.h>

int triangle_behavior_parse(behavior *b, char *options)
{
	int result = 0;

	parsecommandline(1, &options, b, &result);

	return result;
}

context* triangle_context_create()
{
	int result = 0;

	mesh *m = malloc(sizeof *m);
	behavior *b = malloc(sizeof *b);

	context *ctx = malloc(sizeof *ctx);

	ctx->m = m;
	ctx->b = b;
	
	/* Initialize default behavior values. */
	parsecommandline(0, (char **)NULL, b, &result);

	triangleinit(ctx->m);

	m->steinerleft = b->steiner;

	/* Initialize some mesh pointers to zero. */

	m->lastflip = NULL;
	m->infvertex1 = NULL;
	m->infvertex2 = NULL;
	m->infvertex3 = NULL;
	m->dummytri = NULL;
	m->dummytribase = NULL;
	m->dummysub = NULL;
	m->dummysubbase = NULL;

	return ctx;
}

VOID triangle_context_destroy(context* ctx)
{
	triangledeinit(ctx->m, ctx->b);

	free(ctx->b);
	free(ctx->m);

	// TODO: free(ctx)
}

int triangle_mesh_quality(context* ctx, quality *q)
{
	return quality_statistics(ctx->m, ctx->b, q);
}

int triangle_mesh_statistics(context* ctx, statistics *s)
{
	mesh *m = ctx->m;
	behavior *b = ctx->b;

	s->vertices = m->vertices.items;
	s->undeads = m->undeads;
	s->triangles = m->triangles.items;
	s->hullsize = m->hullsize;
	s->edges = m->edges;

	if (b->poly || b->refine) {
		s->subsegs = m->subsegs.items;
	} else {
		s->subsegs = 0;
	}

	s->memory = m->vertices.maxitems * m->vertices.itembytes +
		m->triangles.maxitems * m->triangles.itembytes +
		m->subsegs.maxitems * m->subsegs.itembytes +
		m->viri.maxitems * m->viri.itembytes +
		m->badsubsegs.maxitems * m->badsubsegs.itembytes +
		m->badtriangles.maxitems * m->badtriangles.itembytes +
		m->flipstackers.maxitems * m->flipstackers.itembytes +
		m->splaynodes.maxitems * m->splaynodes.itembytes;

	return 0;
}

int triangle_mesh_create(context* ctx, triangleio *in)
{
	mesh *m = ctx->m;
	behavior *b = ctx->b;

	int result = 0;

	transfernodes(m, b, in->pointlist, in->pointattributelist,
		in->pointmarkerlist, in->numberofpoints,
		in->numberofpointattributes);

	m->hullsize = delaunay(m, b); /* Triangulate the vertices. */

	/* Ensure that no vertex can be mistaken for a triangular bounding */
	/*   box vertex in insertvertex().                                 */
	m->infvertex1 = (vertex) NULL;
	m->infvertex2 = (vertex) NULL;
	m->infvertex3 = (vertex) NULL;

	if (b->usesegments) {
		m->checksegments = 1; /* Segments will be introduced next. */

		/* Insert PSLG segments and/or convex hull segments. */
		formskeleton(m, b, in->segmentlist,
			in->segmentmarkerlist, in->numberofsegments, &result);
		if (result > 0) {
			triangledeinit(m, b); /* TODO: triangledeinit ok? */
			return result;
		}
	}

	if (b->poly && (m->triangles.items > 0)) {
		m->holes = in->numberofholes;
		m->regions = in->numberofregions;

		/* Carve out holes and concavities. */
		carveholes(m, b, in->holelist, m->holes, in->regionlist, m->regions);
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

int triangle_mesh_load(context* ctx, triangleio *in)
{
	mesh *m = ctx->m;
	behavior *b = ctx->b;

	int result = 0;

	//if (!b->refine) {
	//   ... don't need to check. calling this method implies the option.
	//}

	transfernodes(m, b, in->pointlist, in->pointattributelist,
		in->pointmarkerlist, in->numberofpoints,
		in->numberofpointattributes);

	/* Read and reconstruct a mesh. */
	m->hullsize = reconstruct(m, b, in->trianglelist,
		in->triangleattributelist, in->trianglearealist,
		in->numberoftriangles, in->numberofcorners,
		in->numberoftriangleattributes,
		in->segmentlist, in->segmentmarkerlist,
		in->numberofsegments);

	/* Ensure that no vertex can be mistaken for a triangular bounding */
	/*   box vertex in insertvertex().                                 */
	m->infvertex1 = (vertex) NULL;
	m->infvertex2 = (vertex) NULL;
	m->infvertex3 = (vertex) NULL;

	if (b->usesegments) {
		m->checksegments = 1; /* Segments will be introduced next. */
	}

	if (b->poly && (m->triangles.items > 0)) {
		m->holes = in->numberofholes;
		m->regions = in->numberofregions;
	} else {
		/* Without a PSLG, there can be no holes or regional attributes   */
		/*   or area constraints.  The following are set to zero to avoid */
		/*   an accidental free() later.                                  */
		m->holes = 0;
		m->regions = 0;
	}

	/* Calculate the number of edges. */
	m->edges = (3l * m->triangles.items + m->hullsize) / 2l;

	return result;
}

int triangle_mesh_refine(context* ctx)
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

int triangle_mesh_copy(context* ctx, triangleio *out)
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

#ifndef NO_FILE_IO

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

#endif /* NO_FILE_IO */