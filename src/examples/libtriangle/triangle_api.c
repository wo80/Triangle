
#include "triangle_api.h"
#include "triangle_helper.h"
#include <triangle_internal.h>

int triangle_behavior_parse(behavior *b, char *options)
{
	parsecommandline(options, b);

	return check_behavior(b);
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
	parsecommandline("\0", b);

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
	
	m->steinerleft = b->steiner;
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

	// TODO: check for error (hullsize < 0)

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

int triangle_mesh_copy(context* ctx, triangleio *out,
	int edges, int neighbors)
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

	/* writenodes() numbers the vertices too. */
	writenodes(m, b, &out->pointlist, &out->pointattributelist,
		&out->pointmarkerlist);

	/* Always write elements. */
	writeelements(m, b, &out->trianglelist, &out->triangleattributelist);

	/* The -c switch (convex switch) causes a PSLG to be written */
	/*   even if none was read.                                  */
	if (b->poly || b->convex) {
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

	if (edges) {
		writeedges(m, b, &out->edgelist, &out->edgemarkerlist);
	}

	if (neighbors) {
		writeneighbors(m, b, &out->neighborlist);
	}

	return result;
}

#ifndef NO_FILE_IO

int triangle_write_nodes(context *ctx, FILE *file)
{
	return file_writenodes(ctx->m, ctx->b, file);
}

int triangle_write_elements(context *ctx, FILE *file)
{
	return file_writeelements(ctx->m, ctx->b, file);
}

int triangle_write_poly(context *ctx, FILE *file,
						REAL *holelist, int holes, REAL *regionlist, int regions)
{
	return file_writepoly(ctx->m, ctx->b, file,
		holelist, holes, regionlist, regions);
}

int triangle_write_edges(context *ctx, FILE *file)
{
	return file_writeedges(ctx->m, ctx->b, file);
}

int triangle_write_neighbors(context *ctx, FILE *file)
{
	return file_writeneighbors(ctx->m, ctx->b, file);
}

int triangle_read_nodes(const char* filename, triangleio *io, int *firstnode)
{
	int status;
	FILE *file = fopen(filename, "r");

	if (file == (FILE *) NULL) {
		return TRI_FILE_OPEN;
	}

	status = file_readnodes(file, io, firstnode);

    fclose(file);

	return status;
}

int triangle_read_poly(const char* filename, triangleio *io, int *firstnode)
{
	int status;
	FILE *file = fopen(filename, "r");

	if (file == (FILE *) NULL) {
		return TRI_FILE_OPEN;
	}

	status = file_readpoly(file, io, firstnode);

    fclose(file);

	return status;
}

int triangle_read_elements(const char* filename, triangleio *io)
{
	int status;
	FILE *file = fopen(filename, "r");

	if (file == (FILE *) NULL) {
		return TRI_FILE_OPEN;
	}

	status = file_readelements(file, io);

    fclose(file);

	return status;
}

int triangle_read_area(const char* filename, triangleio *io)
{
	int status;
	FILE *file = fopen(filename, "r");

	if (file == (FILE *) NULL) {
		return TRI_FILE_OPEN;
	}

	status = file_readelementsarea(file, io);

    fclose(file);

	return status;
}

#endif /* NO_FILE_IO */