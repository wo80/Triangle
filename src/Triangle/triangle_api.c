
#include "triangle_api.h"
#include "triangle_helper.h"

#include <triangle_internal.h>

context* triangle_context_create()
{
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

void triangle_context_destroy(context* ctx)
{
	if (check_context(ctx) < 0) {
		return;
	}

	triangledeinit(ctx->m, ctx->b);

	free(ctx->b);
	free(ctx->m);

	//free(ctx);
}

int triangle_context_options(context* ctx, char *options)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	parsecommandline(options, ctx->b);

	return check_behavior(ctx->b);
}

void triangle_context_get_behavior(context* ctx, behavior *out)
{
	if (check_context(ctx) < 0) {
		return;
	}

	*out = *ctx->b;
}

int triangle_context_set_behavior(context* ctx, behavior *in)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	*ctx->b = *in;

	behavior_update(ctx->b);

	return check_behavior(in);
}

int triangle_mesh_quality(context* ctx, quality *q)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	return quality_statistics(ctx->m, ctx->b, q);
}

int triangle_mesh_statistics(context* ctx, statistics *s)
{
	mesh *m;
	behavior *b;

	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	m = ctx->m;
	b = ctx->b;

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

	s->bounds.xmin = m->xmin;
	s->bounds.ymin = m->ymin;
	s->bounds.xmax = m->xmax;
	s->bounds.ymax = m->ymax;

	return 0;
}

int triangle_memory(context* ctx)
{
	mesh *m;

	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	m = ctx->m;

	return m->vertices.maxitems * m->vertices.itembytes +
		m->triangles.maxitems * m->triangles.itembytes +
		m->subsegs.maxitems * m->subsegs.itembytes +
		m->viri.maxitems * m->viri.itembytes +
		m->badsubsegs.maxitems * m->badsubsegs.itembytes +
		m->badtriangles.maxitems * m->badtriangles.itembytes +
		m->flipstackers.maxitems * m->flipstackers.itembytes +
		m->splaynodes.maxitems * m->splaynodes.itembytes;
}

int triangle_check_mesh(context *ctx)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	return checkmesh(ctx->m, ctx->b);
}

int triangle_check_delaunay(context *ctx)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	return checkdelaunay(ctx->m, ctx->b);
}

int triangle_mesh_create(context* ctx, triangleio *in)
{
	mesh *m;
	behavior *b;

	int status = 0;

	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	m = ctx->m;
	b = ctx->b;

	/* Mesh already created. */
	if (m->triangles.items > 0) {
		return TRI_FAILURE;
	}

	status = transfernodes(m, b, in->pointlist, in->pointattributelist,
		in->pointmarkerlist, in->numberofpoints,
		in->numberofpointattributes);

	if (status < 0) {
		return status;
	}

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
			in->segmentmarkerlist, in->numberofsegments, &status);
		if (status < 0) {
			return status;
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
		enforcequality(m, b, &status);           
		if (status < 0) {
			return status;
		}
	}
#endif

	/* Calculate the number of edges. */
	m->edges = (3l * m->triangles.items + m->hullsize) / 2l;

	if (b->order > 1) {
		highorder(m, b); /* Promote elements to higher polynomial order. */
	}

	return status;
}

int triangle_mesh_load(context* ctx, triangleio *in)
{
	mesh *m;
	behavior *b;

	int status = 0;

	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	m = ctx->m;
	b = ctx->b;
	
	/* Mesh already created. */
	if (m->triangles.items > 0) {
		return TRI_FAILURE;
	}

	//if (!b->refine) {
	//   ... don't need to check. calling this method implies the option.
	//}

	status = transfernodes(m, b, in->pointlist, in->pointattributelist,
		in->pointmarkerlist, in->numberofpoints,
		in->numberofpointattributes);

	if (status < 0) {
		return status;
	}

	/* Read and reconstruct a mesh. */
	m->hullsize = reconstruct(m, b, in->trianglelist,
		in->triangleattributelist, in->trianglearealist,
		in->numberoftriangles, in->numberofcorners,
		in->numberoftriangleattributes,
		in->segmentlist, in->segmentmarkerlist,
		in->numberofsegments);

	if (m->hullsize < 0l) {
		return (int) m->hullsize;
	}

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

	return status;
}

int triangle_mesh_refine(context* ctx)
{
	mesh *m;
	behavior *b;

	int status = 0;

	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	m = ctx->m;
	b = ctx->b;

#ifndef CDT_ONLY
	if (b->quality && (m->triangles.items > 0)) {
		/* Enforce angle and area constraints. */
		enforcequality(m, b, &status);           
		if (status < 0) {
			return status;
		}

		/* Calculate the number of edges. */
		m->edges = (3l * m->triangles.items + m->hullsize) / 2l;
	}
#endif

	return status;
}

int triangle_mesh_copy(context* ctx, triangleio *out,
					   int edges, int neighbors)
{
	mesh *m;
	behavior *b;

	int status = 0;

	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	m = ctx->m;
	b = ctx->b;

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

	return status;
}

void triangle_free(VOID *memptr)
{
	free(memptr);
}

#ifndef NO_FILE_IO

int triangle_write_nodes(context *ctx, FILE *file)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	return file_writenodes(ctx->m, ctx->b, file);
}

int triangle_write_elements(context *ctx, FILE *file)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	return file_writeelements(ctx->m, ctx->b, file);
}

int triangle_write_poly(context *ctx, FILE *file,
						REAL *holelist, int holes, REAL *regionlist, int regions)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	return file_writepoly(ctx->m, ctx->b, file,
		holelist, holes, regionlist, regions);
}

int triangle_write_edges(context *ctx, FILE *file)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	return file_writeedges(ctx->m, ctx->b, file);
}

int triangle_write_neighbors(context *ctx, FILE *file)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	return file_writeneighbors(ctx->m, ctx->b, file);
}

int triangle_write_eps(context *ctx, FILE *file)
{
	if (check_context(ctx) < 0) {
		return TRI_NULL;
	}

	return file_write_eps(ctx->m, ctx->b, file);
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