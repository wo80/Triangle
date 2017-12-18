
#include <stdio.h>
#include <stdlib.h>

#include "tests.h"
#include "util.h"

void assert(int result, char *message)
{
	if(result) {
		printf("OK: %s\n", message);
	} else {
		printf("FAILED: %s\n", message);
		exit(1);
	}
}

int test_context_create_destroy()
{
	context *ctx = triangle_context_create();

	triangle_context_destroy(ctx);

	return 1;
}

int test_behavior_parse()
{
	context *ctx = triangle_context_create();
	behavior *b;
	
	int result = SUCCESS;
	
	triangle_context_options(ctx, "prq25U100a2.5lcznYYS100");

	b = ctx->b;

	result &= b->poly;
	result &= b->refine;
	result &= b->quality;
	result &= (b->minangle == 25.0);
#ifndef NO_ACUTE
    result &= (b->maxangle == 100.0);
#endif
	result &= b->fixedarea;
	result &= (b->maxarea == 2.5);
	result &= !b->dwyer;
	result &= b->convex;
	result &= (b->firstnumber == 0);
	result &= b->neighbors;
	result &= (b->nobisect == 2);
	result &= (b->steiner == 100);
	
	return result;
}

int test_mesh_create(context *ctx)
{
	triangleio in;
	mesh *m;

	int result = SUCCESS;

	create_rectangle(&in, 0.0, 10.0, 1.0, 0.0);

	triangle_context_options(ctx, "pczAn");

	triangle_mesh_create(ctx, &in);

	m = ctx->m;

	result &= (m->invertices == 4);
	result &= (m->insegments == 0);
	result &= (m->triangles.items == 2);
	result &= (m->vertices.items == 4);
	result &= (m->subsegs.items == 4);
	result &= (m->regions == 1);
	result &= (m->edges == 5);
	result &= (m->nextras == 1);
	result &= (m->eextras == 1);
	result &= (m->hullsize == 4);
	result &= (m->steinerleft == -1);

	triangleio_free(&in);

	return result;
}

int test_mesh_copy(context *ctx)
{
	triangleio out;

	int result = SUCCESS;
	
	triangleio_reset(&out);

	triangle_mesh_copy(ctx, &out, 1, 1);

	result &= (out.numberofpoints == 4);
	result &= (out.numberofpointattributes == 1);
	result &= (out.numberoftriangles == 2);
	result &= (out.numberoftriangleattributes == 1);
	result &= (out.numberofsegments == 4);
	result &= (out.numberofedges == 5);
	result &= (out.numberofregions == 1);

	triangleio_free(&out);

	return result;
}

int test_mesh_statistics(context *ctx)
{
	statistics s;
	
	int result = SUCCESS;
	
	triangle_mesh_statistics(ctx, &s);
	
	result &= (s.vertices == 4);
	result &= (s.triangles == 2);
	result &= (s.hullsize == 4);
	result &= (s.subsegs == 4);
	result &= (s.edges == 5);

	return result;
}

int test_mesh_quality(context *ctx)
{
	quality q;

	int result = SUCCESS;
	
	triangle_mesh_quality(ctx, &q);

	return result;
}

int test_mesh_refine(context *ctx)
{
	mesh *m;
	
	int result = SUCCESS;

	ctx->b->fixedarea = 1;
	ctx->b->maxarea = 2.5;
	ctx->b->quality = 1;

	triangle_mesh_refine(ctx);

	m = ctx->m;

	result &= (m->triangles.items == 4);
	result &= (m->vertices.items == 5);
	result &= (m->subsegs.items == 4);
	result &= (m->edges == 8);
	result &= (m->hullsize == 4);

	return result;
}

int test_mesh_load_refine(context *ctx)
{
	triangleio in;
	mesh *m;

	int result = SUCCESS;

	create_rectangle_mesh(&in, 0.0, 10.0, 1.0, 0.0);
	
	triangle_context_options(ctx, "prazBP");

	triangle_mesh_load(ctx, &in);

	m = ctx->m;

	result &= (m->invertices == 4);
	result &= (m->insegments == 4);
	result &= (m->triangles.items == 2);
	result &= (m->vertices.items == 4);
	result &= (m->subsegs.items == 4);
	result &= (m->regions == 1);
	result &= (m->edges == 5);
	result &= (m->nextras == 1);
	result &= (m->eextras == 1);
	result &= (m->hullsize == 4);
	result &= (m->steinerleft == -1);
	
	// Refine the triangulation according to the attached
	// triangle area constraints.

	result &= (triangle_mesh_refine(ctx) == 0);

	result &= (m->triangles.items == 16);
	result &= (m->vertices.items == 14);
	result &= (m->subsegs.items == 10);
	result &= (m->edges == 29);
	result &= (m->hullsize == 10);

	triangleio_free(&in);

	return result;
}
