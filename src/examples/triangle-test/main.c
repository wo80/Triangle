/*****************************************************************************/
/*                                                                           */
/*  (tricall.c)                                                              */
/*                                                                           */
/*  Example program that demonstrates how to call Triangle.                  */
/*                                                                           */
/*  Accompanies Triangle Version 1.6                                         */
/*  July 19, 1996                                                            */
/*                                                                           */
/*  This file is placed in the public domain (but the file that it calls     */
/*  is still copyrighted!) by                                                */
/*  Jonathan Richard Shewchuk                                                */
/*  2360 Woolsey #H                                                          */
/*  Berkeley, California  94705-1927                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "triangle.h"
#include "triangle_api.h"

#include "tests.h"

/*****************************************************************************/
/*                                                                           */
/*  main()   Create and refine a mesh.                                       */
/*                                                                           */
/*****************************************************************************/

int main()
{
	context *ctx;

	assert(test_version(), "test_version");
	assert(test_context_create_destroy(), "test_context_create_destroy");
	assert(test_behavior_parse(), "test_behavior_parse");

	ctx = triangle_context_create();

	assert(test_mesh_create(ctx), "test_mesh_create");
	assert(test_mesh_copy(ctx), "test_mesh_copy");
    assert(test_mesh_statistics(ctx), "test_mesh_statistics");
    assert(test_mesh_quality(ctx), "test_mesh_quality");
    assert(test_mesh_refine(ctx), "test_mesh_refine");
	
	triangle_context_destroy(ctx);
	
	ctx = triangle_context_create();

    assert(test_mesh_load_refine(ctx), "test_mesh_load_refine");
	
	triangle_context_destroy(ctx);

	return 0;
}
