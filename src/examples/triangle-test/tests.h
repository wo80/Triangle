#ifndef TESTS_H
#define TESTS_H

#include "triangle.h"
#include "triangle_api.h"

#define SUCCESS 1
#define FAILURE 0

void assert(int result, char *message);

int test_context_create_destroy();

int test_behavior_parse();

int test_mesh_create(context *ctx);

int test_mesh_copy(context *ctx);

int test_mesh_statistics(context *ctx);

int test_mesh_quality(context *ctx);

int test_mesh_refine(context *ctx);

int test_mesh_load_refine(context *ctx);

#endif /* TESTS_H */