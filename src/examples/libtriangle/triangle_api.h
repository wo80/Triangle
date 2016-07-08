#ifndef TRIANGLE_API_H
#define TRIANGLE_API_H

#include <triangle.h>

#include <stdio.h>
#include <stdlib.h>

#define EXPORT __declspec(dllexport)

#ifdef __cplusplus
extern "C" {
#endif

	typedef struct context_t {
		mesh *m;
		behavior *b;
	} context;

	typedef struct statistics_t {
		int vertices;
		int undeads;
		int triangles;
		int hullsize;
		int subsegs;
		int edges;
		int memory;
	} statistics;
	
	/**
	 * Create a context struct.
	 * @return Pointer to context struct.
	 */
	EXPORT context* triangle_context_create();
	
	/**
	 * Destroy a context struct.
	 * @param ctx Pointer to context struct.
	 */
	EXPORT VOID triangle_context_destroy(context* ctx);
	
	/**
	 * Parse Triangle options string.
	 * @param b Pointer to behavior struct.
	 * @param options String containing Triangle options.
	 * @return Integer status code.
	 */
	EXPORT int triangle_behavior_parse(behavior *b, char *options);
	
	/**
	 * Triangulate an input polygon.
	 * @param ctx Pointer to context struct.
	 * @param in Pointer to triangleio struct, containing the input polygon.
	 * @return Integer status code.
	 */
	EXPORT int triangle_mesh_create(context* ctx, triangleio *in);
	
	/**
	 * Load a previously computed mesh.
	 * @param ctx Pointer to context struct.
	 * @param in Pointer to triangleio struct, containing a previously computed mesh.
	 * @return Integer status code.
	 */
	EXPORT int triangle_mesh_load(context* ctx, triangleio *in);
	
	/**
	 * Refine the current mesh.
	 * @param ctx Pointer to context struct.
	 * @return Integer status code.
	 */
	EXPORT int triangle_mesh_refine(context* ctx);
	
	/**
	 * Copy mesh to triangleio struct.
	 * @param ctx Pointer to context struct.
	 * @param out Pointer to triangleio struct, containing the mesh.
	 * @param edges If non-zero, write edges to triangleio.
	 * @param neighbors If non-zero, write neighbors to triangleio.
	 * @return Integer status code.
	 */
	EXPORT int triangle_mesh_copy(context* ctx, triangleio *out,
		int edges, int neighbors);
	
	/**
	 * Compute mesh statistics.
	 * @param ctx Pointer to context struct.
	 * @param s Pointer to statistics struct.
	 * @return Integer status code.
	 */
	EXPORT int triangle_mesh_statistics(context *ctx, statistics *s);
	
	/**
	 * Compute mesh quality statistics.
	 * @param ctx Pointer to context struct.
	 * @param q Pointer to quality struct.
	 * @return Integer status code.
	 */
	EXPORT int triangle_mesh_quality(context *ctx, quality *q);

#ifndef NO_FILE_IO
	/**
	 * Write node file.
	 * @param ctx Pointer to context struct.
	 * @param file File pointer (to .node file).
	 * @return Integer status code.
	 */
	EXPORT int triangle_write_nodes(context *ctx, FILE *file);
	
	/**
	 * Write ele file.
	 * @param ctx Pointer to context struct.
	 * @param file File pointer (to .ele file).
	 * @return Integer status code.
	 */
	EXPORT int triangle_write_elements(context *ctx, FILE *file);
	
	/**
	 * Write poly file.
	 * @param ctx Pointer to context struct.
	 * @param file File pointer (to .poly file).
	 * @param holelist Pointer to an array of REALs.
	 * @param holes Number of holes.
	 * @param regionlist Pointer to an array of REALs.
	 * @param regions Number of regions.
	 * @return Integer status code.
	 */
	EXPORT int triangle_write_poly(context *ctx, FILE *file,
		REAL *holelist, int holes, REAL *regionlist, int regions);
	
	/**
	 * Write edge file.
	 * @param ctx Pointer to context struct.
	 * @param file File pointer (to .edge file).
	 * @return Integer status code.
	 */
	EXPORT int triangle_write_edges(context *ctx, FILE *file);
	
	/**
	 * Write neighbor file.
	 * @param ctx Pointer to context struct.
	 * @param file File pointer (to .neigh file).
	 * @return Integer status code.
	 */
	EXPORT int triangle_write_neighbors(context *ctx, FILE *file);
	
	/**
	 * Read nodes file.
	 * @param file File pointer (to .node file).
	 * @param ctx Pointer to triangleio struct.
	 * @return Integer status code.
	 */
	EXPORT int triangle_read_nodes(FILE *file, triangleio *io);
	
	/**
	 * Read poly file.
	 * @param file File pointer (to .poly file).
	 * @param ctx Pointer to triangleio struct.
	 * @return Integer status code.
	 */
	EXPORT int triangle_read_poly(FILE *file, triangleio *io);
	
	/**
	 * Read elements file.
	 * @param file File pointer (to .ele file).
	 * @param ctx Pointer to triangleio struct.
	 * @return Integer status code.
	 */
	EXPORT int triangle_read_elements(FILE *file, triangleio *io);
#endif /* NO_FILE_IO */

#ifdef __cplusplus
}
#endif

#endif /* TRIANGLE_API_H */