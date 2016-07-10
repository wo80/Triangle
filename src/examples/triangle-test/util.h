#ifndef UTIL_H
#define UTIL_H

#include "triangle.h"
#include "triangle_api.h"

void triangleio_reset(triangleio *io);

void triangleio_free(triangleio *io);

void create_rectangle(triangleio *io, REAL left, REAL top, REAL right, REAL bottom);

void create_rectangle_mesh(triangleio *io, REAL left, REAL top, REAL right, REAL bottom);

#endif /* UTIL_H */