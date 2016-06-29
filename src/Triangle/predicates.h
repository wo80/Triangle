#ifndef PREDICATES_H
#define PREDICATES_H

#include "triangle_config.h"
#include "triangle_core.h"

void exactinit();

REAL counterclockwise(mesh *m, behavior *b,
                      vertex pa, vertex pb, vertex pc);

REAL incircle(mesh *m, behavior *b,
              vertex pa, vertex pb, vertex pc, vertex pd);

REAL nonregular(mesh *m, behavior *b,
                vertex pa, vertex pb, vertex pc, vertex pd);

void findcircumcenter(mesh *m, behavior *b,
                      vertex torg, vertex tdest, vertex tapex,
                      vertex circumcenter, REAL *xi, REAL *eta, int offcenter);

#endif /* PREDICATES_H */