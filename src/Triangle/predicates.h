#ifndef PREDICATES_H
#define PREDICATES_H

#include "triangle.h"

void exactinit();

REAL counterclockwise(struct mesh *m, struct behavior *b,
                      vertex pa, vertex pb, vertex pc);

REAL incircle(struct mesh *m, struct behavior *b,
              vertex pa, vertex pb, vertex pc, vertex pd);

REAL nonregular(struct mesh *m, struct behavior *b,
                vertex pa, vertex pb, vertex pc, vertex pd);

void findcircumcenter(struct mesh *m, struct behavior *b,
                      vertex torg, vertex tdest, vertex tapex,
                      vertex circumcenter, REAL *xi, REAL *eta, int offcenter);

#endif /* PREDICATES_H */