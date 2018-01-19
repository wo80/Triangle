#ifndef TRIANGLE_HELPER_H
#define TRIANGLE_HELPER_H

#include "../triangle_api.h"

int check_context(context* c);
int check_behavior(behavior* b);
int check_triangleio(triangleio* io, int firstnumber);

int restore_pointmarkers(context* ctx, int *pointmarkers);

#endif /* TRIANGLE_HELPER_H */
