#ifndef UTIL_H
#define UTIL_H

#include "triangle_api.h"

/* Data structure for command line switches and file names.                  */

typedef struct {

  int firstnumber;
  int edgesout, epsout, neighbors;
  int nobound, nopolywritten, nonodewritten, noelewritten, noiterationnum;
  int noholes;
  int docheck;
  int quiet, verbose;

/* Variables for file names.                                                 */

  char innodefilename[FILENAMESIZE];
  char inelefilename[FILENAMESIZE];
  char inpolyfilename[FILENAMESIZE];
  char areafilename[FILENAMESIZE];
  char outnodefilename[FILENAMESIZE];
  char outelefilename[FILENAMESIZE];
  char outpolyfilename[FILENAMESIZE];
  char edgefilename[FILENAMESIZE];
  char neighborfilename[FILENAMESIZE];
  char epsfilename[FILENAMESIZE];

} behavior_legacy;

void syntax();
void info();
void parsecommandline_legacy(int argc, char **argv, behavior_legacy *b);
void check_behavior(behavior *b, behavior_legacy *lb);

void print_statistics(context *ctx, int verbose);

void reset_triangleio(triangleio *io);
void free_triangleio(triangleio *io);
void check_mesh(context *ctx);

void check(int code, const char *site);

void finishfile(FILE *file, int argc, char **argv);

int write_nodes(context *ctx, char *filename, int argc, char **argv);
int write_elements(context *ctx, char *filename, int argc, char **argv);
int write_poly(context *ctx, char *filename, triangleio *in, int argc, char **argv);
int write_edges(context *ctx, char *filename, int argc, char **argv);
int write_neighbors(context *ctx, char *filename, int argc, char **argv);
int write_postscript(context *ctx, char *filename, int argc, char **argv);

#endif /* UTIL_H */