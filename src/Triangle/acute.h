
// ACUTE MEMORY POOL
struct acutepool {
    int size;
    // getWedgeIntersection (fixed size)
    REAL *initialpoly;
    // getWedgeIntersection (dynamic size)
    REAL *petalx;
    REAL *petaly;
    REAL *petalr;
    REAL *wedges;
    // doSmoothing (fixed size [500])
    REAL *points_p;
    REAL *points_q;
    REAL *points_r;
};

void acutepool_init(int n, behavior *b, struct acutepool *p) {
    p->size = n;
    
    p->initialpoly = (REAL *)malloc(sizeof(REAL)* 500);
    p->petalx = (REAL *)malloc(sizeof(REAL)* 2 * n);
    p->petaly = (REAL *)malloc(sizeof(REAL)* 2 * n);
    p->petalr = (REAL *)malloc(sizeof(REAL)* 2 * n);
    if(b->maxangle == 0.00000){
       p->wedges = (REAL *)malloc(sizeof(REAL)* 2 * n * 16 + 36);
    }else{
       p->wedges = (REAL *)malloc(sizeof(REAL)* 2 * n * 20 + 40);
    }

    p->points_p = (REAL *)malloc(sizeof(REAL)* 500);
    p->points_q = (REAL *)malloc(sizeof(REAL)* 500);
    p->points_r = (REAL *)malloc(sizeof(REAL)* 500);
}

void acutepool_resize(int n, behavior *b, struct acutepool *p) {
    if (p->size < n) {
        p->size = n;

        // Free old memory
        free(p->petalx);
        free(p->petaly);
        free(p->petalr);
        free(p->wedges);

        // Allocate new memory
        p->petalx = (REAL *)malloc(sizeof(REAL)* 2 * n);
        p->petaly = (REAL *)malloc(sizeof(REAL)* 2 * n);
        p->petalr = (REAL *)malloc(sizeof(REAL)* 2 * n);
        if(b->maxangle == 0.00000){
           p->wedges = (REAL *)malloc(sizeof(REAL)* 2 * n * 16 + 36);
        }else{
           p->wedges = (REAL *)malloc(sizeof(REAL)* 2 * n * 20 + 40);
        }
    }
}

void acutepool_deinit(struct acutepool *p) {
  free(p->initialpoly);
  free(p->petalx);
  free(p->petaly);
  free(p->petalr);
  free(p->wedges);

  free(p->points_p);
  free(p->points_q);
  free(p->points_r);
}
// END ACUTE MEMORY POOL