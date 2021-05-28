#include "map.h"

typedef struct
{
    int (*fw)(double* x0, double* x1, void* params);
    int (*bw)(double* x0, double* x1, void* params);
    void* params;
} map;

map_h
map_create()
{
    return (map_h)malloc(sizeof(map));
}

void
map_destroy(map_h m)
{
    free((map*)m);
}

void
map_set_forward(map_h mh, int (*fw)(double* x0, double* x1, void* params)){
    map* m = (map*)mh;
    m->fw = fw;
}

void
map_set_backward(map_h mh, int (*bw)(double* x0, double* x1, void* params)){
    map* m = (map*)mh;
    m->bw = bw;
}

void
map_set_params(map_h mh, void* params){
    map* m = (map*)mh;
    m->params = params;
}

void
map_fw(map_h mh, double* x0, double* x1){
    map* m = (map*)mh;
    m->fw(x0, x1, m->params);
}

void
map_bw(map_h mh, double* x0, double* x1){
    map* m = (map*)mh;
    m->bw(x0, x1, m->params);
}
