#ifndef MAP_H
#define MAP_H

#include <stdio.h>
#include <stdlib.h>

typedef struct map* map_h;

map_h map_create();
void map_destroy(map_h m);

void map_set_forward(map_h m, int(*fw)(double* x0, double* x1, void* params));
void map_set_backward(map_h m, int(*bw)(double* x0, double* x1, void* params));
void map_set_params(map_h m, void* params);

void map_fw(map_h m, double* x0, double* x1);
void map_bw(map_h m, double* x0, double* x1);

#endif // MAP_H