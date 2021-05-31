#ifndef MAP_H
#define MAP_H

#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct map* map_h;

typedef enum
{
    map_dx, // approx. jacobian finite dif.
} map_params;

map_h // allocate map and set user parameters
map_create(size_t dim,
           int (*fw)(const gsl_vector* x0,
                     gsl_vector*       x1,
                     void*             params), // forward map
           int (*bw)(const gsl_vector* x0,
                     gsl_vector*       x1,
                     void*             params), // backwards map
           void* params             // map parameters
);

void // free map memory
map_destroy(map_h m);

void // define jacobian function for the forward map
map_set_jacobian(map_h m,
                 int (*jac)(const gsl_vector* x0,
                            gsl_matrix*       J,
                            void*             params));

int // map forward 'order' times
map_fw(map_h m, size_t order, const gsl_vector* x0, gsl_vector* x1);

int // map backwards 'order' times
map_bw(map_h m, size_t order, const gsl_vector* x0, gsl_vector* x1);

int // evaluate map jacobian of 'order'
map_jac(map_h m, size_t order, const gsl_vector* x0, gsl_matrix* J);

int // search fixed point of 'order'
map_fixed_point(map_h             m,
                size_t            order,
                const gsl_vector* x0,
                gsl_vector*       x_fixed);

#endif // MAP_H