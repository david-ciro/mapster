#include "orbit.h"
#include <math.h>

typedef struct
{
    map_h        m;
    size_t       dim;
    size_t       order;
    size_t       steps;
    bool         forward;
    gsl_vector** x;
} orbit;

orbit_h
orb_create(map_h             m,
           size_t            order,
           bool              forward,
           size_t            steps,
           const gsl_vector* x0)
{
    orbit* orb = (orbit*)malloc(sizeof(orbit));
    orb->m = m;
    orb->dim = map_get_dim(m);
    orb->order = order;
    orb->steps = steps;
    orb->forward = forward;
    // memory allocation
    orb->x = (gsl_vector**)malloc((steps + 1) * sizeof(gsl_vector*));
    for (size_t i = 0; i <= steps; i++)
        orb->x[i] = gsl_vector_alloc(orb->dim);
    // orbit generation
    gsl_vector_memcpy(orb->x[0], x0);
    for (size_t i = 1; i <= steps; i++) {
        if (orb->forward)
            map_fw(orb->m, orb->order, orb->x[i - 1], orb->x[i]);
        else
            map_bw(orb->m, orb->order, orb->x[i - 1], orb->x[i]);
    }
    return (orbit_h)orb;
}

void
orb_destroy(orbit_h orbh)
{
    orbit* orb = (orbit*)orbh;
    for (size_t i = 0; i <= orb->steps; i++)
        gsl_vector_free(orb->x[i]);
    free(orb);
}

double
orb_get_comp(orbit_h orbh, size_t point, size_t comp)
{
    orbit* orb = (orbit*)orbh;
    if (point <= orb->steps && comp < orb->dim) {
        return gsl_vector_get(orb->x[point], comp);
    } else {
        fprintf(stderr, "orb_get_comp: invalid argument");
        return NAN;
    }
}

void
orb_save(orbit_h orbh, char* filename)
{
    orbit* orb = (orbit*)orbh;
    FILE*  f = fopen(filename, "w");
    for (size_t i = 0; i <= orb->steps; i++) {
        for (size_t j = 0; j < orb->dim; j++) {
            fprintf(f, "%.6e  ", gsl_vector_get(orb->x[i], j));
        }
        fprintf(f, "\n");
    }
    fclose(f);
}