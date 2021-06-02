#include "sch_orbit.h"
#include <math.h>

typedef struct
{
    sch_map      map;
    size_t       dim;
    size_t       order;
    size_t       steps;
    bool         forward;
    gsl_vector** x;
} sch_orbit_t;

sch_orbit
sch_orbit_create(sch_map           map,
                 size_t            order,
                 bool              forward,
                 size_t            steps,
                 const gsl_vector* x0)
{
    sch_orbit_t* orbit = (sch_orbit_t*)malloc(sizeof(sch_orbit_t));
    orbit->map = map;
    orbit->dim = sch_map_get_dim(map);
    orbit->order = order;
    orbit->steps = steps;
    orbit->forward = forward;
    // memory allocation
    orbit->x =
      (gsl_vector**)malloc((steps + 1) * sizeof(gsl_vector*));
    for (size_t i = 0; i <= steps; i++)
        orbit->x[i] = gsl_vector_alloc(orbit->dim);
    // orbit generation
    gsl_vector_memcpy(orbit->x[0], x0);
    for (size_t i = 1; i <= steps; i++) {
        if (orbit->forward)
            sch_map_fw(orbit->map, orbit->order, orbit->x[i - 1],
                       orbit->x[i]);
        else
            sch_map_bw(orbit->map, orbit->order, orbit->x[i - 1],
                       orbit->x[i]);
    }
    return (sch_orbit)orbit;
}

void
sch_orbit_destroy(sch_orbit orbit_h)
{
    sch_orbit_t* orbit = (sch_orbit_t*)orbit_h;
    for (size_t i = 0; i <= orbit->steps; i++)
        gsl_vector_free(orbit->x[i]);
    free(orbit);
}

double
sch_orbit_get_comp(sch_orbit orbit_h, size_t point, size_t comp)
{
    sch_orbit_t* orbit = (sch_orbit_t*)orbit_h;
    if (point <= orbit->steps && comp < orbit->dim) {
        return gsl_vector_get(orbit->x[point], comp);
    } else {
        fprintf(stderr, "orb_get_comp: invalid argument");
        return NAN;
    }
}

void
sch_orbit_save(sch_orbit orbit_h, FILE* f)
{
    sch_orbit_t* orbit = (sch_orbit_t*)orbit_h;
    for (size_t i = 0; i <= orbit->steps; i++) {
        for (size_t j = 0; j < orbit->dim; j++) {
            fprintf(f, "%.6e  ", gsl_vector_get(orbit->x[i], j));
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}