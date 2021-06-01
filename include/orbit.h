#ifndef ORBIT_H
#define ORBIT_H

#include "map.h"
#include <stdbool.h>

typedef struct orbit* orbit_h;

// create orbit from map^order, use steps < 0 for backwards orbit
orbit_h
orb_create(map_h             m,
           size_t            order,
           bool              forward,
           size_t            steps,
           const gsl_vector* x0);

void
orb_destroy(orbit_h orb);

// get component 'comp' of a point in the orbit
double
orb_get_comp(orbit_h orb, size_t point, size_t comp);

void
orb_save(orbit_h orb, char* filename);

#endif // ORBIT_H
