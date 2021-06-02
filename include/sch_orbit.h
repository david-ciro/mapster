#ifndef ORBIT_H
#define ORBIT_H

#include "sch_map.h"
#include <stdbool.h>

typedef struct sch_orbit_t* sch_orbit;

// create orbit from map^order, use steps < 0 for backwards orbit
sch_orbit
sch_orbit_create(sch_map           map,
               size_t            order,
               bool              forward,
               size_t            steps,
               const gsl_vector* x0);

void
sch_orbit_destroy(sch_orbit orb);

// get component 'comp' of a point in the orbit
double
sch_orbit_get_comp(sch_orbit orb, size_t point, size_t comp);

void
sch_orbit_save(sch_orbit orb, FILE* f);

#endif // ORBIT_H
