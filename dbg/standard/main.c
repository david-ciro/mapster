#include "sch_map.h"
#include "sch_orbit.h"
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdbool.h>

typedef struct
{
    double k;
} std_params;

int
std_fw(const gsl_vector* x0, gsl_vector* x1, void* params)
{
    std_params* stdp = (std_params*)params;
    double      th0 = gsl_vector_get(x0, 0);
    double      p0 = gsl_vector_get(x0, 1);
    double      p1 = p0 + (stdp->k) * sin(th0);
    double      th1 = th0 + p1;
    p1 = fmod(p1, 2 * M_PI);
    th1 = fmod(th1, 2 * M_PI);
    gsl_vector_set(x1, 0, th1);
    gsl_vector_set(x1, 1, p1);
    return 0;
}

int
std_bw(const gsl_vector* x0, gsl_vector* x1, void* params)
{
    std_params* stdp = (std_params*)params;
    double      th1 = gsl_vector_get(x0, 0);
    double      p1 = gsl_vector_get(x0, 1);
    double      th0 = th1 - p1;
    double      p0 = p1 - (stdp->k) * sin(th0);
    p0 = fmod(p0, 2 * M_PI);
    th0 = fmod(th0, 2 * M_PI);
    gsl_vector_set(x1, 0, th0);
    gsl_vector_set(x1, 1, p0);
    return 0;
}

int
std_jac(const gsl_vector* x0, gsl_matrix* J, void* params)
{
    std_params* stdp = (std_params*)params;
    double      th0 = gsl_vector_get(x0, 0);
    double      dp1_dth0 = (stdp->k) * cos(th0);
    double      dp1_dp0 = 1;
    double      dth1_dth0 = 1 + (stdp->k) * cos(th0);
    double      dth1_dp0 = 1;
    gsl_matrix_set(J, 0, 0, dth1_dth0);
    gsl_matrix_set(J, 0, 1, dth1_dp0);
    gsl_matrix_set(J, 1, 0, dp1_dth0);
    gsl_matrix_set(J, 1, 1, dp1_dp0);
    return 0;
}

int
main()
{
    size_t     dim = 2;
    std_params stdp = { 0.4 };
    sch_map    map = sch_map_create(dim, std_fw, std_bw, &stdp);

    // initial condition
    size_t      n_orbits = 100;
    FILE*       file = fopen("orbits.dat", "w");
    gsl_vector* x0 = gsl_vector_alloc(dim);
    for (size_t i = 0; i < n_orbits; i++) {
        // random initial condition
        for (size_t i = 0; i < dim; i++)
            gsl_vector_set(x0, i, 2 * M_PI * drand48());
        // orbit properties
        size_t order = 1;
        bool   forward = true;
        size_t steps = 10000;
        // create orbit
        sch_orbit orbit =
          sch_orbit_create(map, order, forward, steps, x0);
        // save orbit
        sch_orbit_save(orbit, file);
        // destroy orbit
        sch_orbit_destroy(orbit);
    }

    // free memory
    gsl_vector_free(x0);
    sch_map_destroy(map);
    return 0;
}