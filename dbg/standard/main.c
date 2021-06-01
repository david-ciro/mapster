#include "map.h"
#include "orbit.h"
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
    map_h      m = map_create(dim, std_fw, std_bw, &stdp);

    // initial condition
    double      x[] = { 0.5, 0.5 };
    gsl_vector* x0 = gsl_vector_alloc(dim);
    for (size_t i = 0; i < dim; i++)
        gsl_vector_set(x0, i, x[i]);
    // orbit properties
    size_t order = 1;
    bool   forward = true;
    size_t steps = 10000;
    // create orbit
    orbit_h orb = orb_create(m, order, forward, steps, x0);
    // save orbit
    orb_save(orb, "orbit.dat");

    // free memory
    orb_destroy(orb);
    map_destroy(m);
    return 0;
}