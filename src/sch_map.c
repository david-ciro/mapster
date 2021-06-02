#include "sch_map.h"
#include <stdbool.h>

// MAIN INTERNAL STRUCTURE

typedef struct
{
    size_t dim; // dimension of dynamical system
    int (*fw)(const gsl_vector* x0,
              gsl_vector*       x1,
              void*             params); // forward map
    int (*bw)(const gsl_vector* x0,
              gsl_vector*       x1,
              void*             params); // backwards map
    int (*jac)(const gsl_vector* x0,
               gsl_matrix*       J,
               void*             params); // jacobian of map
    void* params;             // map parameters

    bool   ex_jac; // true for exact jacobian
    double dx;     // finite element for apx jacobian

    gsl_vector* xa; // auxiliary state vector
    gsl_vector* xb; // auxiliary state vector
    gsl_vector* xc; // auxiliary state vector
    gsl_matrix* Ja; // auxiliary state matrix
    gsl_matrix* Jb; // auxiliary state matrix
    gsl_matrix* Jc; // auxiliary state matrix

} sch_map_t;

// INTERNAL METHODS DEFINITIONS

int // approximated jacobian by forward differences
map_jac_apx(sch_map_t* m, const gsl_vector* x0, gsl_matrix* J);

// API IMPLEMENTATION

sch_map // allocate map and initialize parameters
sch_map_create(size_t dim,
               int (*fw)(const gsl_vector* x0,
                         gsl_vector*       x1,
                         void*             params), // forward map
               int (*bw)(const gsl_vector* x0,
                         gsl_vector*       x1,
                         void*             params), // backwards map
               void* params             // map parameters
)
{
    sch_map_t* map = (sch_map_t*)malloc(sizeof(sch_map_t));
    map->dim = dim;
    map->fw = fw;
    map->bw = bw;
    map->jac = NULL;
    map->params = params;

    map->ex_jac = false;
    map->dx = 1e-6;
    map->xa = gsl_vector_alloc(dim);
    map->xb = gsl_vector_alloc(dim);
    map->xc = gsl_vector_alloc(dim);
    map->Ja = gsl_matrix_alloc(dim, dim);
    map->Jb = gsl_matrix_alloc(dim, dim);
    map->Jc = gsl_matrix_alloc(dim, dim);
    return (sch_map)map;
}

void
sch_map_destroy(sch_map map_h)
{
    sch_map_t* map = (sch_map_t*)map_h;
    gsl_vector_free(map->xa);
    gsl_vector_free(map->xb);
    gsl_vector_free(map->xc);
    gsl_matrix_free(map->Ja);
    gsl_matrix_free(map->Jb);
    gsl_matrix_free(map->Jc);
    free(map);
}

void // define jacobian function for the forward map
sch_map_set_jacobian(sch_map map_h,
                     int (*jac)(const gsl_vector* x0,
                                gsl_matrix*       J,
                                void*             params))
{
    sch_map_t* map = (sch_map_t*)map_h;
    map->jac = jac;
}

size_t
sch_map_get_dim(sch_map map_h)
{
    sch_map_t* map = (sch_map_t*)map_h;
    return map->dim;
}

int // map forward 'order' times
sch_map_fw(sch_map           map_h,
           size_t            order,
           const gsl_vector* x0,
           gsl_vector*       x1)
{
    sch_map_t* map = (sch_map_t*)map_h;
    gsl_vector_memcpy(map->xa, x0);
    int status = 0;
    for (size_t i = 0; i < order && status == 0; i++) {
        status = map->fw(map->xa, map->xb, map->params);
        // xa <-> xb
        gsl_vector* tmp = map->xa;
        map->xa = map->xb;
        map->xb = tmp;
    }
    gsl_vector_memcpy(x1, map->xa);
    return status;
}

int // map backwards 'order' times
sch_map_bw(sch_map           map_h,
           size_t            order,
           const gsl_vector* x0,
           gsl_vector*       x1)
{
    sch_map_t* map = (sch_map_t*)map_h;
    gsl_vector_memcpy(map->xa, x0);
    int status = 0;
    for (size_t i = 0; i < order && status == 0; i++) {
        status = map->bw(map->xa, x1, map->params);
        // xa <-> xb
        gsl_vector* tmp = map->xa;
        map->xa = map->xb;
        map->xb = tmp;
    }
    gsl_vector_memcpy(x1, map->xa);
    return status;
}

int // evaluate map jacobian of 'order'
sch_map_jac(sch_map           map_h,
            size_t            order,
            const gsl_vector* x0,
            gsl_matrix*       J)
{
    sch_map_t* map = (sch_map_t*)map_h;
    gsl_vector_memcpy(map->xa, x0);   // xa = x0
    gsl_matrix_set_identity(map->Jb); // Jb = Id

    int status = 0;
    for (size_t i = 0; i < order && status == 0; i++) {
        // Jn(x) = J(Tn-1(x))J(Tn-2(x))...J(x)
        if (map->ex_jac) // exact jacobian
            status = map->jac(map->xa, map->Ja, map->params);
        else // approximated jacobian
            status = map_jac_apx(map, map->xa, map->Ja);
        // xb = T(xa)
        status += map->fw(map->xa, map->xb, map->params);
        // xa <-> xb
        gsl_vector* tmp = map->xa;
        map->xa = map->xb;
        map->xb = tmp;
        // Jc = Ja x Jb
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, map->Ja,
                       map->Jb, 0.0, map->Jc);
        // Jb = Jc
        gsl_matrix_memcpy(map->Jb, map->Jc);
    }
    return status;
}

// INTERNAL METHODS IMPLEMENTATION

int
sch_map_jac_apx(sch_map_t* map, const gsl_vector* x0, gsl_matrix* J)
{
    gsl_vector_memcpy(map->xa, x0);
    int status = 0;
    // T(x0)
    status = map->fw(map->xa, map->xb, map->params);
    for (size_t j = 0; j < map->dim && status == 0; j++) {
        // T(x0 + dxj)
        double xj = gsl_vector_get(map->xa, j);
        gsl_vector_set(map->xa, j, xj + map->dx);
        status = map->fw(map->xa, map->xc, map->params);
        // [Ti(x0 + dxj) - Ti(x0)]/dxj
        for (size_t i = 0; i < map->dim && status == 0; i++) {
            double T1 = gsl_vector_get(map->xc, i);
            double T0 = gsl_vector_get(map->xb, i);
            gsl_matrix_set(J, i, j, (T1 - T0) / map->dx);
        }
        // restore x0
        gsl_vector_set(map->xa, j, xj - map->dx);
    }
    return status;
}