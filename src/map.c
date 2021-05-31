#include "map.h"
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

} map;

// INTERNAL METHODS DEFINITIONS

int // approximated jacobian by forward differences
map_jac_apx(map* m, const gsl_vector* x0, gsl_matrix* J);

// API IMPLEMENTATION

map_h // allocate map and initialize parameters
map_create(size_t dim,
           int (*fw)(const gsl_vector* x0,
                     gsl_vector*       x1,
                     void*             params), // forward map
           int (*bw)(const gsl_vector* x0,
                     gsl_vector*       x1,
                     void*             params), // backwards map
           void* params             // map parameters
)
{
    map* m = (map*)malloc(sizeof(map));
    m->dim = dim;
    m->fw = fw;
    m->bw = bw;
    m->jac = NULL;
    m->params = params;

    m->ex_jac = false;
    m->dx = 1e-6;
    m->xa = gsl_vector_alloc(dim);
    m->xb = gsl_vector_alloc(dim);
    m->xc = gsl_vector_alloc(dim);
    m->Ja = gsl_matrix_alloc(dim, dim);
    m->Jb = gsl_matrix_alloc(dim, dim);
    m->Jc = gsl_matrix_alloc(dim, dim);
    return (map_h)m;
}

void
map_destroy(map_h mh)
{
    map* m = (map*)mh;
    gsl_vector_free(m->xa);
    gsl_vector_free(m->xb);
    gsl_vector_free(m->xc);
    gsl_matrix_free(m->Ja);
    gsl_matrix_free(m->Jb);
    gsl_matrix_free(m->Jc);
    free(m);
}

void // define jacobian function for the forward map
map_set_jacobian(map_h mh,
                 int (*jac)(const gsl_vector* x0,
                            gsl_matrix*       J,
                            void*             params))
{
    map* m = (map*)mh;
    m->jac = jac;
}

int // map forward 'order' times
map_fw(map_h mh, size_t order, const gsl_vector* x0, gsl_vector* x1)
{
    map* m = (map*)mh;
    gsl_vector_memcpy(m->xa, x0);
    int status = 0;
    for (size_t i = 0; i < order && status == 0; i++) {
        status = m->fw(m->xa, m->xb, m->params);
        // xa <-> xb
        gsl_vector* tmp = m->xa;
        m->xa = m->xb;
        m->xb = tmp;
    }
    gsl_vector_memcpy(x1, m->xa);
    return status;
}

int // map backwards 'order' times
map_bw(map_h mh, size_t order, const gsl_vector* x0, gsl_vector* x1)
{
    map* m = (map*)mh;
    gsl_vector_memcpy(m->xa, x0);
    int status = 0;
    for (size_t i = 0; i < order && status == 0; i++) {
        status = m->bw(m->xa, x1, m->params);
        // xa <-> xb
        gsl_vector* tmp = m->xa;
        m->xa = m->xb;
        m->xb = tmp;
    }
    gsl_vector_memcpy(x1, m->xa);
    return status;
}

int // evaluate map jacobian of 'order'
map_jac(map_h mh, size_t order, const gsl_vector* x0, gsl_matrix* J)
{
    map* m = (map*)mh;
    gsl_vector_memcpy(m->xa, x0); // xa = x0
    gsl_matrix_set_identity(m->Jb); // Jb = Id
    
    int status = 0;
    for (size_t i = 0; i < order && status == 0; i++) {
        // Jn(x) = J(Tn-1(x))J(Tn-2(x))...J(x)
        if (m->ex_jac) // exact jacobian
            status = m->jac(m->xa, m->Ja, m->params);
        else // approximated jacobian
            status = map_jac_apx(m, m->xa, m->Ja);
        // xb = T(xa)
        status += m->fw(m->xa, m->xb, m->params);
        // xa <-> xb
        gsl_vector* tmp = m->xa;
        m->xa = m->xb;
        m->xb = tmp;
        // Jc = Ja x Jb
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m->Ja, m->Jb,
                       0.0, m->Jc);
        // Jb = Jc
        gsl_matrix_memcpy(m->Jb, m->Jc);
    }
    return status;
}

// INTERNAL METHODS IMPLEMENTATION

int
map_jac_apx(map* m, const gsl_vector* x0, gsl_matrix* J)
{
    gsl_vector_memcpy(m->xa, x0);
    int status = 0;
    // T(x0)
    status = m->fw(m->xa, m->xb, m->params);
    for (size_t j = 0; j < m->dim && status == 0; j++) {
        // T(x0 + dxj)
        double xj = gsl_vector_get(m->xa, j);
        gsl_vector_set(m->xa, j, xj + m->dx);
        status = m->fw(m->xa, m->xc, m->params);
        // [Ti(x0 + dxj) - Ti(x0)]/dxj
        for (size_t i = 0; i < m->dim && status == 0; i++) {
            double T1 = gsl_vector_get(m->xc, i);
            double T0 = gsl_vector_get(m->xb, i);
            gsl_matrix_set(J, i, j, (T1 - T0) / m->dx);
        }
        // restore x0
        gsl_vector_set(m->xa, j, xj - m->dx);
    }
    return status;
}