#ifndef SCH_sch_map
#define SCH_sch_map

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct sch_map_t* sch_map;

typedef enum
{
    map_dx, // approx. jacobian finite dif.
} map_params;

sch_map // allocate map and set user parameters
sch_map_create(size_t dim,
               int (*fw)(const gsl_vector* x0,
                         gsl_vector*       x1,
                         void*             params), // forward map
               int (*bw)(const gsl_vector* x0,
                         gsl_vector*       x1,
                         void*             params), // backwards map
               void* params             // map parameters
);

void // free map memory
sch_map_destroy(sch_map m);

void // define jacobian function for the forward map
sch_map_set_jacobian(sch_map m,
                     int (*jac)(const gsl_vector* x0,
                                gsl_matrix*       J,
                                void*             params));

size_t // get system dimension
  sch_map_get_dim(sch_map);

int // map forward 'order' times
sch_map_fw(sch_map           m,
           size_t            order,
           const gsl_vector* x0,
           gsl_vector*       x1);

int // map backwards 'order' times
sch_map_bw(sch_map           m,
           size_t            order,
           const gsl_vector* x0,
           gsl_vector*       x1);

int // evaluate map jacobian of 'order'
sch_map_jac(sch_map           m,
            size_t            order,
            const gsl_vector* x0,
            gsl_matrix*       J);

int // search fixed point of 'order'
sch_map_fixed_point(sch_map           m,
                    size_t            order,
                    const gsl_vector* x0,
                    gsl_vector*       x_fixed);

#endif // sch_map