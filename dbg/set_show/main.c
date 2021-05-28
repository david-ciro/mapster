#include "map.h"
#include <stdio.h>
#include<math.h>

typedef struct
{
    double k;
} ch_params;

int
ch_fw(double* x0, double* x1, void* params)
{
    ch_params* par = (ch_params*)params;
    x1[1] = x0[1] + (par->k) * sin(x0[0]);
    x1[0] = x0[0] + x1[1];
    return 0;
}

int
ch_bw(double* x0, double* x1, void* params)
{
    return 0;
}

int
main()
{
    map_h m = map_create();
    map_set_forward(m, ch_fw);
    map_set_backward(m, ch_bw);

    ch_params par = {1.46};
    double x0[] = {0.1, 0.3};
    double x1[2];

    map_set_params(m, &par);

    map_fw(m, x0, x1);

    printf("x1 = %lf  %lf\n", x1[0], x1[1]);

    map_destroy(m);
}