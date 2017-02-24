#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maptools.h"

int main(void)
{
    double B1, L1, S, A1, B2, L2, A2;
    double _A12, _S = 0, _A21;

    B1 = 45.0 * M_PI / 180;
    L1 = 0.0;
    S  = 19500000.0;
    A1 = 265 * M_PI / 180;

    mt_point_t pnt1, pnt2;
    pnt1.b = B1;
    pnt1.l = L1;

    // direct
    pnt2 = mt_point_move(&pnt1, S, A1, MtGeoidKras);
    printf("B2=%f\tL2=%f\n", pnt2.b, pnt2.l);
    printf("B2=%f\tL2=%f\n", pnt2.b * 180.0 / M_PI, pnt2.l * 180.0 / M_PI);

    // reverse
    double new_s = mt_distance(&pnt1, &pnt2, MtGeoidKras);
    double azim  = mt_azimuth(&pnt1, &pnt2, MtGeoidKras);
    printf("S=%f,\tA1=%f\n", new_s, azim);

    exit(0);
}
