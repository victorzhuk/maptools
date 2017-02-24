#include "maptools.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 * Определение закрытых функций
 */

static double radian_normalize(double val)
{
    if (val == 0.0) return 0.0;
    
    double new_val = 0.0;
    
    int cnt = fabs(val / (2 * M_PI));

    if (cnt > 0) {
        if (val > 0) {
            new_val = val - 2 * M_PI * cnt;
        } else {
            new_val = val + 2 * M_PI * (++cnt);
        }
    } else {
        if (val > 0) {
            new_val = val;
        } else {
            new_val = val + 2 * M_PI;
        }
    }

    return new_val;
}


/*
 * Определение экспортируемых функций
 */

mt_geoid_t mt_geoid_create(mt_geoidtypes_t geoid_type)
{
    mt_geoid_t geoid;

    switch (geoid_type) {
        case MtGeoidKras:
            geoid.a        = 6378245.0;
            geoid.b        = 6356863.0188;
            geoid.alpha    = 0.0033523299;
            geoid.e_sq     = 0.006693421623;
            geoid.e2_sq    = 0.006738525415;
            break;
        default:
            geoid.a        = 0.0;
            geoid.b        = 0.0;
            geoid.alpha    = 0.0;
            geoid.e_sq     = 0.0;
            geoid.e2_sq    = 0.0;
            break;
    }

    return geoid;
}

mt_point_t mt_point_move(const mt_point_t *pnt, double dist, double azim, mt_geoidtypes_t geoid_type)
{
    // получаем геоид
    mt_geoid_t geoid = mt_geoid_create(geoid_type);

    // нормализация азимута
    azim = radian_normalize(azim);
    
    // 1. Вычисление вспомогательных функций и коэффициентов
    // Вычисление приведенной широты начальной точки
    double var_W1     = sqrt(1 - geoid.e_sq * pow(sin(pnt->b), 2.0));
    // Вычисление вспомогательных функций
    double sin_u1     = sin(pnt->b) * sqrt(1 - geoid.e_sq) / var_W1;
    double cos_u1     = cos(pnt->b) / var_W1;
    double sin_A0     = cos_u1 * sin(azim);
    double ctg_sig1   = cos_u1 * cos(azim) / sin_u1;
    double sin_2sig1  = 2 * ctg_sig1 / (pow(ctg_sig1, 2.0) + 1);
    double cos_2sig1  = (pow(ctg_sig1, 2.0) - 1) / (pow(ctg_sig1, 2.0) + 1);
    // Вычисление вспомогательных коэффициентов
    double cos_sq_A0  = 1 - pow(sin_A0, 2.0);
    // хардкодные коэффициенты
    double var_A      = (10708.949 - 13.474 * cos_sq_A0) * cos_sq_A0 + geoid.b;
    double var_B      = (5354.469 - 8.978 * cos_sq_A0) * cos_sq_A0;
    double var_C      = (2.238 * cos_sq_A0) * cos_sq_A0 + 0.006;
    double var_alp    = 691.46768 - (0.58143 - 0.00144 * cos_sq_A0) * cos_sq_A0;
    double var_bet    = (0.2907 - 0.0010 * cos_sq_A0) * cos_sq_A0;
    
    // 2. Вычисление сферического расстояния
    double var_sig0   = (dist - (var_B + var_C * cos_2sig1) * sin_2sig1) / var_A;
    double sin_2sig01 = sin_2sig1 * cos(2 * var_sig0) + cos_2sig1 * sin(2 * var_sig0);
    double cos_2sig01 = cos_2sig1 * cos(2 * var_sig0) - sin_2sig1 * sin(2 * var_sig0);
    double var_sig    = var_sig0 + (var_B + 5 * var_C * cos_2sig01) * sin_2sig01 / var_A;
    
    // 3. Вычисление поправки в разность долгот
    double var_dlt    = (var_alp * var_sig + var_bet * (sin_2sig01 - sin_2sig1)) * sin_A0;
    
    // 4. Вычисление геодезических координат в конечной точке
    double sin_u2     = sin_u1 * cos(var_sig) + cos_u1 * cos(azim) * sin(var_sig);
    double var_lam    = atan(sin(azim) * sin(var_sig) / (cos_u1 * cos(var_sig) - sin_u1 * sin(var_sig) * cos(azim)));
    if (sin(azim) > 0) {
        if (tan(var_lam) > 0) {
            var_lam = fabs(var_lam);
        } else {
            var_lam = M_PI - fabs(var_lam);
        }
    } else {
        if (tan(var_lam) > 0) {
            var_lam = fabs(var_lam) - M_PI;
        } else {
            var_lam = -(fabs(var_lam));
        }
    }
    
    double b = atan(sin_u2 / (sqrt(1 - geoid.e_sq) * sqrt(1 - pow(sin_u2, 2.0))));
    double l = radian_normalize(pnt->l + var_lam - var_dlt);

    mt_point_t new_pnt;
    new_pnt.b = b;
    new_pnt.l = l;

    return new_pnt;
}

double mt_distance(const mt_point_t *src_pnt, const mt_point_t *dst_pnt, mt_geoidtypes_t geoid_type)
{
    // Расчет расстояния между двумя точками

    // получаем геоид
    mt_geoid_t geoid = mt_geoid_create(geoid_type);
    
    // 1. Подготовительные вычисления
    double var_W1 = sqrt(1 - geoid.e_sq * pow(sin(src_pnt->b), 2.0));
    double var_W2 = sqrt(1 - geoid.e_sq * pow(sin(dst_pnt->b), 2.0));
    double sin_u1 = sin(src_pnt->b) * sqrt(1 - geoid.e_sq) / var_W1;
    double sin_u2 = sin(dst_pnt->b) * sqrt(1 - geoid.e_sq) / var_W2;
    double cos_u1 = cos(src_pnt->b) / var_W1;
    double cos_u2 = cos(dst_pnt->b) / var_W2;
    double var_l = dst_pnt->l - src_pnt->l;
    if (var_l > M_PI) {
        var_l -= 2 * M_PI;
    }
    if (var_l < -M_PI) {
        var_l += 2 * M_PI;
    }
    double var_a1 = sin_u1 * sin_u2;
    double var_a2 = cos_u1 * cos_u2;
    double var_b1 = cos_u1 * sin_u2;
    double var_b2 = sin_u1 * cos_u2;
    
    // 2. Последовательные приближения
    double  var_lam, 
            var_A1, 
            var_sig, 
            sin_A0, 
            cos_sq_A0, 
            sin_sig, 
            cos_sig, 
            var_x, 
            var_dlt = 0.0;
    int     i, break_next = 0;
    for (i=1;;i++) {
        var_lam = var_l + var_dlt;
        double var_p = cos_u2 * sin(var_lam);
        double var_q = var_b1 - var_b2 * cos(var_lam);
        var_A1 = atan(var_p / var_q);
        if (var_p > 0) {
            if (var_q > 0) {
                var_A1 = fabs(var_A1);
            } else {
                var_A1 = M_PI - fabs(var_A1);
            }
        } else {
            if (var_q > 0) {
                var_A1 = 2 * M_PI - fabs(var_A1);
            } else {
                var_A1 = M_PI + fabs(var_A1);
            }
        }
        double sin_A1 = sin(var_A1);
        double cos_A1 = cos(var_A1);
        sin_sig = var_p * sin(var_A1) + var_q * cos(var_A1);
        cos_sig = var_a1 + var_a2 * cos(var_lam);
        var_sig = atan(sin_sig / cos_sig);
        if (cos_sig > 0) {
            var_sig = fabs(var_sig);
        } else {
            var_sig = M_PI - fabs(var_sig);
        }
        sin_A0 = cos_u1 * sin(var_A1);
        cos_sq_A0 = 1 - pow(sin_A0, 2.0);
        var_x = 2 * var_a1 - cos_sq_A0 * cos_sig;
        if (break_next > 0) {
            break;
        }
        double var_alp = (33523299 - (28189 - 70 * cos_sq_A0) * cos_sq_A0) * 1e-10;
        double var_2bet = (28189 - 94 * cos_sq_A0) * 1e-10;
        double var_2dlt = (var_alp * var_sig - var_2bet * var_x * sin_sig) * sin_A0;
        double var_diff = fabs(var_2dlt - var_dlt);
        var_dlt = var_2dlt;
        if (var_diff < 1e-10) {
            break_next = 1;
        }
    }
    
    // 3. Вычисление коэффициентов A, B' и C'
    double var_A = (10708.949 - 13.474 * cos_sq_A0) * cos_sq_A0 + geoid.b;
    double var_2B = 10708.938 - 17.956 * cos_sq_A0;
    double var_2C = 4.487;
    
    // 4. Вычисление расстояния
    double var_y = (pow(cos_sq_A0, 2.0) - 2 * pow(var_x, 2.0)) * cos_sig;
    double var_s = var_A * var_sig + (var_2B * var_x + var_2C * var_y) * sin_sig;
    
    return var_s;
}

double mt_azimuth(const mt_point_t *src_pnt, const mt_point_t *dst_pnt, mt_geoidtypes_t geoid_type)
{
    // Расчет азимута между двумя точками

    // получаем геоид
    mt_geoid_t geoid = mt_geoid_create(geoid_type);

    // 1. Подготовительные вычисления
    double var_W1 = sqrt(1 - geoid.e_sq * pow(sin(src_pnt->b), 2.0));
    double var_W2 = sqrt(1 - geoid.e_sq * pow(sin(dst_pnt->b), 2.0));
    double sin_u1 = sin(src_pnt->b) * sqrt(1 - geoid.e_sq) / var_W1;
    double sin_u2 = sin(dst_pnt->b) * sqrt(1 - geoid.e_sq) / var_W2;
    double cos_u1 = cos(src_pnt->b) / var_W1;
    double cos_u2 = cos(dst_pnt->b) / var_W2;
    double var_l = dst_pnt->l - src_pnt->l;
    if (var_l > M_PI) {
        var_l -= 2 * M_PI;
    }
    if (var_l < -M_PI) {
        var_l += 2 * M_PI;
    }
    double var_a1 = sin_u1 * sin_u2;
    double var_a2 = cos_u1 * cos_u2;
    double var_b1 = cos_u1 * sin_u2;
    double var_b2 = sin_u1 * cos_u2;
    
    // 2. Последовательные приближения
    double  var_lam, 
    var_A1, 
    var_sig, 
    sin_A0, 
    cos_sq_A0, 
    sin_sig, 
    cos_sig, 
    var_x, 
    var_dlt = 0.0;
    int     i, break_next = 0;
    for (i=1;;i++) {
        var_lam = var_l + var_dlt;
        double var_p = cos_u2 * sin(var_lam);
        double var_q = var_b1 - var_b2 * cos(var_lam);
        var_A1 = atan(var_p / var_q);
        if (var_p > 0) {
            if (var_q > 0) {
                var_A1 = fabs(var_A1);
            } else {
                var_A1 = M_PI - fabs(var_A1);
            }
        } else {
            if (var_q > 0) {
                var_A1 = 2 * M_PI - fabs(var_A1);
            } else {
                var_A1 = M_PI + fabs(var_A1);
            }
        }
        double sin_A1 = sin(var_A1);
        double cos_A1 = cos(var_A1);
        sin_sig = var_p * sin(var_A1) + var_q * cos(var_A1);
        cos_sig = var_a1 + var_a2 * cos(var_lam);
        var_sig = atan(sin_sig / cos_sig);
        if (cos_sig > 0) {
            var_sig = fabs(var_sig);
        } else {
            var_sig = M_PI - fabs(var_sig);
        }
        sin_A0 = cos_u1 * sin(var_A1);
        cos_sq_A0 = 1 - pow(sin_A0, 2.0);
        var_x = 2 * var_a1 - cos_sq_A0 * cos_sig;

        if (break_next > 0) {
            break;
        }

        double var_alp = (33523299 - (28189 - 70 * cos_sq_A0) * cos_sq_A0) * 1e-10;
        double var_2bet = (28189 - 94 * cos_sq_A0) * 1e-10;
        double var_2dlt = (var_alp * var_sig - var_2bet * var_x * sin_sig) * sin_A0;
        double var_diff = fabs(var_2dlt - var_dlt);
        var_dlt = var_2dlt;
        if (var_diff < 1e-10) {
            break_next = 1;
        }
    }

    return var_A1;
}
