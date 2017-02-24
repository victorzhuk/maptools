#ifndef __LIBMAPTOOLS_H__
#define __LIBMAPTOOLS_H__

/*
 * Объявление типов данных
 */

typedef enum {
    MtGeoidKras,
    MtGeoidWSG80
} mt_geoidtypes_t;

typedef struct
{
    // большая полуось
    double a;
    // малая полуось
    double b;
    // сжатие
    double alpha;
    // квадрат эксцентриситета
    double e_sq;
    // квадрат второго эксцентриситета
    double e2_sq;
} mt_geoid_t;

typedef struct
{
    // широта
    double b;
    // долгота
    double l;
} mt_point_t;


/*
 * Объявление экспортируемых функций
 */

mt_point_t mt_point_move(const mt_point_t *, double /*расстояние*/, double /*азимут*/, mt_geoidtypes_t geoid);
double     mt_distance(const mt_point_t* /*от*/, const mt_point_t* /*до*/, mt_geoidtypes_t geoid);
double     mt_azimuth(const mt_point_t* /*от*/, const mt_point_t* /*до*/, mt_geoidtypes_t geoid);

#endif /* __LIBMAPTOOLS_H__ */
