#ifndef _EPSG_H
#define _EPSG_H

#include <proj.h>

struct epsg {
        PJ_CONTEXT *C;
        PJ *P;
};


/**
 * create epsg context from epsg's integer codes
 *
 * @src, @dst: integer epsg codes, e.g. 4326 is the GPS longitude and
 * latitude. src and dst affect the direction of transformation, see
 * also coordinates.h:convert_point.
 */
struct epsg* epsg_init_epsg(int, int);


/**
 * create epsg context from projection string
 *
 * @epsg_src, @epsg_dst: any strings that proj::proj_create_crs_to_crs
 * can understand
 */
struct epsg* epsg_init(const char*, const char*);


void epsg_free(struct epsg*);


/**
 * get utm epsg code from a point
 *
 * @lat, @lon: coordinate in epsg:4326
 */
int determine_utm(double lat, double lon);

#endif
