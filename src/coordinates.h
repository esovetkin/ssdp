#ifndef _COORDINATES_H
#define _COORDINATES_H

#include "epsg.h"


/**
 * point coordinates
 *
 * x is treated as "latitude", y is "longitude"
 */
struct point {
        double x, y;
};


/**
 * 4-nodes polygon
 */
struct poly4 {
        struct point p[4];
};


/**
 * a list of coordinates
 *
 * @p, @nx, @ny, @np: array of points
 *
 * @br: stores bounding box of the coordinates
 */
struct coordinates {
        struct point *p;
        int nx, ny, np;
        struct poly4 br; /* bbox */
        double x1, y1, x2, y2;
};


/**
 * allocate/free memory for coordinates
 */
struct coordinates* coordinates_init(int np);
void coordinates_free(struct coordinates*);


/**
 * convert point between epsg systems
 *
 * @pc: epsg projection context.
 * @p: point to convert.
 * @fwd: positive or negative, defines the direction of conversion.
 *
 * Note, for a given point p, the following is an identify operation.
 *         convert_point(pc, p,  1);
 *         convert_point(pc, p, -1);
 */
void convert_point(struct epsg* pc, struct point* p, int fwd);


/**
 * sample coordinate from a given box
 *
 * @oc: allocate and return point to output coordinates
 *
 * @x1,@y1,@x2,@y2: bbox [lat_min, lon_min, lat_max, lon_max] in
 * epsg:4326.
 *
 * @step: step at which to sample points (in given epsg system)
 *
 * @epsg: epsg system to use for sampling. Use approriate UTM if you
 * want meters.
 */
struct coordinates* box2coordinates(double x1, double y1,
                                    double x2, double y2,
                                    double step, int epsg);


/**
 * write polygon from a bounding box
 *
 * @op: output polygon
 *
 * @a,@b,@c,@d: bbox [lat_min, lon_min, lat_max, lon_max].
 *
 * The order of a,b,c,d is **not** checked and fixed, as the order
 * depends on coordinate system.
 */
void bbox2poly4(struct poly4* op,
                double x1, double y1, double x2, double y2);


/**
 * place template at a location with rotated with azimuth
 *
 * despite a simple formulation, the rotation happens in the given
 * epsg coordinate system. In that coordinate system going up in x
 * coordinate does not necessarily means going north. this function
 * does the check.
 *
 * @pc: epsg projection context.
 * @lat,@lon: latitude and longitude of the origin of the template
 * @azi: azimuth in radians: 0 North, pi/2 East, ...
 * @x,@y,@N: input and output double and its size
 *
 */
void placetemplate(struct epsg* pc, double lat, double lon,
                   double azi, double* x, double* y, int N);

#endif
