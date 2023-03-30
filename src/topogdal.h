#ifndef _TOPOGDAL_H
#define _TOPOGDAL_H

#include <gdal.h>

#include "epsg.h"
#include "coordinates.h"


/**
 * data needed for converting coordinates to pixels and back
 *
 * @d: geotransform (coordinates -> pixel)
 * @i: inverse geotransform (pixel -> coordinates)
 * @nx, @ny: raster shape (x latitude)
 */
struct geotransform {
        double d[6], i[6];
        int nx, ny;
};


/**
 * data for several gdal rasters
 *
 * required for simulations bbox may involve more than one
 * raster. this structure stored only the necessary data allowing to
 * sample those rasters. note, the actual raster data is stored only
 * while sampling and hence is struct `raster` is used for that.
 *
 * @ds: an array of raster bands. We work with elevation rasters (only
 * 1 band), giving rasters with more or less than 1 bands will result
 * in a panic.
 *
 * @pj: coordinate projection for each of the considered rasters. it
 * converts epsg:4326 to the coordinate system of each raster.
 *
 * @br: 4-polygon describing the area covered by raster. in raster's
 * coordinates.
 *
 * @gt: geotransform data
 *
 * @nds: number of raster used
 */
struct gdaldata {
        GDALDatasetH *ds;
        struct epsg **pj;
        struct poly4 *br;
        struct geotransform *gt;
        int nds;
};


/**
 * init/free gdaldata
 *
 * @fns, @nfns: list of paths for filenames
 *
 * returns NULL if:
 *         - any of the raster is not readable
 *         - malloc errors
 *         - rasters have number of bands != 1.
 */
struct gdaldata* gdaldata_init(const char **fns, int nfns);
void gdaldata_free(struct gdaldata*);


/**
 * actual raster data is stored
 *
 * @d, @n: data and its length
 *
 * @xoff, @yoff, @xsize, @ysize: define the subraster that is stored.
 */
struct raster {
        double *d;
        double nodata_value;
        int xoff, yoff, xsize, ysize, n;
};

/**
 * read raster into memory. only the required box part of the raster
 * is read.
 *
 * @gd, @i: gdaldata and the index of the raster.
 *
 * @cb: required polygon in epsg:4326
 *
 * returns NULL when failed to read raster or allocate memory
 */
struct raster* raster_init(struct gdaldata *gd, int i, struct poly4* cb);
void raster_free(struct raster *self);

/**
 * sample georasters at the coordinate locations
 *
 * the function allocated crdnts->np array and returns its
 * pointer. The array contains column-wise points from south-west to
 * north-east.
 */
double* topogrid_from_gdal(struct gdaldata*, struct coordinates*);

#endif
