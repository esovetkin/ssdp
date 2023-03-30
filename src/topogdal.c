#include <cpl_conv.h> /* for CPLMalloc() */
#include <math.h>

#include "error.h"
#include "epsg.h"
#include "topogdal.h"
#include "edt.h"


void save_raster(double *d, int nd)
{
        int i;
        FILE *fd;
        fd = fopen("raster.txt", "w");
        if (NULL == fd) goto efopen;

        for (i=0; i < nd; ++i)
                fprintf(fd, "%.5f\n", d[i]);
        fclose(fd);
efopen:
        return;
}


static struct poly4 raster_corners(double gt[6], int xs, int ys)
{
        struct poly4 br;
        bbox2poly4(&br, gt[3], gt[0], gt[3]+ys*gt[5], gt[0]+xs*gt[1]);
        return br;
}


static int set_geotransform(GDALDatasetH *ds, struct geotransform *gt)
{
        GDALGetGeoTransform(*ds, gt->d);

        if (0 == GDALInvGeoTransform(gt->d, gt->i))
                goto einvtransform;

        gt->nx = GDALGetRasterYSize(*ds);
        gt->ny = GDALGetRasterXSize(*ds);
        return 0;
einvtransform:
        return 1;
}


// ERRORFLAG GDALNCHANNELNOTONE "Provided raster has channels != 1"
// ERRORFLAG GDALREADFAILED "Cannot read provided raster"
struct gdaldata* gdaldata_init(const char *fns[], int nfs)
{
        struct gdaldata* self;
        self = malloc(sizeof(*self));
        if (NULL == self)
                goto self_emalloc;

        int i;
        GDALAllRegister();

        if (NULL == (self->ds = malloc(nfs*sizeof(*self->ds))))
                goto ds_emalloc;
        if (NULL == (self->br = malloc(nfs*sizeof(*self->br))))
                goto br_emalloc;
        if (NULL == (self->gt = malloc(nfs*sizeof(*self->gt))))
                goto gt_emalloc;

        self->nds = 0;
        for (i=0; i < nfs; ++i) {
                if (NULL == (self->ds[i] = GDALOpen(fns[i], GA_ReadOnly))) {
                        AddErr(GDALREADFAILED);
                        goto loop_egdal;
                }

                self->nds = i + 1;

                // raise problem, when #channels != 1. it supposed to
                // be elevation raster after all...
                if (1 != GDALGetRasterCount(self->ds[i])) {
                        AddErr(GDALNCHANNELNOTONE);
                        goto loop_egdal;
                }

                if (set_geotransform(&self->ds[i], &self->gt[i]))
                        goto loop_egdal;

                self->br[i] = raster_corners(self->gt[i].d, self->gt[i].nx, self->gt[i].ny);
        }

        int nepj = 0;
        if (NULL == (self->pj = malloc(nfs*sizeof(*self->pj))))
                goto loop_epj;
        for (i=0; i < nfs; ++i) {
                self->pj[i] = epsg_init("EPSG:4326",GDALGetProjectionRef(self->ds[i]));
                if (NULL == self->pj[i])
                        goto loop_epj;
                ++nepj;
        }

        return self;
loop_epj:
        for (i=0; i < nepj; ++i)
                epsg_free(self->pj[i]);
        free(self->pj);
loop_egdal:
        for (i=0; i < self->nds; ++i)
                GDALClose(self->ds[i]);
        free(self->gt);
gt_emalloc:
        free(self->br);
br_emalloc:
        free(self->ds);
ds_emalloc:
        free(self);
self_emalloc:
        return NULL;
}


void gdaldata_free(struct gdaldata* self) {
        int i;

        if (NULL == self)
                return;

        if (self->ds) {
                for (i=0; i < self->nds; ++i)
                        GDALClose(self->ds[i]);
                free(self->ds);
                self->ds = NULL;
        }

        if (self->pj) {
                for (i=0; i < self->nds; ++i)
                        epsg_free(self->pj[i]);
                free(self->pj);
                self->pj = NULL;
        }

        if (self->br)
                free(self->br);
        self->br = NULL;

        if (self->gt)
                free(self->gt);
        self->gt = NULL;

        self->nds = 0;
        free(self);
        self = NULL;
}


static inline double minp(struct poly4 *p, char what)
{
        double a,b,c,d;
        if ('x' == what) {
                a = p->p[0].x;
                b = p->p[1].x;
                c = p->p[2].x;
                d = p->p[3].x;
        } else {
                a = p->p[0].y;
                b = p->p[1].y;
                c = p->p[2].y;
                d = p->p[3].y;
        }

        a = (a < b) ? a : b;
        c = (c < d) ? c : d;
        return (a < c) ? a : c;
}


static inline double maxp(struct poly4 *p, char what)
{
        double a,b,c,d;
        if ('x' == what) {
                a = p->p[0].x;
                b = p->p[1].x;
                c = p->p[2].x;
                d = p->p[3].x;
        } else {
                a = p->p[0].y;
                b = p->p[1].y;
                c = p->p[2].y;
                d = p->p[3].y;
        }

        a = (a > b) ? a : b;
        c = (c > d) ? c : d;
        return (a > c) ? a : c;
}


static void conv_box(int out[4], struct gdaldata *gd, int i,
                     struct poly4 cb, int xsize, int ysize)
{
        convert_point(gd->pj[i], &cb.p[0], 1);
        convert_point(gd->pj[i], &cb.p[1], 1);
        convert_point(gd->pj[i], &cb.p[2], 1);
        convert_point(gd->pj[i], &cb.p[3], 1);

        struct poly4 pb;
        GDALApplyGeoTransform(gd->gt[i].i, cb.p[0].y, cb.p[0].x, &pb.p[0].y, &pb.p[0].x);
        GDALApplyGeoTransform(gd->gt[i].i, cb.p[1].y, cb.p[1].x, &pb.p[1].y, &pb.p[1].x);
        GDALApplyGeoTransform(gd->gt[i].i, cb.p[2].y, cb.p[2].x, &pb.p[2].y, &pb.p[2].x);
        GDALApplyGeoTransform(gd->gt[i].i, cb.p[3].y, cb.p[3].x, &pb.p[3].y, &pb.p[3].x);

        out[0] = minp(&pb, 'x');
        out[1] = minp(&pb, 'y');
        out[2] = maxp(&pb, 'x') - out[0] + 1;
        out[3] = maxp(&pb, 'y') - out[1] + 1;

        if (out[0] < 0)
                out[0] = 0;
        if (out[0] > ysize)
                out[0] = ysize;
        if (out[1] < 0)
                out[1] = 0;
        if (out[1] > xsize)
                out[1] = xsize;
        if (out[2] + out[0] > ysize)
                out[2] = ysize - out[0];
        if (out[3] + out[1] > xsize)
                out[3] = xsize - out[1];
}


// ERRORFLAG GDALRASTEREMALLOC "Cannot allocate memory for the raster"
// ERRORFLAG GDALRASTEREREAD "Error on GDALRasterIO!"
struct raster* raster_init(struct gdaldata *gd, int i, struct poly4 *cb)
{
        struct raster *self;
        self = malloc(sizeof(*self));
        if (NULL == self)
                goto self_emalloc;

        GDALRasterBandH hband;
        hband = GDALGetRasterBand(gd->ds[i], 1);

        int rb[4];
        conv_box(rb, gd, i, *cb,
                 GDALGetRasterBandXSize(hband),
                 GDALGetRasterBandYSize(hband));

        self->xoff = (int) rb[0];
        self->yoff = (int) rb[1];
        self->xsize = (int) rb[2];
        self->ysize = (int) rb[3];
        self->n = self->xsize * self->ysize;

        // special case when raster has no data
        if (0 == self->n) {
                self->d = NULL;
                return self;
        }

        self->d = (double*)CPLMalloc(self->n*sizeof(*self->d));
        if (NULL == self->d) {
                AddErr(GDALRASTEREMALLOC);
                goto ecplmalloc;
        }

        if (GDALRasterIO(hband, GF_Read,
                         self->yoff, self->xoff,
                         self->ysize, self->xsize,
                         self->d,
                         self->ysize, self->xsize,
                         GDT_Float64, 0, 0 )) {
                AddErr(GDALRASTEREREAD);
                goto erasterio;
        }

        return self;
erasterio:
        CPLFree(self->d);
ecplmalloc:
        free(self);
self_emalloc:
        return NULL;
}


void raster_free(struct raster *self)
{
        if (self->d)
                CPLFree(self->d);
        self->d = NULL;
        self->xsize = 0;
        self->ysize = 0;
        free(self);
        self = NULL;
}


static double get_pixel(struct geotransform *gt, struct raster *r, struct point *p)
{
        double u, v;

        // case when raster has no data
        if (0 == r->n)
                return -9999.0;

        GDALApplyGeoTransform(gt->i, p->y, p->x, &v, &u);
        u -= r->xoff;
        v -= r->yoff;

        if (u < 0 || u >= r->xsize || v < 0 || v >= r->ysize)
                return -9999.0;

        return r->d[(int) u * r->ysize + (int) v];
}


static double lookup_raster(struct gdaldata *gd, struct raster **r, struct point p)
{
        int i;
        double x, val = -9999.0;

        for (i=0; i < gd->nds; ++i) {
                convert_point(gd->pj[i], &p, 1);

                x = get_pixel(&gd->gt[i], r[i], &p);
                val = x > val ? x : val;

                convert_point(gd->pj[i], &p, -1);
        }

        return val;
}


void free_rasters(struct raster **rd, int nrd)
{
        int i;
        for (i=0; i<nrd; ++i)
                raster_free(rd[i]);
        free(rd);
}


double* topogrid_from_gdal(struct gdaldata *gd, struct coordinates *locs)
{
        double *z;
        if (NULL == (z = malloc(locs->np*sizeof(*z))))
                goto z_emalloc;

        struct raster **rd;
        if (NULL == (rd = malloc(gd->nds*sizeof(*rd))))
                goto rd_emalloc;

        // read patches of rasters in memory
        int i, nrd = 0;
        for (i=0; i<gd->nds; ++i) {
                rd[i] = raster_init(gd, i, &locs->br);
                if (NULL == rd[i]) goto eraster_init;
                nrd = i + 1;
        }

        for (i=0; i<locs->np; ++i)
                z[i] = lookup_raster(gd, rd, locs->p[i]);

        free_rasters(rd, nrd);

        // fill the missing values
        struct edt *dt = edt_init(z, locs->np, locs->nx, locs->ny);
        if (NULL == dt) goto eedt_init;
        edt_compute(dt);
        edt_fill(dt);
        edt_free(dt);

        return z;
eedt_init:
eraster_init:
        free_rasters(rd, nrd);
rd_emalloc:
        free(z);
z_emalloc:
        return NULL;
}
