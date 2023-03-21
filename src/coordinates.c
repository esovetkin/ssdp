#include <stdlib.h>
#include <math.h>

#include "coordinates.h"


static void arange(double *a, int na, double from, double to, double step)
{
        int i;
        step = to > from ? step : -step;

        for (i=0; i < na; ++i)
                a[i] = from + (double) i * step;
}


void bbox2poly4(struct poly4 *br, double x1, double y1, double x2, double y2)
{
        br->p[0].x = x1;
        br->p[0].y = y1;

        br->p[1].x = x2;
        br->p[1].y = y1;

        br->p[2].x = x2;
        br->p[2].y = y2;

        br->p[3].x = x1;
        br->p[3].y = y2;
}


void convert_point(struct epsg *pc, struct point *p, int fwd)
{
        PJ_COORD a, b;
        a = proj_coord(p->y, p->x, 0, 0);
        b = proj_trans(pc->P, (fwd > 0) ? PJ_FWD : PJ_INV, a);
        p->y = b.lp.lam;
        p->x = b.lp.phi;
}


static void convert_array(struct epsg *pc, double *a, int na, double b, int iffirst, int fwd)
{
        int i;
        struct point p;

        for (i=0; i<na; ++i) {
                if (iffirst > 0) {
                        p.x = a[i];
                        p.y = b;
                } else {
                        p.x = b;
                        p.y = a[i];
                }

                convert_point(pc, &p, fwd);
                a[i] = (iffirst > 0) ? p.x : p.y;
        }
}


struct coordinates* coordinates_init(int n)
{
        struct coordinates* self;
        self = malloc(sizeof(*self));
        if (NULL == self)
                goto self_emalloc;

        self->p = malloc(n*sizeof(*self->p));
        if (NULL == self->p)
                goto p_emalloc;
        self->np = n;

        return self;
p_emalloc:
        free(self);
self_emalloc:
        return NULL;
}


void coordinates_free(struct coordinates *self)
{
        if (self->p)
                free(self->p);
        self->p = NULL;
        self->np = 0;
        free(self);
        self = NULL;
}


struct mesh {
        double *lat, *lon;
        int nlat, nlon;
};


static struct mesh* mesh_init(struct epsg *pc,
                              double x1, double y1,
                              double x2, double y2,
                              double step)
{
        struct mesh *self;
        self = malloc(sizeof(*self));
        if (NULL == self)
                goto self_emalloc;


        self->nlat = (int) ceil(fabs(x2-x1)/step);
        self->nlon = (int) ceil(fabs(y2-y1)/step);

        self->lat = malloc(self->nlat*sizeof(*self->lat));
        if (NULL == self->lat)
                goto lat_emalloc;

        self->lon = malloc(self->nlon*sizeof(*self->lon));
        if (NULL == self->lon)
                goto lon_emalloc;

        arange(self->lat, self->nlat, x1, x2, step);
        arange(self->lon, self->nlon, y1, y2, step);

        convert_array(pc, self->lat, self->nlat, y1, 1, -1);
        convert_array(pc, self->lon, self->nlon, x1, 0, -1);

        return self;
lon_emalloc:
        free(self->lat);
lat_emalloc:
        free(self);
self_emalloc:
        return NULL;
}


static void mesh_free(struct mesh *self)
{
        if (self->lat)
                free(self->lat);
        self->lat = NULL;
        self->nlat = 0;

        if (self->lon)
                free(self->lon);
        self->lon = NULL;
        self->nlon = 0;
        free(self);
        self = NULL;
}


static struct coordinates* product_mesh(struct mesh *m)
{
        int i;

        struct coordinates *self = coordinates_init(m->nlat * m->nlon);
        if (NULL == self) goto ecoordinates_init;

        // this defines order of the resulting raster. I want: "column
        // major, from the south-west corner to the north-east corner"
        for (i=0; i < self->np; ++i) {
                self->p[i].x = m->lat[i % m->nlat];
                self->p[i].y = m->lon[i / m->nlat];
        }

        bbox2poly4(&self->br, m->lat[0], m->lon[0],
                   m->lat[m->nlat-1], m->lon[m->nlon-1]);
        self->nx = m->nlat;
        self->ny = m->nlon;

        return self;
ecoordinates_init:
        return NULL;
}




struct coordinates* box2coordinates(double x1, double y1,
                                    double x2, double y2,
                                    double step, int epsg)
{
        struct epsg *pc = epsg_init_epsg(4326, epsg);
        if (NULL == pc) goto eprojc;

        struct point bl, tr;
        bl.x = (x1 < x2) ? x1 : x2;
        bl.y = (y1 < y2) ? y1 : y2;
        tr.x = (x2 > x1) ? x2 : x1;
        tr.y = (y2 > y1) ? y2 : y1;

        convert_point(pc, &bl, 1); convert_point(pc, &tr, 1);

        struct mesh *m = mesh_init(pc, bl.x, bl.y, tr.x, tr.y, step);
        if (NULL == m) goto emesh;

        struct coordinates *self = product_mesh(m);
        if (NULL == self) goto eproduct_mesh;

        self->x1 = bl.x;
        self->y1 = bl.y;
        self->x2 = tr.x;
        self->y2 = tr.y;

        mesh_free(m);
        epsg_free(pc);
        return self;
eproduct_mesh:
        mesh_free(m);
emesh:
        epsg_free(pc);
eprojc:
        return NULL;
}
