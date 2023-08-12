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


struct coordinates* coordinates_init(int n)
{
        if (n <= 0)
                goto self_emalloc;

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
        if (NULL == self)
                return;

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


static struct mesh* mesh_init(double x1, double y1,
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


static struct coordinates* product_mesh(struct epsg *pc, struct mesh *m)
{
        int i;

        struct coordinates *self = coordinates_init(m->nlat * m->nlon);
        if (NULL == self) goto ecoordinates_init;

        // this defines order of the resulting raster. I want: "column
        // major, from the south-west corner to the north-east corner"
        for (i=0; i < self->np; ++i) {
                self->p[i].x = m->lat[i % m->nlat];
                self->p[i].y = m->lon[i / m->nlat];
				convert_point(pc, self->p + i, -1);
        }

		struct point p = { .x = m->lat[0], .y = m->lon[0]},
                q = { .x = m->lat[m->nlat-1], .y = m->lon[m->nlon-1]};

		convert_point(pc, &p, -1);
		convert_point(pc, &q, -1);

        bbox2poly4(&self->br, p.x, p.y, q.x, q.y);
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

        struct mesh *m = mesh_init(bl.x, bl.y, tr.x, tr.y, step);
        if (NULL == m) goto emesh;

        struct coordinates *self = product_mesh(pc, m);
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


void placetemplate(struct epsg* pc, double lat, double lon,
                   double azi, double* x, double* y, int N)
{
        // q points north from origin
        struct point p = { .x = lat, .y = lon},
                q = { .x = lat + 0.001, .y = lon};
        convert_point(pc, &p, -1); convert_point(pc, &q, -1);

        // adjust bearing to the coordinate system
        azi += atan2(q.y-p.y, q.x-p.x);

        double c = cos(azi), s = sin(azi);
        for (int i=0; i < N; ++i) {
                q.x = p.x + x[i]*c - y[i]*s;
                q.y = p.y + x[i]*s + y[i]*c;
                convert_point(pc, &q, 1);
                x[i] = q.x; y[i] = q.y;
        }
}


#ifdef RUNTEST

#include <stdio.h>
#include <assert.h>


void test_coords_init(void)
{
        struct coordinates* p;
        p = coordinates_init(100);
        assert(p);
        coordinates_free(p);

        p = coordinates_init(0);
        assert(NULL == p);
        coordinates_free(p);
}


void test_box2coordinates(void)
{
        int epsg;
        epsg = determine_utm(-34.498, 149.491);

        struct coordinates* p;
        p = box2coordinates(-34.498, 149.491, -32.296, 151.754, 100, epsg);
        assert(p);
        coordinates_free(p);

        p = box2coordinates(-34.498, 149.491, -32.296, 151.754, 500, epsg);
        assert(p);
        coordinates_free(p);

        p = box2coordinates(-34.498, 149.491, -32.296, 151.754, 1500, epsg);
        assert(p);
        coordinates_free(p);
}


double _bearing(double lat0, double lon0, double lat1, double lon1)
{
        lon0 *= M_PI/180;
        lon1 *= M_PI/180;
        lat0 *= M_PI/180;
        lat1 *= M_PI/180;

        double dLon = lon1 - lon0;

        return atan2(sin(dLon)*cos(lat1),
                     cos(lat0)*sin(lat1) - sin(lat0)*cos(lat1)*cos(dLon));
}


void _test_placetemplate(double lat, double lon, double azi)
{
        double x[2] = {0.,1.};
        double y[2] = {0.,1.};
        int N = 2, epsg = determine_utm(lat, lon);
        struct epsg *pc = epsg_init_epsg(epsg, 4326);
        if (NULL == pc) goto eepsg;

        double v = azi + atan2(x[1]-x[0], y[1]-y[0]);

        placetemplate(pc,lat,lon,azi,x,y,N);

        v -= _bearing(x[0],y[0],x[1],y[1]);
        assert(fmod(v,2*M_PI) > 2*M_PI - 0.01 || fmod(v,2*M_PI) < 0.01);
eepsg:
        epsg_free(pc);
}


void test_placetemplate(void)
{
        int i,j;
        double b = 0;
        double lat[4] = {5, -5, 50, -50};
        double lon[4] = {5, -5, 50, -50};

        for (i=0; i < 4; ++i)
                for (j=0; j<4; ++j)
                        for (b = 0; b < 2*M_PI; b+=2*M_PI/60)
                                _test_placetemplate(lat[i], lon[j], b);
}


int main(void)
{
        printf("testing coordinates ...");

        test_coords_init();
        test_box2coordinates();
        test_placetemplate();

        printf("PASSED\n");
        return 0;
}


#endif
