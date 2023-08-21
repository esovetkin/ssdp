#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <proj.h>

#include "epsg.h"

#ifdef RUNMEMTEST
#include "random_fail_malloc.h"
#define malloc(x) random_fail_malloc(x)
#endif


struct epsg* epsg_init_epsg(int src, int dst)
{
        char epsg_src[16], epsg_dst[16];
        snprintf(epsg_src, 16, "EPSG:%d", src);
        snprintf(epsg_dst, 16, "EPSG:%d", dst);

        struct epsg *self = epsg_init(epsg_src, epsg_dst);
        if (NULL == self)
                goto eepsg;

        return self;
eepsg:
        return NULL;
}


struct epsg* epsg_init(const char *epsg_src, const char *epsg_dst)
{
        PJ *norm;
        struct epsg *self;
        self = malloc(sizeof(*self));
        if (NULL == self) goto self_emalloc;

        self->C = proj_context_create();
        if (NULL == self->C) goto ecproj;

        self->P = proj_create_crs_to_crs(self->C, epsg_src, epsg_dst, NULL);
        if (NULL == self->P) goto epproj;

        // This will ensure that the order of coordinates for the
        // input CRS will be longitude, latitude, whereas EPSG:4326
        // mandates latitude, longitude
        norm = proj_normalize_for_visualization(self->C, self->P);
        if (NULL == norm) goto enproj;

        proj_destroy(self->P);
        self->P = norm;

        return self;
enproj:
        proj_destroy(self->P);
epproj:
        proj_context_destroy(self->C);
ecproj:
        free(self);
self_emalloc:
        return NULL;
}


void epsg_free(struct epsg *self)
{
        if (NULL == self)
                return;

        if (self->P)
                proj_destroy(self->P);
        self->P = NULL;

        if (self->C)
                proj_context_destroy(self->C);
        self->C = NULL;
        free(self);
        self = NULL;
}


int determine_utm(double lat, double lon)
{
        int i;

        i = lat >= 0 ? 32600 : 32700;
        i += (int) (fmod(floor((lon + 180) / 6),60) + 1);
        return i;
}


#ifdef RUNTEST

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "coordinates.h"
#include "epsg.h"


double runif(void)
{
        return (double)rand() / (double)RAND_MAX;
}


void test_conversion(struct epsg *pc, struct point p)
{
        struct point p0 = {.x = p.x, .y = p.y};
        convert_point(pc, &p, 1);
        convert_point(pc, &p, -1);
        assert(fabs(p0.x-p.x) < 1e-8 && fabs(p0.y-p.y) < 1e-8);
}


void test_determine_epsg(void)
{
        assert(32632 == determine_utm(50,6.5));
        assert(32732 == determine_utm(-50,6.5));
        assert(32729 == determine_utm(-50,-6.5));
        assert(32631 == determine_utm(0.5,0.5));
        assert(32730 == determine_utm(-0.5,-0.5));
        assert(32660 == determine_utm(31,175));
}


void test_1()
{
        struct epsg *pc = epsg_init_epsg(4326, 3857);
        assert(pc);

        struct point p;
        int i;

        for (i=0; i < 1e6; ++i){
                p.x = runif()*180 - 90;
                p.y = runif()*360 - 180;
                test_conversion(pc, p);
        }

        epsg_free(pc);
}


void test_2()
{
        struct epsg *pc = epsg_init_epsg(4326, 4326);
        assert(pc);

        struct point p, p0;
        int i;

        for (i=0; i<1e6; ++i){
                p.x = runif()*180 - 90;
                p.y = runif()*360 - 180;
                p0.x = p.x;
                p0.y = p.y;
                convert_point(pc, &p, 1);
                assert(fabs(p0.x-p.x) < 1e-8 && fabs(p0.y-p.y) < 1e-8);
        }

        epsg_free(pc);
}


void test_epsg()
{
        struct epsg *pc;
        pc = epsg_init_epsg(-42,4326);
        assert(NULL == pc);
        epsg_free(pc);
        pc = epsg_init("EPSG:4326","");
        assert(NULL == pc);
        epsg_free(pc);
}


int main(void)
{
        printf("testing epsg ...");

        test_determine_epsg();
        test_1();
        test_2();
        test_epsg();

        printf("PASSED\n");
        return 0;
}

#endif
