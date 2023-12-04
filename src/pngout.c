#include <png.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>
#include "pngout.h"


struct pixel {
        uint8_t red;
        uint8_t green;
        uint8_t blue;
};


struct bitmap {
        struct pixel *pixels;
        size_t width;
        size_t height;
};


static struct pixel* pixel_at(struct bitmap* bm, int x, int y)
{
        return bm->pixels + bm->width * y + x;
}


static uint8_t d2i(double d)
{
        if (d < 0)
                return (uint8_t)0;

        if (d > 255)
                return (uint8_t)255;

        return (uint8_t)d;
}


static void norm(double *src, struct bitmap* dst, double m, double M)
{
        int i,x,y;
        uint8_t v;
        struct pixel *p;

        for (i = 0; i < (int)(dst->height*dst->width); ++i) {
                v = d2i(round(255*(src[i]-m)/(M-m)));
                // transpose to maintain north-west is top-left corner
                // (debugged into existence...)
                x = dst->height - 1 - i % dst->height;
                y = i / dst->height;
                p = pixel_at(dst, y, x);
                p->red = v;
                p->green = v;
                p->blue = v;
        }
}


static void norm_maxmin(double *src, struct bitmap* dst)
{
        int i, N = dst->width * dst->height;
        double M = -DBL_MAX, m = DBL_MAX;
        for (i=0; i < N; ++i) {
                if (src[i] < m)
                        m = src[i];
                if (src[i] > M)
                        M = src[i];
        }

        norm(src, dst, m, M);
}


static struct bitmap* bitmap_init(double *z, int nx, int ny, enum normalisation normalisation)
{
        struct bitmap* self;
        self = malloc(sizeof(*self));
        if (NULL == self) goto eself;

        self->width = nx;
        self->height = ny;
        self->pixels = malloc(nx*ny*sizeof(*self->pixels));
        if (NULL == self->pixels) goto epixels;

        switch (normalisation) {
        case NORM_NONE:
                norm(z, self, 0.0, 255.0);
                break;
        case NORM_MAXMIN:
                norm_maxmin(z, self);
                break;
        default:
                norm_maxmin(z, self);
                break;
        }

        return self;
        free(self->pixels);
epixels:
        free(self);
eself:
        return NULL;
}


static void bitmap_free(struct bitmap* self)
{
        if (NULL == self)
                return;
        if (self->pixels)
                free(self->pixels);
        self->pixels = NULL;
        self->width = 0;
        self->height = 0;
        free(self);
        self = NULL;
}


static int save_png(struct bitmap *bm, const char *ofn)
{
        FILE * fp;
        png_structp png_ptr = NULL;
        png_infop info_ptr = NULL;
        size_t x, y;
        png_byte ** row_pointers = NULL;

        /* The following number is set by trial and error only. I cannot
           see where it it is documented in the libpng manual.
        */
        int pixel_size = 3;
        int depth = 8;

        fp = fopen(ofn, "wb");
        if (!fp) goto efopen;

        png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (NULL == png_ptr) goto epngptr;

        info_ptr = png_create_info_struct(png_ptr);
        if (NULL == info_ptr) goto einfoptr;

        /* Set up error handling. */
        if (setjmp (png_jmpbuf (png_ptr))) goto epng_failure;

        /* Set image attributes. */
        png_set_IHDR(png_ptr, info_ptr, bm->width, bm->height,
                     depth, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

        /* Initialize rows of PNG. */
        row_pointers = png_malloc(png_ptr, bm->height*sizeof(png_byte*));
        for (y = 0; y < bm->height; y++) {
                png_byte *row = png_malloc(png_ptr, sizeof(uint8_t)*bm->width*pixel_size);
                row_pointers[y] = row;
                for (x = 0; x < bm->width; x++) {
                        struct pixel * pixel = pixel_at(bm, x, y);
                        *row++ = pixel->red;
                        *row++ = pixel->green;
                        *row++ = pixel->blue;
                }
        }

        /* Write the image data to "fp". */
        png_init_io(png_ptr, fp);
        png_set_rows(png_ptr, info_ptr, row_pointers);
        png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

        // TODO shouldn't this be caught after epng_failure
        for (y = 0; y < bm->height; y++)
                png_free(png_ptr, row_pointers[y]);
        png_free(png_ptr, row_pointers);

        png_destroy_write_struct (&png_ptr, &info_ptr);
        fclose(fp);
        return 0;
epng_failure:
einfoptr:
        png_destroy_write_struct(&png_ptr, &info_ptr);
epngptr:
        fclose(fp);
efopen:
        return -1;
}


int write_png(const char* ofn, double *z, int nx, int ny,
              enum normalisation norm)
{
        struct bitmap* bm = bitmap_init(z,nx,ny,norm);
        if (NULL == bm) goto ebm;

        if (save_png(bm, ofn)) goto esave;

        bitmap_free(bm);
        return 0;
esave:
        bitmap_free(bm);
ebm:
        return -1;
}


#ifdef RUNTEST

#include <stdio.h>
#include <assert.h>


void test_0(int nx, int ny)
{
        int i;
        double *z;
        z = malloc(nx*ny*sizeof(*z));
        if (NULL == z) goto ez;

        for (i=0; i < nx*ny; ++i)
                z[i] = (double)i;

        write_png("test_0.png",z,nx,ny,NORM_NONE);
        write_png("test_1.png",z,nx,ny,NORM_MAXMIN);

        free(z);
        return;
        free(z);
ez:
        assert(0);
}


int main(void)
{
        printf("testing pngout ...");
        test_0(100,100);
        printf("PASSED\n");
}

#endif
