#include <stdlib.h>

#include "edt.h"


struct edt* edt_init(double *z, int n, int nx, int ny)
{
        struct edt *self;
        self = malloc(sizeof(*self));
        if (NULL == self) goto edt_eself;

        self->z = z;
        self->n = n;
        self->nx = nx;
        self->ny = ny;
        self->missing_value = -9000;

        self->x = malloc(n*sizeof(*self->x));
        if (NULL == self->x) goto edt_ex;

        self->xg = malloc(n*sizeof(*self->xg));
        if (NULL == self->xg) goto edt_exg;

        self->G = malloc(n*sizeof(*self->G));
        if (NULL == self->G) goto edt_eG;

        self->s = malloc(nx*sizeof(*self->s));
        if (NULL == self->s) goto edt_es;

        self->t = malloc(nx*sizeof(*self->t));
        if (NULL == self->t) goto edt_et;

        return self;
edt_et:
        free(self->s);
edt_es:
        free(self->G);
edt_eG:
        free(self->xg);
edt_exg:
        free(self->x);
edt_ex:
        free(self);
edt_eself:
        return NULL;
}


void edt_free(struct edt *self)
{
        free(self->t);
        free(self->s);
        free(self->G);
        free(self->xg);
        free(self->x);
        free(self);
        self = NULL;
}


void edt_fill(struct edt *self)
{
        int i;
        for (i=0; i<self->n; ++i) {
                if (self->z[i] >= self->missing_value)
                        continue;

                if (self->x[i] < 0 || self->x[i] >= self->n)
                        // TODO print bug warning
                        continue;

                self->z[i] = self->z[self->x[i]];
        }
}


static void scan1(struct edt *self, int x)
{
        if (self->z[x] >= self->missing_value) {
                self->G[x] = 0;
                self->xg[x] = x;
        } else {
                self->G[x] = -1;
                self->xg[x] = -1;
        }

        int y, i, j;
        for (y=1; y < self->ny; ++y) {
                i = x + y*self->nx;

                if (self->z[i] >= self->missing_value) {
                        self->G[i] = 0;
                        self->xg[i] = i;
                        continue;
                }

                j = x+(y-1)*self->nx;
                if (-1 == self->G[j]) {
                        self->G[i] = -1;
                        self->xg[i] = -1;
                        continue;
                }

                self->G[i] = 1 + self->G[j];
                self->xg[i] = self->xg[j];
        }
}


static void scan2(struct edt *self, int x)
{
        int y, i, j;
        for (y=self->ny-2; y>=0; --y) {
                i = x+y*self->nx;
                j = x+(y+1)*self->nx;

                if (-1 != self->G[j] && -1 == self->G[i]) {
                        self->G[i] = 1 + self->G[j];
                        self->xg[i] = self->xg[j];
                        continue;
                }

                if (-1 != self->G[j] && self->G[j] < self->G[i]) {
                        self->G[i] = 1 + self->G[j];
                        self->xg[i] = self->xg[j];
                        continue;
                }
        }
}


static void computeG(struct edt *self)
{
        // compute xg, G
        int x;
        for (x=0; x < self->nx; ++x) {
                scan1(self, x);
                scan2(self, x);
        }
}


static int f(struct edt *self, int y, int x, int i)
{
        int gi = self->G[i + y*self->nx];

        if (-1 == gi)
                return -1;

        return (x-i)*(x-i) + gi*gi;
}


static int sep(struct edt *self, int y, int i, int u)
{
        int gu = self->G[u + y*self->nx];
        int gi = self->G[i + y*self->nx];

        if (-1 == gu || -1 == gi)
                return -1;

        return (u*u - i*i + gu*gu - gi*gi) / (2*(u-i));
}


static int scan3(struct edt *self, int y)
{
        int w, u, q = 0, *s=self->s, *t=self->t;
        s[0] = 0; t[0] = 0;

        for (u=1; u < self->nx; ++u) {
                while (q >= 0) {
                        if (-1 == f(self,y,t[q],u))
                                break;

                        if (f(self,y,t[q],s[q]) <= f(self,y,t[q],u)
                            && -1 != f(self,y,t[q],u)
                            && -1 != f(self,y,t[q],s[q]))
                                break;

                        --q;
                }

                if (q<0) {
                        q = 0;
                        s[0] = u;
                } else {
                        w = 1 + sep(self,y,s[q],u);
                        if (0 != w && w < self->nx) {
                                ++q;
                                s[q] = u;
                                t[q] = w;
                        }
                }
        }

        return q;
}


static void scan4(struct edt *self, int y, int q)
{
        int u, i;
        int *s = self->s, *t = self->t;

        for (u = self->nx-1; u >= 0; --u) {
                i = u + y*self->nx;
                // self->D[i] = f(self,y,u,s[q]);
                self->x[i] = self->xg[s[q] + y*self->nx];
                if (u == t[q])
                        --q;
        }
}


static void computeX(struct edt *self)
{
        int y, q;

        for (y=0; y < self->ny; ++y) {
                q = scan3(self, y);
                scan4(self, y, q);
        }
}


void edt_compute(struct edt* self)
{
        computeG(self);
        computeX(self);
}


#ifdef RUNTEST

#include <stdio.h>
#include <math.h>
#include <assert.h>


void print_a(void *a, int sa, int nx, int ny, const char *name)
{
        int i, row, n = nx*ny;
        printf("%s:\n", name);
        for (row=nx-1; row >= 0; --row) {
                for (i=0; i<n; ++i)
                        if (i % nx == row) {
                                if (sizeof(double) == sa)
                                        printf(" %8.0f", *(double*)(a + i*sa));
                                else
                                        printf(" %8d", *(int*)(a + i*sa));
                        }
                printf("\n");
        }
        printf("\n");
}


void assert_zero(struct edt *dt)
{
        int i;
        for (i=0; i<dt->n; ++i)
                assert(fabs(dt->z[i]) < 1e-8);
}


void assert_nonmissing(struct edt *dt)
{
        int i;
        for (i=0; i<dt->n; ++i) {
                assert(fabs(dt->z[i]+1) > 1e-8);
                assert(-1 != dt->x[i]);
        }
}


void test_1(int nx, int ny, int every)
{
        int i, n = nx*ny;
        double *z = malloc(nx*ny*sizeof(*z));
        assert(z);

        struct edt *dt = edt_init(z, n, nx, ny);
        assert(dt);

        for (i=0; i<n; ++i)
                z[i] = (i % every ? -9999.0 : 0.0);

        edt_compute(dt);
        edt_fill(dt);
        assert_zero(dt);

        for (i=0; i<n; ++i)
                z[i] = (i % every ? 0.0 : -9999.0);

        edt_compute(dt);
        edt_fill(dt);
        assert_zero(dt);

        edt_free(dt);
        free(z);
}


void test_2()
{
        int i, nx = 4, ny = 4, every = 2, n = nx*ny;
        double *z = malloc(nx*ny*sizeof(*z));
        assert(z);

        for (i=0; i<n; ++i)
                z[i] = (i % every ? (double)i : -9999.0);

        struct edt *dt = edt_init(z, n, nx, ny);
        assert(dt);

        edt_compute(dt);
        edt_fill(dt);

        assert_nonmissing(dt);

        double sum = 0;
        for (i=0; i<n; ++i)
                sum += z[i];

        assert((sum - 120) < 1e-8);
        edt_free(dt);
        free(z);
}


void test_corners(int nx, int ny)
{
        nx *= 2; ny *= 2;
        int i, n = nx*ny;
        double *z, *ndx;
        z = malloc(nx*ny*sizeof(*z));
        assert(z);
        ndx = malloc(nx*ny*sizeof(*ndx));
        assert(ndx);

        struct edt *dt = edt_init(z, n, nx, ny);
        assert(dt);

        for (i=0; i<n; ++i) {
                ndx[i] = (double) i;

                if (0 == i || i == nx - 1 || i == n - nx || i == n - 1) {
                        z[i] = (double) i;
                        continue;
                }

                z[i] = -9999.0;
        }

        edt_compute(dt);
        edt_fill(dt);

        assert_nonmissing(dt);

        double sum = 0;
        for (i=0; i<n; ++i)
                sum += z[i];
        assert(fabs(sum-(n-1 + n-nx + nx-1)*((double)(n/4))) < 1e-8);

        edt_free(dt);
        free(ndx);
        free(z);
}

int main(void)
{
        printf("testing fillmissing ...\n");

        test_1(6,6,2);
        test_1(6,6,7);
        test_1(5,5,3);
        test_1(100,5,2);
        test_1(2,100,2);
        test_1(100,100,3);

        test_corners(5,5);
        test_corners(100,2);
        test_corners(2,100);
        test_corners(100,100);

        test_2();

        printf("PASSED\n");
}

#endif
