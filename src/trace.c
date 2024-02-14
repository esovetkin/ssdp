/*
    Simple Sky-Dome Projector Library
    Copyright (C) 2021  B. E. Pieters, 
    IEK-5 Photovoltaik, Forschunszentrum Juelich

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vector.h"
#include "sky_dome.h"
#include "trace.h"
#include "error.h"
#include "print.h"
#include "config.h"
#include "fatan2.h"
#include "rtree.h"


sky_transfer InitSkyTransfer(int N, sky_transfer* src)
{
	sky_transfer T;
	int i;
	T.N = (src) ? src->N : N;
	// ERRORFLAG MALLOCFAILSKYTRANS  "Error memory allocation failed in creating a sky-transfer map"
	if (NULL == (T.t=malloc(T.N*sizeof(*T.t)))) {
		AddErr(MALLOCFAILSKYTRANS);
		T.N=0;
		return T;
	}

	if (src)
			memcpy(T.t, src->t, N*sizeof(*(src->t)));
	else
			for (i=0; i < N; ++i) T.t[i] = (double) 1.0;

	T.g = (src) ? src->g : 0.0;
	return T;
}


void FreeSkyTransfer(sky_transfer *T)
{
	if (T->t)
		free(T->t);
	T->t=NULL;
	T->g=0.0;
	T->N=0;
}


horizon InitHorizon(int Nz)
{
	horizon H;
	horizon H0={
			.zen    = NULL,
			.astep  = 0.0,
			.minz   = 0.0,
			.maxz   = M_PI/2.0,
			.N      = 0
	};
	int i;
	H.N=6*Nz;
	H.astep=2*M_PI/((double)H.N);
	H.minz=0.0;
	H.maxz=M_PI/2.0;
	// ERRORFLAG MALLOCFAILHORIZON  "Error memory allocation failed horizon initialization"
	H.zen=malloc(H.N*sizeof(double));
	if (H.zen==NULL)
	{
		AddErr(MALLOCFAILHORIZON);
		return H0;
	}
	for (i=0;i<H.N;i++)
		H.zen[i]=INFINITY; // initialized at infinity, AtanHorizon must be called once after initializing
	return H;
}


void AtanHorizon(horizon *H)
{
		int i;
		double min = M_PI/2.0, max = 0.0;
		for (i=0; i < H->N; ++i) {
				H->zen[i]=ATAN(H->zen[i]);

				if (min > H->zen[i]) min = H->zen[i];
				if (max < H->zen[i]) max = H->zen[i];
		}

		H->minz = min;
		H->maxz = max;
}


void FreeHorizon(horizon *H)
{
	if (H->zen)
		free(H->zen);
	H->N=0;
	H->zen=NULL;
	H->astep=0;
}

void TransTrans(const sky_transfer *a, const sky_transfer *b, sky_transfer *c)
{
	int i;
	if ((a->N!=b->N)||(a->N!=c->N))
	{
		// ERRORFLAG TRANSFERSIZEMISMATCH  "Error cannot compine different sized transfer functions"
		AddErr(TRANSFERSIZEMISMATCH );
		return;
	}
	for (i=0;i<a->N;i++)
		c->t[i]=a->t[i]*b->t[i];
}

int BelowHorizon(const horizon *H, sky_pos p)
{
	int i,j;
	double a1,a2;
	if (p.a<0)
		p.a+=M_PI*2;
	i=(int)(p.a/H->astep);
	i=i%H->N;
	if (adiff(i*H->astep, p.a)>0)
	{
		i--;
		if (i<0)
			i=H->N-1;
	}
	j=(i+1)%H->N;
	a1=fabs(adiff(p.a, ((double)i)*H->astep));
	a2=H->astep-a1;
	return (p.z>(H->zen[i]*a2+H->zen[j]*a1)/(H->astep));
}


void HorizTrans(const sky_grid *sky, const horizon *a, sky_transfer *c)
{
		int i, imin, imax;

		ZenithRange(sky, a->minz, a->maxz, &imin, &imax);

		// here I assume that sky patches are sorted by:
		//
		//        zenith, azimuth (increasing, first zenith)
		//
		// everything below max is below the horizon for all azi
		//
		// assumption: sizeof(double) zero bytes is 0 in double
		if (imax < c->N)
				memset(c->t+imax, 0, (c->N-imax)*sizeof(*c->t));

		// only call BelowHorizon for patches with zenith angle
		// inside the interval [minz, maxz]
		for (i=imin+1; i < imax; ++i)
				if (BelowHorizon(a, sky->P[i].p))
						c->t[i]=0.0;
}


struct rtreecache* rtreecache_init(double xydelta, double zdelta)
{
		struct rtreecache *self;
		if (NULL==(self=malloc(sizeof(*self)))) goto eself;
		if (NULL==(self->rtree=rtree_new())) goto etree;
		self->xydelta = xydelta/2;
		self->zdelta = zdelta/2;

		return self;
etree:
		free(self);
eself:
		return NULL;
}

#define UNUSED(x) (void)(x)

static bool free_horizon(const double *min, const double *max, const void *item, void *udata)
{
		UNUSED(min); UNUSED(max); UNUSED(udata);
		horizon *h = (horizon *)item;
		FreeHorizon(h);
		free(h);
		return true;
}


void horizoncache_free(struct rtreecache* self)
{
		if (NULL == self)
				return;
		rtree_scan(self->rtree, free_horizon, NULL);
		rtree_free(self->rtree);
		free(self);
}


static bool hz_iter(const double *min, const double *max, const void *item, void *udata)
{
		UNUSED(min); UNUSED(max);
		*((const horizon **) udata) = item;
		return false;
}


horizon* horizoncache_get(struct rtreecache* self, double x, double y, double z)
{
		horizon* res = NULL;

		rtree_search(self->rtree,
					 (double[3]){
							 x-self->xydelta,
							 y-self->xydelta,
							 z-self->zdelta},
					 (double[3]){
							 x+self->xydelta,
							 y+self->xydelta,
							 z+self->zdelta},
					 hz_iter, &res);

		if (NULL!=res)
				return res;

		if (NULL==(res=malloc(sizeof(*res))))
				return NULL;

		res->N = 0;
		res->zen=NULL;
		res->astep=0;

		rtree_insert(self->rtree, (double[3]){x,y,z}, NULL, res);

		return res;
}


static bool free_skytransfer(const double *min, const double *max, const void *item, void *udata)
{
		UNUSED(min); UNUSED(max); UNUSED(udata);
		sky_transfer *st = (sky_transfer*)item;
		FreeSkyTransfer(st);
		free(st);
		return true;
}


void stcache_free(struct rtreecache *self)
{
		if (NULL == self) return;
		rtree_scan(self->rtree, free_skytransfer, NULL);
		rtree_free(self->rtree);
		free(self);
}


static bool st_iter(const double *min, const double *max, const void *item, void *udata)
{
		UNUSED(min); UNUSED(max);
		*((const sky_transfer **) udata) = item;
		return false;
}


sky_transfer* stcache_get(struct rtreecache* self, double x, double y, double z)
{
		sky_transfer* res = NULL;

		rtree_search(self->rtree,
					 (double[3]){
							 x-self->xydelta,
							 y-self->xydelta,
							 z-self->zdelta},
					 (double[3]){
							 x+self->xydelta,
							 y+self->xydelta,
							 z+self->zdelta},
					 st_iter, &res);

		if (NULL!=res)
				return res;

		if (NULL==(res=malloc(sizeof(*res))))
				return NULL;

		res->t=NULL;
		res->N = 0;
		res->g=1.0;

		rtree_insert(self->rtree, (double[3]){x,y,z}, NULL, res);

		return res;
}
