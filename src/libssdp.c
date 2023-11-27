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
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "vector.h"
#include "sky_dome.h"
#include "trace.h"
#include "sky_model.h"
#include "project.h"
#include "delaunay.h"
#include "ground.h"
#include "print.h"
#include "error.h"
#include "sunpos.h"
#include <config.h>

#ifdef RUNMEMTEST
#include "random_fail_malloc.h"
#define malloc(x) random_fail_malloc(x)
#endif

/* libssdp entry points */
void ssdp_print_version(void)
{
	printf("%s\n", PACKAGE_STRING);
}
void ssdp_print_error_messages()
{
	E_Messages();
}
void ssdp_reset_errors()
{
	ResetErrors();
}
	
/* skydome routines */
int ssdp_find_skypatch(sky_grid *sky, sky_pos p)
{
	return FindPatch(sky, p);
}
sky_grid* ssdp_init_sky(int Nz, int ncopies)
{
		int i;
		sky_grid* self;
		if (NULL==(self=malloc(ncopies * sizeof(*self))))
				goto eself;

		// TODO not optimal, deep copy instead
		for (i=0; i < ncopies; ++i)
				self[i] = InitSky(Nz);

		return self;
eself:
		return NULL;
}
void ssdp_free_sky(sky_grid *sky, int nskies)
{
		if (NULL == sky)
				return;

		int i;
		for (i=0; i < nskies; ++i)
				free_sky_grid(sky+i);

		free(sky);
}
/* skymodel routines */
void ssdp_make_uniform_sky(sky_grid *sky, sky_pos sun, double GHI, double DHI)
{
	UniformSky(sky, sun, GHI, DHI);
}
void ssdp_make_sky_sunonly(sky_grid *sky, sky_pos sun, double GHI, double DHI)
{
	SkySunOnly(sky, sun, GHI, DHI);
}
void ssdp_make_perez_all_weather_sky(sky_grid * sky, sky_pos sun, double GHI, double DHI, double dayofyear)
{
	PerezSky(sky, sun, GHI, DHI, dayofyear);
}
void ssdp_make_uniform_sky_coordinate(sky_grid *sky, time_t t, double lon, double lat, double E, double p, double T, double GHI, double DHI)
{
	sky_pos sun;	
	sun=sunpos(t, lat, lon, E, p, T);
	UniformSky(sky, sun, GHI, DHI);
}
void ssdp_make_skysunonly_coordinate(sky_grid *sky, time_t t, double lon, double lat, double E, double p, double T, double GHI, double DHI)
{
	sky_pos sun;	
	sun=sunpos(t, lat, lon, E, p, T);
	SkySunOnly(sky, sun, GHI, DHI);
}
void ssdp_make_perez_all_weather_sky_coordinate(sky_grid * sky, time_t t, double lon, double lat, double E, double p, double T, double GHI, double DHI)
{
	sky_pos sun;
	struct tm * ut;
	ut=gmtime(&t);
	if (!ut)
	{
		AddErr(GMTIMENULL); //GMTIMENULL initialized in sunpos.c
		return;
	}
	sun=sunpos(t, lat, lon, E, p, T);
	PerezSky(sky, sun, GHI, DHI, (double)ut->tm_yday);
}
void ssdp_make_perez_cumulative_sky_coordinate(sky_grid * sky, double *t, double lon, double lat, double E, double *p, int np, double *T, int nT, double *GHI, double *DHI, int N)
{
	int i;
	sky_pos *sun;
	double *dayofyear;
	time_t tt;
	struct tm * ut;
	sun=malloc(N*sizeof(sky_pos));
	dayofyear=malloc(N*sizeof(double));
	for (i=0;i<N;i++)
	{
		tt=(time_t)t[i];
		ut=gmtime(&tt);
		if (!ut)
		{
			AddErr(GMTIMENULL); //GMTIMENULL initialized in sunpos.c
			return;
		}
		sun[i]=sunpos(tt, lat, lon, E, p[i%np], T[i%nT]);
		dayofyear[i]=(double)ut->tm_yday;
	}
	CumulativePerezSky(sky, sun, t, GHI, DHI, dayofyear, N);
	free(sun);
	free(dayofyear);
}
// todo add cumulative sky routine
/* project routines */
// orient module w.r.t. ground orientation
void ssdp_poa_to_surface_normal(sky_pos pn0, sky_pos sn, sky_pos *pn)
{
	(*pn)=pn0;
	POA_to_SurfaceNormal(pn, sn);
}
AOI_Model_Data ssdp_init_aoi_model(AOI_Model model,double nf, double nar,double *theta, double *effT, int N)
{
	return InitAOIModel(model,nf,nar,theta,effT,N);
}
//BEGIN_SSDP_EXPORT
typedef struct location {
	sky_transfer T;
	horizon* H; // computed horizons are stored in horizoncache
	double difftrans; // diffuse light transmission efficiency for a uniform sky
} location;
//END_SSDP_EXPORT

void ssdp_free_location(location *l)
{
	if (l)
	{
		FreeSkyTransfer(&(l->T));
		// horizoncache takes care of that
		// FreeHorizon(l->H);
	}
}

#define ALBEPS 1e-10


int ssdp_init_transfer(
		sky_transfer* st, sky_grid *sky, double albedo,
		sky_pos pn, AOI_Model_Data *M)
{
		// corresponds to empty locations in the uST array
		if (NULL == st) return 0;

		// sky_transfer has been initialised
		if (NULL != st->t) return 0;

		// setup transfer, sky direct
		*st = InitSkyTransfer(sky->N, NULL);

		if (ssdp_error_state) goto einit;
		POA_Sky_Transfer(sky, st, pn, M);
		if (ssdp_error_state) goto error;

		int i;
		if (albedo>ALBEPS) {
				st->g=albedo*POA_Albedo_Transfer(sky, pn, M);
				if (ssdp_error_state) goto error;

				for (i=0; i<st->N; ++i)
						st->t[i]*=(1.0+st->g);
		}

		return 0;
error:
		FreeSkyTransfer(st);
einit:
		return -1;
}


static int transfer_sky(location *l, sky_grid *sky, sky_transfer *initsky)
{
		int i;

		l->T=InitSkyTransfer(sky->N, initsky);

		HorizTrans(sky, l->H, &(l->T), &(l->T));
		if (ssdp_error_state) {
				Print(WARNING, "ERROR: Horizon does not match the sky\n");
				goto error;
		}

		l->difftrans=0;
		for (i=0;i<l->T.N;i++)
				l->difftrans+=sky->sa[i]*l->T.t[i];
		l->difftrans/=sky->icosz;

		return 0;
error:
		FreeSkyTransfer(&(l->T));
		return -1;
}

int ssdp_setup_horizon(
		horizon *h, sky_grid *sky, topology *T,
		double xoff, double yoff, double zoff)
{
		// corresponds to empty locations in the uH array
		if (NULL == h)
				return 0;

		// horizon has been initialised
		if (NULL != h->zen)
				return 0;

		*h = InitHorizon(sky->Nz);
		if (ssdp_error_state) goto error;

		if (NULL == T) return 0;

		// TODO: should make minzen a setting
		ComputeHorizon(h, T, 0, xoff, yoff, zoff);
		AtanHorizon(h);
		if (ssdp_error_state) goto error;

		return 0;
error:
		return -1;
}


int ssdp_setup_transfer(location *l, sky_grid *sky, sky_transfer *initsky)
{
		if (transfer_sky(l, sky, initsky)) goto estate;

		return 0;
estate:
		ssdp_free_location(l);
		return -1;
}


int ssdp_topogrid_approxhorizon(topogrid *T, int nT, int nsample,
								enum SampleType stype)
{
		int i, x, ec = 0;

		for (i=0; i < nT; ++i) {
				x = HorizonSet(T+i, nsample, stype);
				if (-1 == x) goto esobol;
				if (0 != x) ec = x;
				if (x > 0)
						printf("Approximate horizon |sample set|=%d\n", x);
		}

		return ec;
esobol:
		return -1;
}


int ssdp_topogrid_approxlaw(topogrid *T, double *q, int nq)
{
		int r = HorizonSetDstr(T, q, nq);

		if (T->horizon_sample) {
				free(T->horizon_sample);
				T->horizon_sample=NULL;
		}

		return r;
}


int ssdp_setup_grid_horizon(
		horizon *h, sky_grid *sky, topogrid *T, int nT,
		double xoff, double yoff, double zoff)
{
		// corresponds to empty locations in the uH array
		if (NULL == h)
				return 0;

		// horizon has been initialised
		if (NULL != h->zen)
				return 0;

		*h = InitHorizon(sky->Nz);
		if (ssdp_error_state) goto error;

		if (NULL == T) return 0;

		int i;
		for (i=0; i < nT; ++i)
				ComputeGridHorizon
						(h, T+i, M_PI/4.0/((double)sky->Nz),
						 xoff, yoff, zoff);
		AtanHorizon(h);
		if (ssdp_error_state) goto error;

		return 0;
error:
		return -1;
}


double ssdp_diffuse_poa(sky_grid *sky, location *l)
{
	return DiffusePlaneOfArray(sky, &(l->T));
}
double ssdp_direct_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, location *l)
{
	double g;
	g=DirectGHI(sky, l->H)*l->T.g;
	return DirectPlaneOfArray(sky, l->H, pn, M)+g;
}

double ssdp_total_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, location *l)
{
	double POA;
	POA=ssdp_diffuse_poa(sky, l);
	POA+=ssdp_direct_poa(sky, pn, M, l);
	fflush(stdout);
	return POA;
}
int ssdp_below_horizon(location *l, sky_pos p)
{
	return BelowHorizon(l->H, p);
}
/* topology routines */
topology ssdp_make_topology(double *x, double *y, double *z, int N)
{
	return MakeTopology(x,y,z, N);
}
topogrid ssdp_make_topogrid(double *z, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	return MakeTopogrid(z,x1,y1,x2,y2,Nx,Ny);
}


static inline int Ta2Tb(int ia, topogrid *Ta, topogrid *Tb)
{
		double x,y;
		int m,n;

		x = Ta->dx*(ia / Ta->Ny) + Ta->x1;
		y = Ta->dy*(ia % Ta->Ny) + Ta->y1;

		m = (int) ((x - Tb->x1)/Tb->dx);
		n = (int) ((y - Tb->y1)/Tb->dy);

		if ((m < 0) || (m >= Tb->Nx) || (n < 0) || (n >= Tb->Ny))
				return -1;

		return m*Tb->Ny + n;
}


void ssdp_min_topogrids(topogrid *T, int nT)
{
		if (nT == 1) return;

		int i, j, t, N = T[0].Nx*T[0].Ny;
		for (j=0; j < N; ++j)
				for (t=1; t < nT; ++t) {
						i = Ta2Tb(j, T, T+t);
						if (i < 0) continue;
						if (T[0].z[j] >= T[t].z[i]) continue;
						T[t].z[i] = T[0].z[j];
				}
}


topogrid ssdp_make_topogdal(double x1, double y1, double x2, double y2, char **fns, int nfns, double step, int epsg)
{
        return MakeTopoGDAL(x1, y1, x2, y2, fns, nfns, step, epsg);
}
topology ssdp_make_rand_topology(double dx, double dy, double dz, int N1, int N2)
{
	return CreateRandomTopology(dx, dy, dz, N1, N2);
}
void ssdp_free_topology(topology *T)
{
	free_topo (T);
}
void ssdp_free_topogrid(topogrid *T, int nT)
{
		int i;
		for (i=0; i < nT; ++i)
				free_topogrid(T+i);
}
double ssdp_sample_topology(double x, double y, topology *T, sky_pos *sn)
{
	return SampleTopo(x, y, T, sn);
}
double ssdp_sample_topogrid(double x, double y, topogrid *T, sky_pos *sn)
{
	return SampleTopoGrid(x, y, T, sn);
}
int ssdp_fillmissing_topogrid(topogrid *T, double na, int maxwalk)
{
		return FillMissingTopoGrid(T, na, maxwalk);
}
int ssdp_addheight_topogrid(topogrid *T,double *x, double *y, double *z, int nx, int nz)
{
		return AddHeightTopoGrid(T,x,y,z,nx,nz);
}
int ssdp_blurtopo_topogrid(topogrid *T, int size)
{
		return BlurTopoGrid(T, size);
}
// solar position
sky_pos ssdp_sunpos(time_t t, double lat, double lon, double E, double p, double T)
{
	return sunpos(t, lat, lon, E, p, T);
}
int ssdp_suntimes(time_t t, double lat, double lon, double E, double p, double T, time_t *sunrise, time_t *transit, time_t *sunset)
{
	return suntimes(t, lat, lon, E, p, T, sunrise, transit, sunset);
}


struct rtreecache* ssdp_rtreecache_init(double xy, double z)
{
		return rtreecache_init(xy, z);
}


horizon* ssdp_horizoncache_get(struct rtreecache* hc, double x, double y, double z)
{
		return horizoncache_get(hc, x, y, z);
}


void ssdp_horizoncache_free(struct rtreecache* hc)
{
		horizoncache_free(hc);
}


sky_transfer* ssdp_stcache_get(struct rtreecache* hc, double x, double y, double z)
{
		return stcache_get(hc, x, y, z);
}


void ssdp_stcache_free(struct rtreecache* hc)
{
		stcache_free(hc);
}


int ssdp_rtreecache_reset(struct rtreecache** hc)
{
		if (NULL == *hc)
				return 0;

		double xydelta = (*hc)->xydelta, zdelta = (*hc)->zdelta;
		ssdp_horizoncache_free(*hc);
		(*hc) = ssdp_rtreecache_init(xydelta, zdelta);
		if (NULL==(*hc))
				return -1;

		return 0;
}


int ssdp_write_sky(sky_grid* sky, const char* ofn, const char* dataset)
{
		return h5write_sky_grid(sky, ofn, dataset);
}


int ssdp_read_sky(sky_grid* sky, int nskies, const char* ifn, const char* dataset)
{
		int i, ret=0;

		for (i=0; i < nskies; ++i) {
				free_sky_grid(sky+i);
				ret+= h5read_sky_grid(sky+i, ifn, dataset);
		}

		return ret;
}
