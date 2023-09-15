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

static int transfer_sky(
		location *l, sky_grid *sky,
		double albedo, sky_pos pn, AOI_Model_Data *M)
{
		int i;

		// setup transfer, sky direct
		l->T=InitSkyTransfer(sky->N);
		if (ssdp_error_state) goto einit;
		POA_Sky_Transfer(sky, &(l->T), pn, M);
		if (ssdp_error_state) goto error;

		if (albedo>ALBEPS) {
				l->T.g=albedo*POA_Albedo_Transfer(sky, pn, M);
				if (ssdp_error_state) goto error;

				for (i=0;i<l->T.N;i++)
						l->T.t[i]*=(1.0+l->T.g);
		}

		HorizTrans(sky, l->H, &(l->T), &(l->T));
		if (ssdp_error_state) {
				// TODO Error, not warning
				Print(WARNING, "Warning: Horizon does not match the sky\n");
				goto error;
		}


		l->difftrans=0;
		for (i=0;i<l->T.N;i++)
				l->difftrans+=sky->sa[i]*l->T.t[i];
		l->difftrans/=sky->icosz;

		return 0;
error:
		FreeSkyTransfer(&(l->T));
einit:
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

int ssdp_setup_transfer(
		location *l, sky_grid *sky, double albedo, sky_pos pn,
		AOI_Model_Data *M)
{
		if (transfer_sky(l, sky, albedo, pn, M)) goto estate;

		return 0;
estate:
		ssdp_free_location(l);
		return -1;
}


int ssdp_setup_grid_horizon(
		horizon *h, sky_grid *sky, topogrid *T,
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

		ComputeGridHorizon
				(h, T, M_PI/4.0/((double)sky->Nz),
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
void ssdp_free_topogrid(topogrid *T)
{
	free_topogrid (T);
}
double ssdp_sample_topology(double x, double y, topology *T, sky_pos *sn)
{
	return SampleTopo(x, y, T, sn);
}
double ssdp_sample_topogrid(double x, double y, topogrid *T, sky_pos *sn)
{
	return SampleTopoGrid(x, y, T, sn);
}
int ssdp_fillmissing_topogrid(topogrid *T, double na)
{
		return FillMissingTopoGrid(T, na);
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

struct horizoncache* ssdp_horizoncache_init(double xy, double z)
{
		return horizoncache_init(xy, z);
}
horizon* ssdp_horizoncache_get(struct horizoncache* hc, double x, double y, double z)
{
		return horizoncache_get(hc, x, y, z);
}
void ssdp_horizoncache_free(struct horizoncache* hc)
{
		horizoncache_free(hc);
}
int ssdp_horizoncache_reset(struct horizoncache** hc)
{
		if (NULL == *hc)
				return 0;

		double xydelta = (*hc)->xydelta, zdelta = (*hc)->zdelta;
		ssdp_horizoncache_free(*hc);
		(*hc) = ssdp_horizoncache_init(xydelta, zdelta);
		if (NULL==(*hc))
				return -1;

		return 0;
}
