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
sky_grid ssdp_init_sky(int Nz)
{
	return InitSky(Nz);
}
void ssdp_free_sky(sky_grid *sky)
{
	free_sky_grid(sky);
}
/* skymodel routines */
void ssdp_make_uniform_sky(sky_grid *sky, sky_pos sun, double GHI, double DHI)
{
	UniformSky(sky, sun, GHI, DHI);
}
void ssdp_make_perez_all_weather_sky(sky_grid * sky, sky_pos sun, double GHI, double DHI, double dayofyear)
{
	PerezSky(sky, sun, GHI, DHI, dayofyear);
}
void ssdp_make_uniform_sky_coordinate(sky_grid *sky, time_t t, double lon, double lat, double GHI, double DHI)
{
	sky_pos sun;	
	sun=sunpos(t, lat, lon);
	UniformSky(sky, sun, GHI, DHI);
}
void ssdp_make_perez_all_weather_sky_coordinate(sky_grid * sky, time_t t, double lon, double lat, double GHI, double DHI)
{
	sky_pos sun;
	struct tm * ut;
	ut=gmtime(&t);
	if (!ut)
	{
		AddErr(GMTIMENULL); //GMTIMENULL initialized in sunpos.c
		return;
	}
	sun=sunpos(t, lat, lon);
	PerezSky(sky, sun, GHI, DHI, ut->tm_yday);
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
	horizon H;
} location;
//END_SSDP_EXPORT
void ssdp_free_location(location *l)
{
	FreeSkyTransfer(&(l->T));
	FreeHorizon(&(l->H));
}
#define ALBEPS 1e-10
location ssdp_setup_location(sky_grid *sky, topology *T, double albedo, sky_pos pn, double xoff, double yoff, double zoff, AOI_Model_Data *M)
{
	int i;
	location l;
	location l0={{NULL,0,0},{0,0,NULL}};
	// setup horizon
	l.H=InitHorizon(sky->Nz);
	l.T=InitSkyTransfer(sky->N);
	if (ssdp_error_state)
	{
		ssdp_free_location(&l);
		return l0;
	}
	if (T)
		ComputeHorizon(&l.H, T, xoff, yoff, zoff);
	AtanHorizon(&l.H);
	if (ssdp_error_state)
	{
		ssdp_free_location(&l);
		return l0;
	}
	// setup transfer
	// sky direct
	POA_Sky_Transfer(sky, &l.T, pn, M);
	if (ssdp_error_state)
	{
		ssdp_free_location(&l);
		return l0;
	}
	if (albedo>ALBEPS)
	{
		if (albedo>1)
			Print(WARNING, "Warning: albedo larger than one\n");
		l.T.g=albedo*POA_Albedo_Transfer(sky, pn, M);
		if (ssdp_error_state)
		{
			ssdp_free_location(&l);
			return l0;
		}
		for (i=0;i<l.T.N;i++)
			l.T.t[i]*=(1.0+l.T.g);
	}
	
	HorizTrans(sky, &(l.H), &(l.T), &(l.T));
	return l;
}
double ssdp_diffuse_poa(sky_grid *sky, location *l)
{
	return DiffusePlaneOfArray(sky, &(l->T));
}
double ssdp_direct_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, location *l)
{
	double g=0;
	if ((sky->sp.z<M_PI/2)&&(sky->sp.z>=0))
		g=cos(sky->sp.z)*sky->sI*l->T.g;
	return DirectPlaneOfArray(sky, &(l->H), pn, M)+g;
}

double ssdp_total_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, location *l)
{
	double POA;
	POA=ssdp_diffuse_poa(sky, l);
	//if (sky->suni>=0) /* this is an ionaccurate solution, better treat the sun separately */
	//	POA+=sky->sI*l->T.t[sky->suni];
	POA+=ssdp_direct_poa(sky, pn, M, l);
	return POA;
}
int ssdp_below_horizon(location *l, sky_pos p)
{
	return BelowHorizon(&(l->H), p);
}
/* topology routines */
topology ssdp_make_topology(double *x, double *y, double *z, int N)
{
	return MakeTopology(x,y,z, N);
}
topology ssdp_make_rand_topology(double dx, double dy, double dz, double fN, int N)
{
	return CreateRandomTopology(dx, dy, dz, fN, N);
}
void ssdp_free_topology(topology *T)
{
	free_topo (T);
}
double ssdp_sample_topology(double x, double y, topology *T, sky_pos *sn)
{
	return SampleTopo(x, y, T, sn);
}
// solar position
sky_pos ssdp_sunpos(time_t t, double lat, double lon)
{
	return sunpos(t, lat, lon);
}
