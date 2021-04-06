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
#include "sky_model.h"
#include "project.h"
#include "delaunay.h"
#include "ground.h"
#include "error.h"
#include "sunpos.h"
/* libssdp entry points */
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
/* AOI accptance equal for all angles */
double ssdp_diffuse_sky_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, sky_mask *mask)
{
	return DiffusePlaneOfArray(sky, pn, M, mask);
}
double ssdp_direct_sky_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, sky_mask *mask)
{
	return DirectPlaneOfArray(sky, pn, M, mask);
}
double ssdp_total_sky_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, sky_mask *mask)
{
	double POA;
	POA=DiffusePlaneOfArray(sky, pn, M, mask);
	POA+=DirectPlaneOfArray(sky, pn, M, mask);
	return POA;
}
double ssdp_groundalbedo_poa(sky_grid *sky, double albedo, sky_pos pn, AOI_Model_Data *M, sky_mask *mask)
{
	return POA_Albedo(sky, albedo, pn, M, mask);
}
double ssdp_total_poa(sky_grid *sky, double albedo, sky_pos pn, AOI_Model_Data *M, sky_mask *mask)
{
	double POA;
	POA=DiffusePlaneOfArray(sky, pn, M, mask);
	POA+=DirectPlaneOfArray(sky, pn, M, mask);
	POA+=POA_Albedo(sky, albedo, pn, M, mask);
	return POA;
}

double ssdp_diffuse_sky_horizontal(sky_grid *sky, AOI_Model_Data *M, sky_mask *mask)
{
	return DiffuseHorizontal(sky, M, mask);
}
double ssdp_direct_sky_horizontal(sky_grid *sky, AOI_Model_Data *M, sky_mask *mask)
{
	return DirectHorizontal(sky, M, mask);
}
double ssdp_total_sky_horizontal(sky_grid *sky, AOI_Model_Data *M, sky_mask *mask)
{
	double GHI;
	GHI=DiffuseHorizontal(sky, M, mask);
	GHI+=DirectHorizontal(sky, M, mask);
	return GHI;
}


/* topology routines */
sky_mask ssdp_mask_horizon(sky_grid *sky, topology *T, double Ox, double Oy, double Oz)
{
	return MaskHorizon(sky, T, Ox, Oy, Oz);
}
sky_mask ssdp_mask_horizon_z_to_ground(sky_grid *sky, topology *T, double Ox, double Oy, double deltaz, sky_pos *sn)
{
	double Oz;
	Oz=SampleTopo(Ox, Oy, T, sn)+deltaz;
	return MaskHorizon(sky, T, Ox, Oy, Oz);
}
void ssdp_unmask_horizon(sky_mask *mask)
{
	ClearSkyMask(mask);
}
void ssdp_free_sky_mask(sky_mask *mask)
{
	FreeSkyMask(mask);
}
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
