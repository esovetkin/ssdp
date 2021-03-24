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
#include "io.h"
#include "util.h"
/* libssdp entry points */

/* skydome routines */
int ssdp_find_skypatch(sky_grid sky, sky_pos p)
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

/* project routines */

double ssdp_diffuse_sky_poa(sky_grid sky, double tilt, double a, int mask)
{
	return DiffusePlaneOfArray(sky, tilt, a, mask);
}
double ssdp_direct_sky_poa(sky_grid sky, double tilt, double a, int mask)
{
	return DirectPlaneOfArray(sky, tilt, a, mask);
}
double ssdp_total_sky_poa(sky_grid sky, double tilt, double a, int mask)
{
	double POA;
	POA=DiffusePlaneOfArray(sky, tilt, a, mask);
	POA+=DirectPlaneOfArray(sky, tilt, a, mask);
	return POA;
}
double ssdp_groundalbedo_poa(sky_grid sky, double albedo, double tilt, double a, int mask)
{
	return POA_Albedo(sky, albedo, tilt, a, mask);
}
double ssdp_total_poa(sky_grid sky, double albedo, double tilt, double a, int mask)
{
	double POA;
	POA=DiffusePlaneOfArray(sky, tilt, a, mask);
	POA+=DirectPlaneOfArray(sky, tilt, a, mask);
	POA+=POA_Albedo(sky, albedo, tilt, a, mask);
	return POA;
}

void ssdp_poa_to_surface_normal(double tilt, double a, sky_pos sn, double *tilt_out, double *a_out)
{
	(*tilt_out)=tilt;
	(*a_out)=a;
	POA_to_SurfaceNormal(tilt_out, a_out, sn);
}

double ssdp_diffuse_sky_horizontal(sky_grid sky, int mask)
{
	return DiffuseHorizontal(sky, mask);
}
double ssdp_direct_sky_horizontal(sky_grid sky, int mask)
{
	return DirectHorizontal(sky, mask);
}
double ssdp_total_sky_horizontal(sky_grid sky, int mask)
{
	double GHI;
	GHI=DiffuseHorizontal(sky, mask);
	GHI+=DirectHorizontal(sky, mask);
	return GHI;
}


/* topology routines */
void ssdp_mask_horizon(sky_grid *sky, topology T, double Ox, double Oy, double Oz, sky_pos *sn)
{
	MakeHorizon(sky, T, Ox, Oy, Oz, sn);
}
void ssdp_unmask_horizon(sky_grid *sky)
{
	ClearHorizon(sky);
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

double ssdp_sample_topology(double x, double y, topology T, sky_pos *sn)
{
	return SampleTopo(x, y, T, sn);
}
