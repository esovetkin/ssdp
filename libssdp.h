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
typedef struct sky_pos {
	double z,a;
} sky_pos;

typedef struct hexpatch {
	double I;
	sky_pos p;
	char mask; // mark elements behind the horizon
	int NL[7]; // next level neighbors
	int PL[3]; // previous level neighbors
	int NI;    // iso level next
	int PI;    // iso level previous
} hexpatch;

typedef struct sky_grid {
	hexpatch *P;
	// get sun out of the hexpatch dome so we can separate contributions
	// if we want
	sky_pos sp;	// solar position
	double sI;	// solar intensity
	char smask;
	int N;		
	int Nz;
} sky_grid;

typedef struct topology {
	double *x, *y, *z;  // 3D coordinates
	int N;		 		// number of points
	double d;	 		// diameter of the points (all points have the same diameter)
						// you'll have to create more than one topology to have different sizes
} topology;

int ssdp_find_skypatch(sky_grid sky, sky_pos p);
sky_grid ssdp_init_sky(int Nz);
void ssdp_free_sky(sky_grid *sky);

void ssdp_make_uniform_sky(sky_grid *sky, sky_pos sun, double GHI, double DHI);
void ssdp_make_perez_all_weather_sky(sky_grid * sky, sky_pos sun, double GHI, double DHI, double dayofyear);

double ssdp_diffuse_sky_poa(sky_grid sky, double tilt, double a, int mask);
double ssdp_direct_sky_poa(sky_grid sky, double tilt, double a, int mask);
double ssdp_total_sky_poa(sky_grid sky, double tilt, double a, int mask);
double ssdp_groundalbedo_poa(sky_grid sky, double albedo, double tilt, double a, int mask);
double ssdp_total_poa(sky_grid sky, double albedo, double tilt, double a, int mask);

double ssdp_diffuse_sky_horizontal(sky_grid sky, int mask);
double ssdp_direct_sky_horizontal(sky_grid sky, int mask);
double ssdp_total_sky_horizontal(sky_grid sky, int mask);

void ssdp_mask_horizon(sky_grid *sky, topology T, double Ox, double Oy, double Oz);
void ssdp_free_topology(topology *T);
