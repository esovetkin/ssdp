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

typedef struct triangles {
	int i, j, k;
	double ccx, ccy, ccr;
} triangles;

typedef struct nodetree {
	int *leaves;
	double bb[4];
	struct nodetree *N1;
	struct nodetree *N2;
	struct nodetree *N3;
	struct nodetree *N4;
} nodetree;

typedef struct topology {
	double *x, *y, *z;  // 3D coordinates
	int N;		 		// number of points
	triangles *T;
	int Nt;		 		// number of triangles
	nodetree *P;
} topology;


typedef enum {QUIET, VERBOSE, VVERBOSE} VERB;
extern VERB ssdp_verbosity;

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

void ssdp_mask_horizon(sky_grid *sky, topology T, double Ox, double Oy, double Oz, sky_pos *sn);
void ssdp_unmask_horizon(sky_grid *sky);

topology ssdp_make_topology(double *x, double *y, double *z, int N);
topology ssdp_make_rand_topology(double dx, double dy, double dz, double fN, int N);
void ssdp_free_topology(topology *T);
