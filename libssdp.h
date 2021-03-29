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

/* spherical coordinate system (only direction, no radius, hardly ever need that) */
typedef struct sky_pos {
	double z,a;
} sky_pos;

/* hexagonal patch of sky */
typedef struct hexpatch {
	double I;	// light intensity
	sky_pos p;	// sky coordinate
	char mask;	// mark elements behind the horizon
	int NL[7];	// next level neighbors (at larger zenith angle)
	int PL[3];	// previous level neighbors (at smaller senith angle)
	int NI;		// iso level next (same zenith angle, larger azimuth)
	int PI;		// iso level previous (same zenith angle, smaller azimuth)
} hexpatch;


/* sky mesh */
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

/* some structures for the topology */
/* triangulation data */
typedef struct triangles {
	int i, j, k;
	double ccx, ccy;
} triangles;
/* quaternary search tree */
typedef struct nodetree {
	int *leafs;
	double bb[4];
	struct nodetree *N1;
	struct nodetree *N2;
	struct nodetree *N3;
	struct nodetree *N4;
} nodetree;

/* topology structure including triangulation data */
typedef struct topology {
	double *x, *y, *z;  // 3D coordinates
	int N;		 		// number of points
	triangles *T;
	int Nt;		 		// number of triangles
	nodetree *P;
} topology;


/* verbosity of the library */
typedef enum {QUIET, VERBOSE, VVERBOSE} VERB;
extern VERB ssdp_verbosity;

extern char *ssdp_error_messages[];
extern int ssdp_error_state;	/* signals error occurred, if 0 everything is OK */
void ssdp_print_error_messages();
void ssdp_reset_errors();

/* locate the sky patch enclosing a certain point */
int ssdp_find_skypatch(sky_grid *sky, sky_pos p);
/* initialize a sky mesh */
sky_grid ssdp_init_sky(int Nz);
/* free a sky mesh */
void ssdp_free_sky(sky_grid *sky);

/* create a sky with a uniform diffuse light distribution */
void ssdp_make_uniform_sky(sky_grid *sky, sky_pos sun, double GHI, double DHI);
/* create a sky with the Perez all weather model */
void ssdp_make_perez_all_weather_sky(sky_grid * sky, sky_pos sun, double GHI, double DHI, double dayofyear);

/* projection routings for plane of array irradiance */
void ssdp_poa_to_surface_normal(double tilt, double a, sky_pos sn, double *tilt_out, double *a_out); // orient module w.r.t grounbd orientation

/* AOI accptance equal for all angles */
double ssdp_diffuse_sky_poa(sky_grid * sky, double tilt, double a, int mask);					// diffuse contribution
double ssdp_direct_sky_poa(sky_grid * sky, double tilt, double a, int mask);					// direct contribution
double ssdp_total_sky_poa(sky_grid * sky, double tilt, double a, int mask);					// all sky contributions together
double ssdp_groundalbedo_poa(sky_grid * sky, double albedo, double tilt, double a, int mask);	// ground albedo contribution (with crude assumptions)
double ssdp_total_poa(sky_grid * sky, double albedo, double tilt, double a, int mask);		// sky+ground
double ssdp_diffuse_sky_horizontal(sky_grid * sky, int mask);
double ssdp_direct_sky_horizontal(sky_grid * sky, int mask);
double ssdp_total_sky_horizontal(sky_grid * sky, int mask);

/* AOI accptance: assume plain front cover */
double ssdp_diffuse_sky_poa_effective(sky_grid *sky, double tilt, double a, double nf, int mask);
double ssdp_direct_sky_poa_effective(sky_grid *sky, double tilt, double a, double nf, int mask);
double ssdp_total_sky_poa_effective(sky_grid *sky, double tilt, double a, double nf, int mask);
double ssdp_groundalbedo_poa_effective(sky_grid *sky, double albedo, double tilt, double a, double nf, int mask);
double ssdp_total_poa_effective(sky_grid *sky, double albedo, double tilt, double a, double nf, int mask);
double ssdp_diffuse_sky_horizontal_effective(sky_grid *sky, double nf, int mask);
double ssdp_direct_sky_horizontal_effective(sky_grid *sky, double nf, int mask);
double ssdp_total_sky_horizontal_effective(sky_grid *sky, double nf, int mask);

/* AOI acceptance: assume front cover with AR coating */
double ssdp_diffuse_sky_poa_effective_ar(sky_grid *sky, double tilt, double a, double nf, double nar, int mask);
double ssdp_direct_sky_poa_effective_ar(sky_grid *sky, double tilt, double a, double nf, double nar, int mask);
double ssdp_total_sky_poa_effective_ar(sky_grid *sky, double tilt, double a, double nf, double nar, int mask);
double ssdp_groundalbedo_poa_effective_ar(sky_grid *sky, double albedo, double tilt, double a, double nf, double nar, int mask);
double ssdp_total_poa_effective_ar(sky_grid *sky, double albedo, double tilt, double a, double nf, double nar, int mask);
double ssdp_diffuse_sky_horizontal_effective_ar(sky_grid *sky, double nf, double nar, int mask);
double ssdp_direct_sky_horizontal_effective_ar(sky_grid *sky, double nf, double nar, int mask);
double ssdp_total_sky_horizontal_effective_ar(sky_grid *sky, double nf, double nar, int mask);


/* compute the horizon */
void ssdp_mask_horizon(sky_grid *sky, topology *T, double Ox, double Oy, double Oz); // absolute x,y,z coordinates
void ssdp_mask_horizon_z_to_ground(sky_grid *sky, topology *T, double Ox, double Oy, double deltaz, sky_pos *sn); 
// z w.r.t. ground, also sets sn, the local surface normal. Can be used to adapt module orientation (ssdp_poa_to_surface_normal)
void ssdp_unmask_horizon(sky_grid *sky); // clear a horizon from a sky dome


// create a topology from a point cloud
topology ssdp_make_topology(double *x, double *y, double *z, int N);
// create random topologies for testing
topology ssdp_make_rand_topology(double dx, double dy, double dz, double fN, int N);
void ssdp_free_topology(topology *T);
// compute elevation (z) at any point x and y in the topology
// beware: will extrapolate from closest triangle to points outside the hull without warning
// Also computes the surface normal if you pass it a non NULL pointer to a sky_pos
// you can use this to rotate the POA using ssdp_poa_to_surface_normal(...)
double ssdp_sample_topology(double x, double y, topology T, sky_pos *sn);
