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
#include <time.h>
#include "libssdp_structs.h"

/* verbosity of the library */
typedef enum {QUIET, VERBOSE, VVERBOSE} VERB;
extern VERB ssdp_verbosity;

extern char *ssdp_error_messages[];
extern int ssdp_error_state;	/* signals error occurred, if 0 everything is OK */
void ssdp_print_error_messages();
void ssdp_reset_errors();

void ssdp_print_version(void);
/* locate the sky patch enclosing a certain point */
int ssdp_find_skypatch(sky_grid *sky, sky_pos p);
/* initialize a sky mesh */
sky_grid ssdp_init_sky(int Nz);
/* free a sky mesh */
void ssdp_free_sky(sky_grid *sky);

/* create a sky with a uniform diffuse light distribution */
void ssdp_make_uniform_sky(sky_grid *sky, sky_pos sun, double GHI, double DHI);
/* only update solar position and intensity, no diffuse light changes */
void ssdp_make_sky_sunonly(sky_grid *sky, sky_pos sun, double GHI, double DHI);
/* create a sky with the Perez all weather model */
void ssdp_make_perez_all_weather_sky(sky_grid * sky, sky_pos sun, double GHI, double DHI, double dayofyear);

void ssdp_make_perez_cumulative_sky_coordinate(sky_grid * sky, double *t, double lon, double lat, double E, double *p, double *T, double *GHI, double *DHI, int N);

/* same as above but now specifyiong time and spatial coordinates to compute the solar position and the dayofyear */
void ssdp_make_uniform_sky_coordinate(sky_grid *sky, time_t t, double lon, double lat, double E, double p, double T, double GHI, double DHI);
void ssdp_make_skysunonly_coordinate(sky_grid *sky, time_t t, double lon, double lat, double E, double p, double T, double GHI, double DHI);
void ssdp_make_perez_all_weather_sky_coordinate(sky_grid * sky, time_t t, double lon, double lat, double E, double p, double T, double GHI, double DHI);

/* projection routings for plane of array irradiance */
void ssdp_poa_to_surface_normal(sky_pos pn0, sky_pos sn, sky_pos *pn); // orient module w.r.t ground orientation

AOI_Model_Data ssdp_init_aoi_model(AOI_Model model,double nf, double nar,double *theta, double *effT, int N);

void ssdp_free_location(location *l);
int ssdp_setup_location(location *l, sky_grid *sky, topology *T, double albedo, sky_pos pn, double xoff, double yoff, double zoff, AOI_Model_Data *M);
int ssdp_setup_grid_location(location *l, sky_grid *sky, topogrid *T, double albedo, sky_pos pn, double xoff, double yoff, double zoff, AOI_Model_Data *M);
double ssdp_diffuse_poa(sky_grid *sky, location *l);
double ssdp_direct_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, location *l);
double ssdp_total_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, location *l);
int ssdp_below_horizon(location *l, sky_pos p);

// create a topology from a point cloud
topology ssdp_make_topology(double *x, double *y, double *z, int N);
topogrid ssdp_make_topogrid(double *z, double x1, double y1, double x2, double y2, int Nx, int Ny);
topogrid ssdp_make_topogdal(double x1, double y1, double x2, double y2, char **fns, int nfns, double step, int epsg);
// create random topologies for testing
topology ssdp_make_rand_topology(double dx, double dy, double dz, int N1, int N2);
void ssdp_free_topology(topology *T);
void ssdp_free_topogrid(topogrid *T);
// compute elevation (z) at any point x and y in the topology
// beware: will extrapolate from closest triangle to points outside the hull without warning
// Also computes the surface normal if you pass it a non NULL pointer to a sky_pos
// you can use this to rotate the POA using ssdp_poa_to_surface_normal(...)
double ssdp_sample_topology(double x, double y, topology *T, sky_pos *sn);
double ssdp_sample_topogrid(double x, double y, topogrid *T, sky_pos *sn);
int ssdp_fillmissing_topogrid(topogrid *T, double na);
int ssdp_addheight_topogrid(topogrid *T,double *x, double *y, double *z, int nx, int nz);

sky_pos ssdp_sunpos(time_t t, double lat, double lon, double E, double p, double T); // lat & lon in radians
int ssdp_suntimes(time_t t, double lat, double lon, double e, double p, double T, time_t *sunrise, time_t *transit, time_t *sunset);
struct horizoncache* ssdp_horizoncache_init(double xy, double z);
horizon* ssdp_horizoncache_get(struct horizoncache* hc, double x, double y, double z);
void ssdp_horizoncache_free(struct horizoncache* hc);
int ssdp_horizoncache_reset(struct horizoncache** hc);
