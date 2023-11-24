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
sky_grid* ssdp_init_sky(int Nz, int nskies);
/* free a sky mesh */
void ssdp_free_sky(sky_grid *sky, int nskies);

/* create a sky with a uniform diffuse light distribution */
void ssdp_make_uniform_sky(sky_grid *sky, sky_pos sun, double GHI, double DHI);
/* only update solar position and intensity, no diffuse light changes */
void ssdp_make_sky_sunonly(sky_grid *sky, sky_pos sun, double GHI, double DHI);
/* create a sky with the Perez all weather model */
void ssdp_make_perez_all_weather_sky(sky_grid * sky, sky_pos sun, double GHI, double DHI, double dayofyear);

void ssdp_make_perez_cumulative_sky_coordinate(sky_grid * sky, double *t, double lon, double lat, double E, double *p, int np, double *T, int nT, double *GHI, double *DHI, int N);

/* same as above but now specifyiong time and spatial coordinates to compute the solar position and the dayofyear */
void ssdp_make_uniform_sky_coordinate(sky_grid *sky, time_t t, double lon, double lat, double E, double p, double T, double GHI, double DHI);
void ssdp_make_skysunonly_coordinate(sky_grid *sky, time_t t, double lon, double lat, double E, double p, double T, double GHI, double DHI);
void ssdp_make_perez_all_weather_sky_coordinate(sky_grid * sky, time_t t, double lon, double lat, double E, double p, double T, double GHI, double DHI);

/* projection routings for plane of array irradiance */
void ssdp_poa_to_surface_normal(sky_pos pn0, sky_pos sn, sky_pos *pn); // orient module w.r.t ground orientation

AOI_Model_Data ssdp_init_aoi_model(AOI_Model model,double nf, double nar,double *theta, double *effT, int N);

void ssdp_free_location(location *l);
int ssdp_init_transfer(sky_transfer* st, sky_grid *sky, double albedo, sky_pos pn, AOI_Model_Data *M);
int ssdp_setup_transfer(location *l, sky_grid *sky, sky_transfer *initsky);
int ssdp_setup_horizon(horizon *h, sky_grid *sky, topology *T, double xoff, double yoff, double zoff);
int ssdp_setup_grid_horizon(horizon *h, sky_grid *sky, topogrid *T, int nT, double xoff, double yoff, double zoff);
int ssdp_topogrid_approxhorizon(topogrid *T, int nT, int nsample);
int ssdp_topogrid_approxlaw(topogrid *T, double *q, int nq);
double ssdp_diffuse_poa(sky_grid *sky, location *l);
double ssdp_direct_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, location *l);
double ssdp_total_poa(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, location *l);
int ssdp_below_horizon(location *l, sky_pos p);

// create a topology from a point cloud
topology ssdp_make_topology(double *x, double *y, double *z, int N);
topogrid ssdp_make_topogrid(double *z, double x1, double y1, double x2, double y2, int Nx, int Ny);
void ssdp_min_topogrids(topogrid *T, int nT);
topogrid ssdp_make_topogdal(double x1, double y1, double x2, double y2, char **fns, int nfns, double step, int epsg);
// create random topologies for testing
topology ssdp_make_rand_topology(double dx, double dy, double dz, int N1, int N2);
void ssdp_free_topology(topology *T);
void ssdp_free_topogrid(topogrid *T, int nT);
// compute elevation (z) at any point x and y in the topology
// beware: will extrapolate from closest triangle to points outside the hull without warning
// Also computes the surface normal if you pass it a non NULL pointer to a sky_pos
// you can use this to rotate the POA using ssdp_poa_to_surface_normal(...)
double ssdp_sample_topology(double x, double y, topology *T, sky_pos *sn);
double ssdp_sample_topogrid(double x, double y, topogrid *T, sky_pos *sn);
int ssdp_fillmissing_topogrid(topogrid *T, double na, int maxwalk);
int ssdp_addheight_topogrid(topogrid *T,double *x, double *y, double *z, int nx, int nz);
int ssdp_blurtopo_topogrid(topogrid *T,int size);

sky_pos ssdp_sunpos(time_t t, double lat, double lon, double E, double p, double T); // lat & lon in radians
int ssdp_suntimes(time_t t, double lat, double lon, double e, double p, double T, time_t *sunrise, time_t *transit, time_t *sunset);

struct rtreecache* ssdp_rtreecache_init(double xy, double z);
int ssdp_rtreecache_reset(struct rtreecache** hc);

horizon* ssdp_horizoncache_get(struct rtreecache* hc, double x, double y, double z);
void ssdp_horizoncache_free(struct rtreecache* hc);

sky_transfer* ssdp_stcache_get(struct rtreecache* hc, double x, double y, double z);
void ssdp_stcache_free(struct rtreecache* hc);

int ssdp_write_sky(sky_grid* sky, const char* ofn, const char *dataset);
int ssdp_read_sky(sky_grid* sky, int nskies, const char* ifn, const char* dataset);
