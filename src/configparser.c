#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#ifdef OPENMP
  #include <omp.h>
#endif
/* local includes */
#include "libssdp.h"
#include "lio.h"
#include "util.h"
#include "variables.h"
#include "parser.h"
#include "parserutil.h"
#include "epsg.h"
#include "coordinates.h"

#ifdef RUNMEMTEST
#include "random_fail_malloc.h"
#define malloc(x) random_fail_malloc(x)
#endif

/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG init_sim_config InitConfig "C=<out-config>"
DESCRIPTION Create a configuration variable.
OUTPUT C output configuration variable
END_DESCRIPTION
*/
void InitConfig(char *in)
{
	simulation_config C;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (GetArg(in, "C", word))
	{
		C=InitConf();
		printf("Defining simulation configuration \"%s\"\n", word);
		if (AddSimConf(word, C)) // note InitConf does not allocate anything, no need to free
			free(word);
		return;
	}
	free(word);
}
/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_coord ConfigCoord "C=<out-config> lat=<float-value> lon=<float-value> E=<float-value>"
DESCRIPTION Setup the coordinate in the configuration variable.
ARGUMENT lat latitude (in radians)
ARGUMENT lon longitude (in radians)
ARGUMENT E elevation (in m, default 0)
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigCoord (char *in)
{
	simulation_config *C;
    double l;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "lon", word, &l))
	{
		free(word);
		return;
	}
	C->lon=l;
	printf("set longitude to %e degrees (%e rad)\n", rad2deg(C->lon), C->lon);
	if (FetchFloat(in, "lat", word, &l))
	{
		free(word);
		return;
	}
	C->lat=l;
	printf("set latitude to %e degrees (%e rad)\n", rad2deg(C->lat), C->lat);
	
	if (FetchOptFloat(in, "E", word, &l))
        l = 0;

    C->E=l;
	printf("set elevation to %e m\n", C->E);	
	free(word);
}

/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_aoi ConfigAOI "C=<out-config> model=<none/front-cover/anti-reflect/user> [nf=<float-value>] [nar=<float-value>] [file=<in-file>]"
DESCRIPTION Setup the Angle of Incidence model to model angular dependent reflection. The model can be one of "none" (no angularly dependent reflection), "front-cover" (simple refractive index), "anti-relect" (two layer front), and "user" (tabular data).
ARGUMENT model string to identify which model to use
ARGUMENT nf front-cover refractive index
ARGUMENT nar refective index of an antireflection coating on the front cover
ARGUMENT file load tabular data of angular dependent reflection (do not use, it is half implemented)
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigAOI(char *in)
{
	simulation_config *C;
	AOI_Model M;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
		
	if (!GetArg(in, "model", word))
	{
		free(word);
		return;
	}
	if (strncmp(word,"none",12)==0)
		M=AOI_NONE;
	else if (strncmp(word,"front-cover",12)==0)
		M=AOI_GLASS;
	else if (strncmp(word,"anti-reflect",12)==0)
		M=AOI_GLASS_AR;
	else if (strncmp(word,"user",12)==0)
		M=AOI_USER;
	else
	{
		Warning("Unknown AOI model %s\n",word); 
		free(word);
		return;
	}
	
	switch (M)
	{
		case AOI_GLASS:
		{
			if (FetchFloat(in, "nf", word, &(C->M.ng)))
			{
				free(word);
				return;
			}
			printf("Setting AOI model to \"front-cover\" with a refractive index of %e\n", C->M.ng);
			C->M.M=M;
			return;
		}
		case AOI_GLASS_AR:
		{
			if (FetchFloat(in, "nf", word, &(C->M.ng)))
			{
				free(word);
				return;
			}
			if (FetchFloat(in, "nar", word, &(C->M.nar)))
			{
				free(word);
				return;
			}
			printf("Setting AOI model to \"anti-reflect\" with a cover refractive index of %e\n", C->M.ng);
			printf("and a anti-reflecxtion refrective index of %e\n", C->M.nar);
			C->M.M=M;
			return;
		}
		case AOI_USER:
		{
			if (!GetArg(in, "file", word))
			{
				free(word);
				return;
			}
			// read and allocate theta-effT data
			Warning("Not implemented yet\n"); 
			C->M.M=M;
			return;
		}
		default:
			break;		
	}
	free(word);
}

void FreeConfigMask(simulation_config *C)
{
	int i;
	if (C->L)
	{
		for (i=0;i<C->Nl;i++)
			ssdp_free_location(&(C->L[i]));
		free(C->L);
	}
	C->L=NULL;
	if (C->uH) {
			free(C->uH);
	}
	C->uH=NULL;
	if (C->uHi) {
			free(C->uHi);
	}
	C->uHi=NULL;
}

void FreeConfigLocation(simulation_config *C)
{
	if (C->loc_init)
	{
		if (C->x)
			free(C->x);
		if (C->y)
			free(C->y);
		if (C->z)
			free(C->z);
		if (C->o)
			free(C->o);
	}
	FreeConfigMask(C);
	C->Nl=0;
	C->loc_init=0;
}


static void set_uH(simulation_config *C, int i)
{
		// basic hashmap. number of entries in hash equals to its
		// length, so we expect a lot of collisions. however, we do
		// not need to search in hash, only set an element once. If
		// collisions become critical, C->uH length can always be
		// extended, and hashv value can be randomised

		uintptr_t j, k;
		j = k = (uintptr_t) C->L[i].H;
		while (NULL != C->uH[j % C->Nl]) {
				// horizon has been visited
				if (k == (uintptr_t) C->uH[j % C->Nl])
						return;

				++j;
		}

		C->uH[j % C->Nl] = C->L[i].H;
		C->uHi[j % C->Nl] = i;
}


static int init_hcache(simulation_config *C)
{
		FreeConfigMask(C); // make sure we are clear to allocate new memory
		if (NULL==(C->L=calloc(C->Nl,sizeof(*(C->L))))) goto ecl;
		if (NULL==(C->uH=calloc(C->Nl,sizeof(*(C->uH))))) goto euH;
		if (NULL==(C->uHi=calloc(C->Nl,sizeof(*(C->uHi))))) goto euHi;

		TIC();
		int i, hits=0;
		// assume C->hcache initialised. associate location with the
		// allocated space for that in the hcache.
		for (i=0; i < C->Nl; ++i) {
				C->L[i].H = ssdp_horizoncache_get(C->hcache, C->x[i], C->y[i], C->z[i]);
				set_uH(C, i);
				if (NULL == C->L[i].H) goto ehcache;
				if (NULL != C->L[i].H->zen) ++hits;
		}
		printf("Checked in locations cache in %g s, hits/misses: %d/%d\n", TOC(), hits, C->Nl-hits);

		return 0;
		// the hcache keeps track of allocated horizon. even if we
		// fail and end up here, just need to cleanup the
		// C->L. ssdp_horizoncache_free at some point will cleanup
		// everything.
ehcache:
		free(C->uHi); C->uHi = NULL;
euHi:
		free(C->uH); C->uH = NULL;
euH:
		free(C->L); C->L = NULL;
ecl:
		return -1;
}


static int init_transfer(simulation_config *C)
{
		double dt;
		int i, pco=0;

		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->Nl; ++i) {
						ssdp_setup_transfer(&(C->L[i]), &(C->S),
											C->albedo, C->o[i], &(C->M));
						ProgressBar((100*(i+1))/C->Nl, &pco, ProgressLen, ProgressTics);
				}
		}
		ProgressBar(100, &pco, ProgressLen, ProgressTics);
		if (ssdp_error_state) goto error;
		dt=TOC();
		printf("Initialised sky transfer in %g s (%e s/locations)\n", dt, dt/((double)C->Nl));
		return 0;
error:
		return -1;
}


void InitConfigMask(simulation_config *C)
{
		double dt;
		int i, pco = 0;

		if (! ((C->sky_init)&&(C->topo_init)&&(C->loc_init))) goto egtl;
		if (init_hcache(C)) goto einitL;

		printf("Tracing %d locations\n", C->Nl);
		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->Nl; ++i) {
						ssdp_setup_horizon(C->uH[i], &(C->S), &(C->T),
										   C->x[C->uHi[i]],
										   C->y[C->uHi[i]],
										   C->z[C->uHi[i]]);
						ProgressBar((100*(i+1))/C->Nl, &pco, ProgressLen, ProgressTics);
				}
		}
		ProgressBar(100, &pco, ProgressLen, ProgressTics);
		dt = TOC();
		printf("Traced %d horizons in %g s (%e s/horizons)\n", C->Nl, dt, dt/((double)C->Nl));

		if (init_transfer(C)) goto estate;
		return;
estate:
		FreeConfigMask(C); // make sure we are clear to allocate new memory
einitL:
egtl:
		return;
}

void InitConfigGridMask(simulation_config *C)
{
		double dt;
		int i, pco=0;

		if (! ((C->sky_init)&&(C->grid_init)&&(C->loc_init))) goto esgl;
		if (init_hcache(C)) goto einitL;

		// process all necessary horizons separately
		printf("Tracing %d locations\n", C->Nl);
		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->Nl; ++i) {
						ssdp_setup_grid_horizon(
								C->uH[i], &(C->S), &(C->Tx),
								C->x[C->uHi[i]],
								C->y[C->uHi[i]],
								C->z[C->uHi[i]]);
						ProgressBar((100*(i+1))/C->Nl, &pco, ProgressLen, ProgressTics);
				}
		}
		ProgressBar(100, &pco, ProgressLen, ProgressTics);
		dt = TOC();
		printf("Traced %d horizons in %g s (%e s/horizons)\n", C->Nl, dt, dt/((double)C->Nl));

		if (init_transfer(C)) goto estate;
		return;
estate:
		FreeConfigMask(C); // make sure we are clear to allocate new memory
einitL:
esgl:
		return;
}

void InitConfigMaskNoH(simulation_config *C) // same as above but without horizon calculation
{
		double dt;
		int i, pco=0;

		if (! ((C->sky_init)&&(C->loc_init))) goto esl;
		if (init_hcache(C)) goto einitL;

		// some horizon needs to be initialised
		printf("Initialising %d horizons\n", C->Nl);
		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->Nl; ++i) {
						ssdp_setup_horizon(C->uH[i], &(C->S), NULL,
										   C->x[C->uHi[i]],
										   C->y[C->uHi[i]],
										   C->z[C->uHi[i]]);
						ProgressBar((100*(i+1))/C->Nl, &pco, ProgressLen, ProgressTics);
				}
		}
		ProgressBar(100, &pco, ProgressLen, ProgressTics);
		dt = TOC();
		printf("Initialised %d horizons in %g s (%e s/horizons)\n", C->Nl, dt, dt/((double)C->Nl));

		if (init_transfer(C)) goto estate;
		return;
estate:
		FreeConfigMask(C); // make sure we are clear to allocate new memory
einitL:
esl:
		return;
}


/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_sky ConfigSKY "C=<out-config> N=<int-value>"
DESCRIPTION Setup the sky. This commands allocates space and initializes the sky data.
ARGUMENT N The number of zenith discretizations. The total number of sky patches equals Ntotal(N)=3N(N-1)+1, e.g. with Ntotal(7)=127
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigSKY(char *in)
{
	simulation_config *C;
	char *word;
	int N;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
		
	if (FetchInt(in, "N", word, &N))
	{
		free(word);
		return;
	}
	if (N>1)
	{
		if (C->sky_init)
			ssdp_free_sky(&C->S);
		else
			C->sky_init=1;
		C->S=ssdp_init_sky(N);
		ssdp_horizoncache_reset(&(C->hcache));
		InitConfigMask(C);		
		printf("Configuring sky with %d zenith discretizations\n", N);
	}
	else
		Warning("Number of zenith discretizations must me larger than 1\n");
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_sky(&C->S);
		C->sky_init=0;
		ssdp_reset_errors();
	}
	free(word);
	return;
}
/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_topology ConfigTOPO "C=<out-config> x=<in-array> y=<in-array> z=<in-array>"
DESCRIPTION Setup the topography. Load the x, y, and z data of the unstructured topography mesh into the configuration data.
ARGUMENT x x coordinates
ARGUMENT y y coordinates
ARGUMENT z z coordinates
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigTOPO (char *in)
{
	simulation_config *C;
	array *x, *y, *z;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
		
	if (FetchArray(in, "x", word, &x))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "y", word, &y))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "z", word, &z))
	{
		free(word);
		return;
	}
	free(word);
	if ((x->N!=y->N)||(x->N!=z->N))
	{
		Warning("x- y- and z-arrays must be of equal length\n"); 
		return;
	}
	
	if (C->topo_init)
	{
		ssdp_free_topology(&C->T);
	}
	else
		C->topo_init=1;
	printf("Configuring topology with %d points\n", x->N);
	C->T=ssdp_make_topology(x->D, y->D, z->D, x->N);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_topology(&C->T);
		C->topo_init=0;
		ssdp_reset_errors();
	}
	ssdp_horizoncache_reset(&(C->hcache));
	InitConfigMask(C);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_topology(&C->T);
		C->topo_init=0;
		ssdp_reset_errors();
	}
	return;
}
/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_topogrid ConfigTOPOGrid "C=<out-config> z=<in-array> Nx=<int-value> Ny=<int-value> x1=<float-value> y1=<float-value> x2=<float-value> y2=<float-value>"
DESCRIPTION Setup the topogrid. Load the z data (column major, from the south-west corner to the north-east corner).
ARGUMENT z z coordinates
ARGUMENT Nx number of x steps.
ARGUMENT Ny number of y steps.
ARGUMENT x1 x-coordinate of the south-west corner
ARGUMENT y1 y-coordinate of the south-west corner
ARGUMENT x2 x-coordinate of the north-east corner
ARGUMENT y2 y-coordinate of the north-east corner
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigTOPOGrid (char *in)
{
	simulation_config *C;
	array *z;
	double x1, y1, x2, y2;
	int Nx, Ny;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
		
	if (FetchFloat(in, "x1", word, &x1))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "y1", word, &y1))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "x2", word, &x2))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "y2", word, &y2))
	{
		free(word);
		return;
	}
	if (FetchInt(in, "Nx", word, &Nx))
	{
		free(word);
		return;
	}
	if (FetchInt(in, "Ny", word, &Ny))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "z", word, &z))
	{
		free(word);
		return;
	}
	free(word);
	if (z->N!=Nx*Ny)
	{
		Warning("Number of steps in x- and y- directions do not match the number of elements in the z array (%d != %d x %d)\n", z->N, Nx, Ny); 
		return;
	}
	
	if (C->grid_init)
	{
		ssdp_free_topogrid(&C->Tx);
	}
	else
		C->grid_init=1;
	printf("Configuring topogrid with %d points\n", z->N);
	C->Tx=ssdp_make_topogrid(z->D, x1, y1, x2, y2, Nx, Ny);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_topogrid(&C->Tx);
		C->grid_init=0;
		ssdp_reset_errors();
	}
	ssdp_horizoncache_reset(&(C->hcache));
	InitConfigGridMask(C);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_topogrid(&C->Tx);
		C->grid_init=0;
		ssdp_reset_errors();
	}
	return;
}


/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_topogdal ConfigTOPOGDAL "C=<out-config> lat1=<float-value> lon1=<float-value> lat2=<float-value> lon2=<float-value> step=<float-value> f0=<file-str> f1=<file-str> .. fN=<file-str> [flist=<file-str>] [epsg=<int-value>]"
DESCRIPTION Setup the topogrid using Georasters
ARGUMENT lat1 latitude of the south-west corner (in WGS84, epsg:4326)
ARGUMENT lon1 longitude of the south-west corner (in WGS84, epsg:4326)
ARGUMENT lat2 latitude of the north-east corner (in WGS84, epsg:4326)
ARGUMENT lon2 longitude of the north-east corner (in WGS84, epsg:4326)
ARGUMENT fi input i-th raster file
ARGUMENT step step at which rasters are sampled (units of step depend on epsg)
ARGUMENT flist optional path to a file containing paths to rasters
ARGUMENT epsg optional coordinate system where to resample rasters. By default an approritate UTM system is used, which is determinted with the centre of the box
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigTOPOGDAL (char *in)
{
        simulation_config *C;
        double x1, y1, x2, y2, step;
        int i = 0, epsg = -1;
        char *word;
        word=malloc((strlen(in)+1)*sizeof(*word));
        if (NULL == word) goto eword;

        if (FetchConfig(in, "C", word, &C)) goto epars;
        if (FetchFloat(in, "lat1", word, &x1)) goto epars;
        if (FetchFloat(in, "lon1", word, &y1)) goto epars;
        if (FetchFloat(in, "lat2", word, &x2)) goto epars;
        if (FetchFloat(in, "lon2", word, &y2)) goto epars;
        if (FetchFloat(in, "step", word, &step)) goto epars;

        struct cvec *fns = cvec_init(4);
        if (NULL == fns) goto ecvec;
        while(GetNumOption(in, "f", i, word)) {
                if (cvec_push(fns, word)) goto efnspush;
                word=malloc((strlen(in)+1)*sizeof(*word));
                if (NULL == word) goto efnspush;
                ++i;
        }

        if (GetOption(in, "flist", word))
                if (read_filelist(word, fns)) {
                        Warning("Error: Failed to read path list from the file!\n");
                        goto efnlist;
                }

        if (FetchOptInt(in, "epsg", word, &epsg))
                epsg = determine_utm((x1+x2)/2, (y1+y2)/2);

        printf("Sampling with step=%.3f epsg=%d\n"
               "\tbox: %.5f %.5f %.5f %.5f\n"
               "\tfiles:\n",
               step, epsg, x1, y1, x2, y2);
        for (i=0; i < fns->n; ++i)
                printf("\t%s\n", fns->s[i]);

        if (C->grid_init)
                ssdp_free_topogrid(&C->Tx);
        else
                C->grid_init = 1;

        TIC();
        C->Tx = ssdp_make_topogdal(x1, y1, x2, y2, fns->s, fns->n, step, epsg);
        if (ssdp_error_state) goto emaketopogdal;
        if (ssdp_horizoncache_reset(&(C->hcache))) goto emaketopogdal;
        InitConfigGridMask(C);
        if (ssdp_error_state) goto emaketopogdal;
        printf("Initialised topogrid in %g s\n", TOC());

        cvec_free(fns);
        free(word);
        return;
emaketopogdal:
        ssdp_print_error_messages();
        ssdp_free_topogrid(&C->Tx);
        C->grid_init=0;
        ssdp_reset_errors();
efnlist:
efnspush:
        cvec_free(fns);
ecvec:
epars:
        free(word);
eword:
        Warning("Error: config_topogdal failed!\n");
        return;
}

/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_cleanlocations ConfigCleanLocations  "C=<out-config>"
DESCRIPTION Clean stored locations from cache
OUTPUT C output configuration variable
END_DESCRIPTION
*/
void ConfigCleanLocations(char *in)
{
		simulation_config *C;
		char *word;
		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (FetchConfig(in, "C", word, &C)) goto eargs;

		ssdp_horizoncache_free(C->hcache);
		C->hcache = NULL;

		free(word);
		return;
eargs:
		free(word);
eword:
        Warning("Error: config_cleanlocations failed!\n");
		return;
}

/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_locations ConfigLoc "C=<out-config> x=<in-array> y=<in-array> z=<in-array> azimuth=<in-array> zenith=<in-array> [type=<topology/topogrid>] [albedo=<float-value>] [xydelta=<float-value>] [zdelta=<float-value]"
DESCRIPTION Setup the topography. Load the x, y, and z data of the unstructured topography mesh into the configuration data.
ARGUMENT x x coordinates
ARGUMENT y y coordinates
ARGUMENT z z coordinates
ARGUMENT type Optional argument to select the unstructured topology mesh or the structured topogrid. Valid values are \"topology\" or \"topogrid\" (default topology unless only a topogrid is defined)
ARGUMENT azimuth azimuth angle of tilted surface
ARGUMENT zenith zenith angle of tilted surface
ARGUMENT albedo optionally provide an albedo value between 0-1
ARGUMENT xydelta,zdelta the coordinates within xydelta in xy plane and zdelta within z direction are considered the same (default: 0.05)
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigLoc(char *in)
{
		simulation_config *C;
		double xydelta, zdelta;
		array *x, *y, *z, *az, *ze;
		int N, i;
		char type='t';
		char *word;
		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (FetchArray(in, "x", word, &x)) goto eargs;
		if (FetchArray(in, "y", word, &y)) goto eargs;
		if (FetchArray(in, "z", word, &z)) goto eargs;
		if (FetchArray(in, "azimuth", word, &az)) goto eargs;
		if (FetchArray(in, "zenith", word, &ze)) goto eargs;
		if (GetOption(in, "type", word)) {
				if (strncmp(word,"topology",10)==0) {
						if (C->topo_init==0) {
								Warning("No topology available\n");
								goto eargs;
						}
				}
				if (strncmp(word,"topogrid",10)==0) {
						if (C->grid_init==0) {
								Warning("No topogrid available\n");
								goto eargs;
						}
						type='g';
				}
		}
		if (C->topo_init==0) {
				type='g';
				if (C->grid_init==0)
						Warning("No topology or topogrid available\n");
		}
		if (FetchOptFloat(in, "albedo", word, &(C->albedo))) C->albedo=0.0;
		if (FetchOptFloat(in, "xydelta", word, &(xydelta))) xydelta=0.05;
		if (FetchOptFloat(in, "zdelta", word, &(zdelta))) zdelta=0.05;

		N = check_shapes(5,(array*[]){x,y,z,az,ze});

		FreeConfigLocation(C);
		printf("Configuring locations with %d points\n", x->N);

		C->Nl=N;
		if (NULL==(C->x=malloc(C->Nl*sizeof(*(C->x))))) goto ex;
		if (NULL==(C->y=malloc(C->Nl*sizeof(*(C->y))))) goto ey;
		if (NULL==(C->z=malloc(C->Nl*sizeof(*(C->z))))) goto ez;
		if (NULL==(C->o=malloc(C->Nl*sizeof(*(C->o))))) goto eo;

		for (i=0; i<C->Nl; ++i) {
				C->x[i]=AT(x,i);
				C->y[i]=AT(y,i);
				C->z[i]=AT(z,i);
				C->o[i].a=AT(az,i);
				C->o[i].z=AT(ze,i);
		}
		C->loc_init=1;

		if (NULL==C->hcache)
				if (NULL==(C->hcache=ssdp_horizoncache_init(xydelta, zdelta)))
						goto ehcache;
		C->hcache->xydelta = xydelta / 2.0;
		C->hcache->zdelta = zdelta / 2.0;

		if ((C->topo_init==1)&&(type=='t'))
				InitConfigMask(C);
		if ((C->grid_init==1)&&(type=='g'))
				InitConfigGridMask(C);

		if (ssdp_error_state)
		{
				ssdp_print_error_messages();
				FreeConfigLocation(C);
				C->loc_init=0;
				ssdp_reset_errors();
				goto elocs;
		}

		free(word);
		return;
elocs:
ehcache:
		free(C->x);
ex:
		free(C->y);
ey:
		free(C->z);
ez:
		free(C->o);
eo:
eargs:
		free(word);
eword:
		Warning("config_locations: failed!\n");
		return;
}


/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG rand_topology RandTOPO "C=<out-config> dx=<float-value> dy=<float-value> dz=<float-value> N1=<int-value> N2=<int-value>"
DESCRIPTION Create a random topography in the configuration variable.
ARGUMENT dx Range in x direction
ARGUMENT dy Range in y-direction
ARGUMENT dz Range in z-direction
ARGUMENT N1 Initial random points
ARGUMENT N2 Resampled number of points (N2>N1)
OUTPUT C configuration variable
END_DESCRIPTION
*/
void RandTOPO (char *in)
{
	simulation_config *C;
	double dx, dy, dz;
	int N1, N2;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
		
	if (FetchFloat(in, "dx", word, &dx))
	{
		free(word);
		return;
	}
		
	if (FetchFloat(in, "dy", word, &dy))
	{
		free(word);
		return;
	}
		
	if (FetchFloat(in, "dz", word, &dz))
	{
		free(word);
		return;
	}
		
	if (FetchInt(in, "N1", word, &N1))
	{
		free(word);
		return;
	}
	if (FetchInt(in, "N2", word, &N2))
	{
		free(word);
		return;
	}
	free(word);
	
	if ((N1<3)||(N2<3))
	{
		Warning("N1 and N2 must be larger than 3 in rand_topology\n"); 
		return;
	}
	if ((dx<1e-10)||(dy<1e-10)||(dz<1e-10))
	{
		Warning("Please provide valid positive ranges for the random topology\n"); 
		return;
	}
	
	if (C->topo_init)
	{
		ssdp_free_topology(&C->T);
	}
	else
		C->topo_init=1;
	printf("Configuring topology with %d points\n", N2);
	C->T=ssdp_make_rand_topology(dx, dy, dz, N1, N2);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_topology(&C->T);
		C->topo_init=0;
		ssdp_reset_errors();
	}
	InitConfigMask(C);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_topology(&C->T);
		C->topo_init=0;
		ssdp_reset_errors();
	}
	return;
}
/*
BEGIN_DESCRIPTION
SECTION Coordinate System
PARSEFLAG determine_utm DetermineUTM "epsg=<out-array> lat=<float-value> lon=<float-value>"
DESCRIPTION Get EPSG code of the UTM system at a given location
ARGUMENT lat latitude coordinates (in WGS84, epsg:4326)
ARGUMENT lon longitude coordinates (in WGS84, epsg:4326)
OUTPUT epsg integer epsg code
END_DESCRIPTION
*/
void DetermineUTM(char *in)
{
        array a;
        double lat, lon;
        int epsg;
        char *word;
        word=malloc((strlen(in)+1)*sizeof(*word));
        if (NULL == word) goto eword;

        if (FetchFloat(in, "lat", word, &lat)) goto epars;
        if (FetchFloat(in, "lon", word, &lon)) goto epars;
        if (!GetArg(in, "epsg", word)) goto epars;

        epsg = determine_utm(lat, lon);
        a.D = malloc(sizeof(*a.D));
        if (NULL == a.D) goto epars;
        a.N = 1;
        a.D[0] = (double) epsg;

        printf("UTM at (%.5f, %.5f) %s := %d\n", lat, lon, word, epsg);
        if (AddArray(word, a)) goto eaddarray;

        return;
eaddarray:
        free(a.D);
epars:
        free(word);
eword:
        return;
}
/*
BEGIN_DESCRIPTION
SECTION Coordinate System
PARSEFLAG convert_epsg ConvertEPSG "x=<in/out-array> y=<in/out-array> epsg_src=<int-value> epsg_dst=<int-value>"
DESCRIPTION Convert coordinate from one projection to another
ARGUMENT x first (latitude) coordinate (in epsg_src)
ARGUMENT y second (longitude) coordinate (in epsg_dst)
ARGUMENT epsg_src source coordinate system (e.g. "4326" corresponds to the WGS84 or GPS system)
ARGUMENT epsg_dst destination coordinate system
OUTPUT x,y output coordinates
END_DESCRIPTION
*/
void ConvertEPSG(char *in)
{
        int i;
        char *word;
        word=malloc((strlen(in)+1)*sizeof(*word));
        if (NULL == word) goto eword;

        int es, ed;
        array *x, *y;
        if (FetchArray(in, "x", word, &x)) goto epars;
        if (FetchArray(in, "y", word, &y)) goto epars;
        if (FetchInt(in, "epsg_src", word, &es)) goto epars;
        if (FetchInt(in, "epsg_dst", word, &ed)) goto epars;

        if (x->N != y->N) {
                Warning("arrays have different lengths!\n"
                        "\tlen(x) = %d\n"
                        "\tlen(y) = %d\n",
                        x->N, y->N);
                goto epars;
        }

        struct epsg *pc = epsg_init_epsg(es, ed);
        if (NULL == pc) {
                Warning("failed to init epsg context"
                       "\t epsg_src = %d"
                       "\t epsg_dst = %d", es, ed);
                goto epars;
        }

        struct point p;
        printf("Converting coordinates from epsg:%d to epsg:%d\n", es, ed);
        for (i=0; i<x->N;++i) {
                p.x = x->D[i];
                p.y = y->D[i];
                convert_point(pc,&p,1);
                x->D[i] = p.x;
                y->D[i] = p.y;
        }

        free(word);
        epsg_free(pc);
        return;
epars:
        free(word);
eword:
        return;
}
/*
BEGIN_DESCRIPTION
SECTION Coordinate System
PARSEFLAG place_template PlaceTemplate "x=<in/out-array> y=<in/out-array> lat=<float-value> lon=<float-value> [azimuth=<float-value>] [epsg=<int-value>]"
DESCRIPTION Rotate template coordinates and place them at a given location. The given coordinates are rotated wrt to the origin in the provided coordinates. The units of the input coordinates are interpreted in the units of the selected epsg coordinate system. The output coordinates are given in WGS84 (epsg:4326) system.
ARGUMENT x,y first and second input coordinates. Output coordiantes are in WGS84, epsg:4326
ARGUMENT lat latitude (in WGS84, epsg:4326)
ARGUMENT lon longitude (in WGS84, epsg:4326)
ARGUMENT azimuth azimuth angle in degrees(default: 0, North)
ARGUMENT epsg optional coordinate system (default UTM according to lat and lon is selected)
OUTPUT x,y output coordinates
END_DESCRIPTION
*/
void PlaceTemplate(char *in)
{
        int epsg;
        double lat, lon, azi;
        array *x, *y;
        char *word;
        word=malloc((strlen(in)+1)*sizeof(*word));
        if (NULL == word) goto eword;

        if (FetchArray(in, "x", word, &x)) goto epars;
        if (FetchArray(in, "y", word, &y)) goto epars;
        if (x->N != y->N) {
                Warning("arrays have different lengths!\n"
                        "\tlen(x) = %d\n"
                        "\tlen(y) = %d\n",
                        x->N, y->N);
                goto epars;
        }

        if (FetchFloat(in, "lat", word, &lat)) goto epars;
        if (FetchFloat(in, "lon", word, &lon)) goto epars;

        if (FetchOptFloat(in, "azimuth", word, &azi))
                azi = 0;
        azi = deg2rad(azi);

        if (FetchOptInt(in, "epsg", word, &epsg))
                epsg = determine_utm(lat, lon);
        struct epsg *pc = epsg_init_epsg(epsg, 4326);
        if (NULL == pc) {
                Warning("failed to init epsg context"
                        "\t epsg = %d", epsg);
                goto epars;
        }

        placetemplate(pc, lat, lon, azi, x->D, y->D, x->N);

        epsg_free(pc);
epars:
        free(word);
eword:
        return;
}


static double* array_copy_at(array *src, int N, double *dst, int at)
{
		int i;
		for (i=0; i < N; ++i)
				dst[at + i] = src->D[i % src->N];

		return dst + at;
}

/*
BEGIN_DESCRIPTION
SECTION Coordinate System
PARSEFLAG place_body PlaceBody "lat=<in-array> lon=<in-array> bearing=<in-array> x/ox=<in/out-array> y/oy=<in/out-array> [z/oz=<in/out-array>] [azimuth/oazimuth=<in/out-array>] [zenith/ozenith=<in/out-array>] [epsg=<int-value>]"
DESCRIPTION Rotate a template of a body coordinates at a given location and direction (bearing). The units of the input coordinates are interpreted in the units of the selected epsg coordinate system. The output coordinates are given in WGS84 (epsg:4326) system. The oazimuth vector is replicated and adjusted by the given bearing. The oz, ozenith vectors are copied without a change. z,oz,azimuth,oazimuth,zenith,ozenith are optional vectors.
ARGUMENT lat,lon latitude and longitude arrays (in WGS84, epsg:4326)
ARGUMENT bearing azimuth at the specific location (in degrees)
ARGUMENT x,y,z,azimuth,zenith input arrays for 3d coordinate and observer orientation (azimuth and zenith are given in degrees)
ARGUMENT ox,oy,oz,oazimuth,ozenith output arrays. len(ox)=len(x)*len(lat). ox[0:len(x)] gives locations for the whole body for the first coordinate. len(ox)=len(oy)=len(oz)=len(oazimuth)=len(ozenith)
ARGUMENT epsg optional coordinate system (default UTM according to lat[0] and lon[0] is selected)
OUTPUT x,y output coordinates
END_DESCRIPTION
*/
void PlaceBody(char *in)
{
		int i, j, epsg;
		double *x, *y;
		char *word, *nox, *noy, *noz, *noazi, *nozen;
		array *lat, *lon, *beta, *ix, *iy, *iz, *iazi, *izen;
		array ox, oy, oz, oazi, ozen;
		ox.D=NULL; oy.D=NULL; oz.D=NULL; oazi.D=NULL; ozen.D=NULL;
		nox=NULL; noy=NULL; noz=NULL; noazi=NULL; nozen=NULL;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (FetchArray(in, "lat", word, &lat)) goto epars;
		if (FetchArray(in, "lon", word, &lon)) goto epars;
		if (FetchArray(in, "bearing", word, &beta))
				goto epars;
		if (lat->N != lon->N || 0==lat->N
			|| (beta->N != lat->N && 1 != beta->N)) {
				Warning("Error: invalid lat,lon,bearing array lengths\n");
				goto epars;
		}

		if (FetchArray(in, "x", word, &ix)) goto epars;
		if (FetchArray(in, "y", word, &iy)) goto epars;
		if (FetchOptArray(in, "z", word, &iz))
				iz = NULL;
		if (FetchOptArray(in, "azimuth", word, &iazi))
				iazi = NULL;
		if (FetchOptArray(in, "zenith", word, &izen))
				izen = NULL;
		if (ix->N!=iy->N || 0==ix->N
			|| (iz && iz->N!=ix->N && 1!=iz->N)
			|| (iazi && iazi->N!=ix->N && 1!=iazi->N)
			|| (izen && izen->N!=ix->N && 1!=izen->N)) {
				Warning("Error: invalid x,y,z,azi,zen array length\n");
				goto epars;
		}

		if (NULL==(nox=malloc((strlen(in)+1)*sizeof(*nox)))) goto enox;
		if (NULL==(noy=malloc((strlen(in)+1)*sizeof(*noy)))) goto enoy;
		if (iz && NULL==(noz=malloc((strlen(in)+1)*sizeof(*noz)))) goto enoz;
		if (iazi && NULL==(noazi=malloc((strlen(in)+1)*sizeof(*noazi))))
				goto enoazi;
		if (izen && NULL==(nozen=malloc((strlen(in)+1)*sizeof(*nozen))))
				goto enozen;

		if (!GetArg(in, "ox", nox)) goto eopars;
		if (!GetArg(in, "oy", noy)) goto eopars;
		if (iz && !GetArg(in, "oz", noz)) goto eopars;
		if (iazi && !GetArg(in, "oazimuth", noazi)) goto eopars;
		if (izen && !GetArg(in, "ozenith", nozen)) goto eopars;

		if (FetchOptInt(in, "epsg", word, &epsg))
				epsg = determine_utm(lat->D[0], lon->D[0]);

		struct epsg *pc = epsg_init_epsg(epsg, 4326);
		if (NULL == pc) {
				Warning("Error: failed to init epsg context"
						"\t epsg = %d", epsg);
				goto eepsg;
		}

		if (NULL==(ox.D=malloc(lat->N*ix->N*sizeof(*ox.D)))) goto eox;
		if (NULL==(oy.D=malloc(lat->N*ix->N*sizeof(*oy.D)))) goto eoy;
		if (iz && NULL==(oz.D=malloc(lat->N*ix->N*sizeof(*oz.D)))) goto eoz;
		if (iazi && NULL==(oazi.D=malloc(lat->N*ix->N*sizeof(*oazi.D))))
				goto eoazi;
		if (izen && NULL==(ozen.D=malloc(lat->N*ix->N*sizeof(*ozen.D))))
				goto eozen;
		ox.N = oy.N = oz.N = oazi.N = ozen.N = lat->N*ix->N;

		for (i=0; i<lat->N; ++i) {
				x = array_copy_at(ix, ix->N, ox.D, i * ix->N);
				y = array_copy_at(iy, ix->N, oy.D, i * ix->N);
				if (iz)
						array_copy_at(iz, ix->N, oz.D, i * ix->N);
				if (iazi)
						array_copy_at(iazi, ix->N, oazi.D, i * ix->N);
				if (izen)
						array_copy_at(izen, ix->N, ozen.D, i * ix->N);

				for (j=0; iazi && j < ix->N; ++j)
						oazi.D[i*ix->N + j] += beta->D[i%beta->N];

				placetemplate(pc, lat->D[i], lon->D[i],
							  deg2rad(beta->D[i%beta->N]), x, y, ix->N);
		}

		if (izen && AddArray(nozen, ozen)) goto eozen;
		if (iazi && AddArray(noazi, oazi)) goto enozen;
		if (iz && AddArray(noz, oz)) goto enoazi;
		if (AddArray(noy, oy)) goto enoz;
		if (AddArray(nox, ox)) goto enoy;

		epsg_free(pc);
		free(word);
		return;
		free(ozen.D); ozen.D=NULL;
eozen:
		free(oazi.D); oazi.D=NULL;
eoazi:
		free(oz.D); oz.D=NULL;
eoz:
		free(oy.D); oy.D=NULL;
eoy:
		free(ox.D); ox.D=NULL;
eox:
		epsg_free(pc);
eepsg:
eopars:
		free(nozen);
		free(ozen.D);
enozen:
		free(noazi);
		free(oazi.D);
enoazi:
		free(noz);
		free(oz.D);
enoz:
		free(noy);
		free(oy.D);
enoy:
		free(nox);
		free(ox.D);
enox:
		free(word);
epars:
eword:
		Warning("Error: place_body failed!\n");
		return;
}
