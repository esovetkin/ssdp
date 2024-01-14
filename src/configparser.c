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


static void set_chunk(simulation_config *C, int chunkid)
{
		if (C->chunked < 0) {
				C->Nl_eff = C->Nl;
				C->Nl_o = 0;
				return;
		}

		if (chunkid < 0) {
				C->Nl_eff = 0;
				C->Nl_o = 0;
				return;
		}

		C->Nl_eff = C->chunked;
		C->Nl_o = chunkid * C->chunked;

		if (C->Nl_o + C->chunked > C->Nl)
				C->Nl_eff = C->Nl - C->Nl_o;
}


int NextChunk(simulation_config *C, int *chunkid)
{
		++(*chunkid);
		printf("Processing locations chunk: %d\n", *chunkid);

		if (C->chunked < 0) {
				if (0 == (*chunkid)) return 1;
				return 0;
		}

		if (C->chunked * (*chunkid) >= C->Nl)
				return 0;

		return 1;
}


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
		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (!GetArg(in, "C", word)) goto eargs;

		C=InitConf();
		printf("Defining simulation configuration \"%s\"\n", word);
		// note InitConf does not allocate anything, no need to free
		if (AddSimConf(word, C)) {
				free(word);
				goto eargs;
		}

		return;
eargs:
		free(word);
eword:
		Warning("Error: init_sim_config failed!\n");
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
		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (FetchFloat(in, "lon", word, &l)) goto eargs;
		C->lon=l;
		printf("set longitude to %e degrees (%e rad)\n", rad2deg(C->lon), C->lon);

		if (FetchFloat(in, "lat", word, &l)) goto eargs;
		C->lat=l;
		printf("set latitude to %e degrees (%e rad)\n", rad2deg(C->lat), C->lat);

		if (FetchOptFloat(in, "E", word, &l)) l = 0;
		C->E=l;
		printf("set elevation to %e m\n", C->E);

		free(word);
		return;
eargs:
		free(word);
eword:
		Warning("Error: config_coord failed!\n");
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
		if (C->L) {
				for (i=0; i < C->Nl_eff; ++i)
						ssdp_free_location(&(C->L[i]));
				free(C->L);
				C->L=NULL;
		}
		if (C->l_hcache) {
				hashmap_free(C->l_hcache, free);
				C->l_hcache=NULL;
		}
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


static int init_stcache(simulation_config *C)
{
		int i, hits=0, uh=0;

		if (C->l_stcache) {hashmap_free(C->l_stcache, free); C->l_stcache=NULL;}
		if (C->l_st) {free(C->l_st); C->l_st=NULL;}
		if (NULL==(C->l_stcache=hashmap_init(2*C->Nl_eff))) goto el_stcache;
		if (NULL==(C->l_st=calloc(C->Nl_eff, sizeof(*C->l_st)))) goto el_st;

		TIC();
		// assume C->stcache initialised. associate location with the
		sky_transfer *st;
		uintptr_t k;
		struct l_stcache_data *lst;

		// allocated space for that in the hcache.
		for (i=0; i < C->Nl_eff; ++i) {
				st = ssdp_stcache_get
						(C->stcache,
						 fmod(C->o[i+C->Nl_o].z+2*M_PI, 2*M_PI),
						 fmod(C->o[i+C->Nl_o].a+2*M_PI, 2*M_PI), C->albedo);
				if (NULL == st) goto estcache;
				if (NULL != st->t) ++hits;
				C->l_st[i] = st;

				k = (uintptr_t) st;
				if (hashmap_isin(C->l_stcache, &k, sizeof(k))) {
						++uh;
						continue;
				}

				if (NULL == (lst=malloc(sizeof(*lst)))) goto estcache;
				*lst = (struct l_stcache_data){.ST=st,.Li=i};

				if (hashmap_insert(C->l_stcache, &k, sizeof(k), lst)) goto estcache;
		}
		printf("Checked in sky transfer cache in %g s, hits/duplicates/missed: %d/%d/%d\n", TOC(), hits, uh, C->Nl-hits);

		return 0;
estcache:
		// the stcache keeps track of allocated horizon. even if we
		// fail and end up here, just need to cleanup the
		// C->l_st. ssdp_stcache_free at some point will cleanup
		// everything.
		free(C->l_st); C->l_st = NULL;
el_st:
		hashmap_free(C->l_stcache, free); C->l_stcache=NULL;
el_stcache:
		return -1;
}


static int init_hcache(simulation_config *C)
{
		int i, hits=0, uh = 0;

		if (NULL==(C->L=calloc(C->Nl_eff,sizeof(*(C->L))))) goto ecl;
		if (NULL==(C->l_hcache=hashmap_init(2*C->Nl_eff))) goto el_hcache;

		horizon *H;
		uintptr_t k;
		struct l_hcache_data *lh;

		TIC();
		// assume C->hcache initialised. associate location with the
		// allocated space for that in the hcache.
		for (i=0; i < C->Nl_eff; ++i) {
				H = ssdp_horizoncache_get(C->hcache,
										  C->x[i+C->Nl_o],
										  C->y[i+C->Nl_o],
										  C->z[i+C->Nl_o]);
				if (NULL == H) goto ehcache;
				if (NULL != H->zen) ++hits;
				C->L[i].H = H;

				k = (uintptr_t) H;
				if (hashmap_isin(C->l_hcache, &k, sizeof(k))) {
						++uh;
						continue;
				}

				if (NULL == (lh=malloc(sizeof(*lh)))) goto ehcache;
				*lh = (struct l_hcache_data){.H=H,.Li=i};

				if (hashmap_insert(C->l_hcache, &k, sizeof(k), lh)) goto ehcache;
		}
		printf("Checked in locations cache in %g s, hits/duplicates/misses: %d/%d/%d\n",
			   TOC(), hits, uh, C->Nl_eff-hits);

		return 0;
		// the hcache keeps track of allocated horizon. even if we
		// fail and end up here, just need to cleanup the
		// C->L. ssdp_horizoncache_free at some point will cleanup
		// everything.
ehcache:
		hashmap_free(C->l_hcache, free);
		C->l_hcache = NULL;
el_hcache:
		free(C->L); C->L = NULL;
ecl:
		return -1;
}


static double init_initst(simulation_config *C)
{
		int i, pco;
		if (NULL==C->stcache)
				// xydelta - sky angles, zdelta - albedo
				if (NULL==(C->stcache=ssdp_rtreecache_init(0.001, 0.001)))
						goto estcache;
		if (init_stcache(C)) goto estcache;

		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->l_stcache->cap; ++i) {
						struct l_stcache_data *lst = (struct l_stcache_data*) C->l_stcache->values[i];
						if (NULL == lst) continue;
						if (NULL != lst->ST->t)  continue;

						ssdp_init_transfer
								(lst->ST, C->S, C->albedo,
								 C->o[lst->Li+C->Nl_o], &(C->M));
						ProgressBar((100*(i+1))/C->l_stcache->cap, &pco, ProgressLen, ProgressTics);
				}
		}
		ProgressBar(50, &pco, ProgressLen, ProgressTics);
		if (ssdp_error_state) goto error;

		return TOC();
error:
estcache:
		return -1.0;
}


static int init_transfer(simulation_config *C)
{
		double dt;
		int i, pco=0;

		if (0 > (dt=init_initst(C))) goto einitst;

		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->Nl_eff; ++i) {
						ssdp_setup_transfer
								(&(C->L[i]), C->S, C->l_st[i]);
						ProgressBar((100*(i+1+C->Nl_eff))/(2*C->Nl_eff), &pco, ProgressLen, ProgressTics);
				}
		}
		ProgressBar(100, &pco, ProgressLen, ProgressTics);
		if (ssdp_error_state) goto error;
		dt+=TOC();
		printf("Initialised sky transfer in %g s (%e s/locations)\n", dt, dt/((double)C->Nl_eff));

		return 0;
error:
einitst:
		return -1;
}


int InitConfigMask(simulation_config *C, int chunkid)
{
		double dt;
		int i, pco = 0;
		set_chunk(C, chunkid);

		if (!C->Nl_eff) return 0;
		if (!((C->sky_init)&&(C->topo_init)&&(C->loc_init))) goto egtl;
		if (init_hcache(C)) goto einitL;

		printf("Tracing %d locations\n", C->Nl_eff);
		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->l_hcache->cap; ++i) {
						struct l_hcache_data *lh = (struct l_hcache_data*) C->l_hcache->values[i];
						if (NULL == lh) continue;
						if (NULL != lh->H->zen) continue;

						ssdp_setup_horizon(lh->H, C->S, &(C->T),
										   C->x[lh->Li + C->Nl_o],
										   C->y[lh->Li + C->Nl_o],
										   C->z[lh->Li + C->Nl_o]);
						ProgressBar((100*(i+1))/C->l_hcache->cap, &pco, ProgressLen, ProgressTics);
				}
		}
		ProgressBar(100, &pco, ProgressLen, ProgressTics);
		dt = TOC();
		printf("Traced %d horizons in %g s (%e s/horizons)\n", C->Nl_eff, dt, dt/((double)C->Nl_eff));

		if (init_transfer(C)) goto estate;
		return 0;
estate:
		FreeConfigMask(C); // make sure we are clear to allocate new memory
einitL:
egtl:
		return -1;
}


int InitConfigGridMask(simulation_config *C, int chunkid)
{
		double dt;
		int i, pco=0;
		set_chunk(C, chunkid);

		if (!C->Nl_eff) return 0;
		if (!((C->sky_init)&&(C->grid_init)&&(C->loc_init))) goto esgl;

		TIC();
		i = ssdp_topogrid_approxhorizon
				(C->Tx, C->nTx, C->approx_n, C->approx_stype, C->S);
		dt=TOC();

		switch (i) {
		case -1:
				goto enapprox;
				break;
		case 0:
				break;
		case -2:
				ssdp_rtreecache_reset(&(C->hcache));
				break;
		default:
				printf("Initialised topography sample set in %g s\n", dt);
				ssdp_rtreecache_reset(&(C->hcache));
				break;
		}


		if (init_hcache(C)) goto einitL;

		// process all necessary horizons separately
		printf("Tracing %d locations\n", C->Nl_eff);
		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->l_hcache->cap; ++i) {
						struct l_hcache_data *lh = (struct l_hcache_data *) C->l_hcache->values[i];
						if (NULL == lh) continue;
						if (NULL != lh->H->zen) continue;

						ssdp_setup_grid_horizon(
								lh->H, C->S, C->Tx, C->nTx,
								C->x[lh->Li + C->Nl_o],
								C->y[lh->Li + C->Nl_o],
								C->z[lh->Li + C->Nl_o]);
						ProgressBar((100*(i+1))/C->l_hcache->cap, &pco, ProgressLen, ProgressTics);
				}
		}
		ProgressBar(100, &pco, ProgressLen, ProgressTics);
		dt = TOC();
		printf("Traced %d horizons in %g s (%e s/horizons)\n", C->Nl_eff, dt, dt/((double)C->Nl_eff));

		if (init_transfer(C)) goto estate;
		return 0;
estate:
		FreeConfigMask(C); // make sure we are clear to allocate new memory
einitL:
enapprox:
esgl:
		return -1;
}


// same as above but without horizon calculation
int InitConfigMaskNoH(simulation_config *C, int chunkid)
{
		double dt;
		int i, pco=0;
		set_chunk(C, chunkid);

		if (!C->Nl_eff) return 0;
		if (!((C->sky_init)&&(C->loc_init))) goto esl;
		if (init_hcache(C)) goto einitL;

		// some horizon needs to be initialised
		printf("Initialising %d horizons\n", C->Nl_eff);
		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->l_hcache->cap; ++i) {
						struct l_hcache_data *lh = (struct l_hcache_data *) C->l_hcache->values[i];
						if (NULL == lh) continue;
						if (NULL != lh->H->zen) continue;

						ssdp_setup_horizon(lh->H, C->S, NULL,
										   C->x[lh->Li + C->Nl_o],
										   C->y[lh->Li + C->Nl_o],
										   C->z[lh->Li + C->Nl_o]);
						ProgressBar((100*(i+1))/C->l_hcache->cap, &pco, ProgressLen, ProgressTics);
				}
		}
		ProgressBar(100, &pco, ProgressLen, ProgressTics);
		dt = TOC();
		printf("Initialised %d horizons in %g s (%e s/horizons)\n", C->Nl_eff, dt, dt/((double)C->Nl_eff));

		if (init_transfer(C)) goto estate;
		return 0;
estate:
		FreeConfigMask(C); // make sure we are clear to allocate new memory
einitL:
esl:
		return -1;
}


/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_sky ConfigSKY "C=<out-config> N=<int-value>"
DESCRIPTION Setup the sky. This commands allocates space and initializes the sky data.
ARGUMENT N The number of zenith discretizations. The total number of sky patches equals Ntotal(N)=3N(N+1)+1, e.g. with Ntotal(7)=127. The first layer of sky induces the azimuth discretisation, which equals 6N.
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigSKY(char *in)
{
		simulation_config *C;
		char *word;
		int N;
		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (FetchInt(in, "N", word, &N)) goto eargs;
		if (N<=1) {
				Warning("Error: number of zenith discretizations must me larger than 1\n");
				goto eargs;
		}

		if (C->sky_init)
				ssdp_free_sky(C->S, C->nS);
		else
				C->sky_init=1;

		// locations and horizon are no longer valid with the new sky
		ssdp_rtreecache_reset(&(C->hcache));
		ssdp_rtreecache_reset(&(C->stcache));
		FreeConfigLocation(C);
		printf("Cleared configured locations, new sky renders them invalid\n");

#ifdef OPENMP
		C->nS = omp_get_max_threads();
#else
		C->nS = 1;
#endif

		C->S=ssdp_init_sky(N, C->nS);
		printf("Configuring sky with %d zenith discretizations\n", N);

		if (ssdp_error_state) {
				ssdp_print_error_messages();
				ssdp_free_sky(C->S, C->nS);
				C->sky_init=0;
				C->nS=0;
				ssdp_reset_errors();
		}

		free(word);
		return;
eargs:
		free(word);
eword:
		Warning("Error: config_sky failed!\n");
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

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (FetchArray(in, "x", word, &x)) goto eargs;
		if (FetchArray(in, "y", word, &y)) goto eargs;
		if (FetchArray(in, "z", word, &z)) goto eargs;

		if ((x->N!=y->N)||(x->N!=z->N)) {
				Warning("x- y- and z-arrays must be of equal length\n");
				goto eargs;
		}

		if (C->topo_init) ssdp_free_topology(&C->T);
		C->topo_init=1;

		printf("Configuring topology with %d points\n", x->N);
		C->T=ssdp_make_topology(x->D, y->D, z->D, x->N);
		if (ssdp_error_state) goto emktopo;

		ssdp_rtreecache_reset(&(C->hcache));
		if (InitLocations(C, -1)) goto einitloc;

		free(word);
		return;
einitloc:
emktopo:
		ssdp_free_topology(&C->T);
		C->topo_init=0;
		ssdp_print_error_messages();
		ssdp_reset_errors();
eargs:
		free(word);
eword:
		printf("Error: config_topology failed!\n");
}


static int resetadd_topogrid(simulation_config *C, int ifconfig)
{
		if (ifconfig) {
				if (C->grid_init) {
						ssdp_free_topogrid(C->Tx, C->nTx);
						free(C->Tx); C->Tx=NULL; C->nTx=0;
				} else {
						C->grid_init=1;
				}
		}

		if (0 == C->nTx) {
				if (NULL==(C->Tx=malloc(sizeof(*C->Tx)))) goto emalloc;
				C->nTx = 1;
				return 0;
		}

		topogrid *tmp = realloc(C->Tx, (C->nTx+1)*sizeof(*tmp));
		if (NULL == tmp) goto emalloc;
		C->Tx = tmp;
		++C->nTx;

		return 0;
emalloc:
		return -1;
}


static void configadd_topogrid(char *in, int ifconfig)
{
		simulation_config *C;
		array *z;
		double x1, y1, x2, y2;
		int Nx, Ny;
		char *word;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (FetchFloat(in, "x1", word, &x1)) goto eargs;
		if (FetchFloat(in, "y1", word, &y1)) goto eargs;
		if (FetchFloat(in, "x2", word, &x2)) goto eargs;
		if (FetchFloat(in, "y2", word, &y2)) goto eargs;
		if (FetchInt(in, "Nx", word, &Nx)) goto eargs;
		if (FetchInt(in, "Ny", word, &Ny)) goto eargs;
		if (FetchArray(in, "z", word, &z)) goto eargs;

		if (z->N!=Nx*Ny)
		{
				Warning("Number of steps in x- and y- directions "
						"do not match the number of elements in "
						"the z array (%d != %d x %d)\n",
						z->N, Nx, Ny);
				goto eargs;
		}

		if (resetadd_topogrid(C, ifconfig)) goto eTx;

		printf("Configuring topogrid with %d points\n", z->N);
		C->Tx[C->nTx-1]=ssdp_make_topogrid(z->D, x1, y1, x2, y2, Nx, Ny);
		ssdp_min_topogrids(C->Tx, C->nTx);
		if (ssdp_error_state) goto essdp;
		if (ssdp_rtreecache_reset(&(C->hcache))) goto essdp;
		if (InitLocations(C, -1)) goto essdp;

		free(word);
		return;
essdp:
		ssdp_print_error_messages();
		ssdp_free_topogrid(C->Tx, C->nTx);
		free(C->Tx); C->Tx=NULL; C->nTx=0;
		C->grid_init=0;
		ssdp_reset_errors();
eTx:
eargs:
		free(word);
eword:
		if (ifconfig)
				Warning("Error: config_topogrid failed!\n");
		else
				Warning("Error: add_topogrid failed!\n");
		return;
}


static void configadd_topogdal(char *in, int ifconfig)
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
		if (step < 0) {
				Warning("Error: step must be positive!\n");
				goto epars;
		}

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

		if (resetadd_topogrid(C, ifconfig)) goto eTx;

		TIC();
		C->Tx[C->nTx-1] = ssdp_make_topogdal
				(x1, y1, x2, y2, fns->s, fns->n, step, epsg);
		ssdp_min_topogrids(C->Tx, C->nTx);
		if (ssdp_error_state) goto emaketopogdal;
		if (ssdp_rtreecache_reset(&(C->hcache))) goto emaketopogdal;
		printf("Initialised topogrid in %g s\n", TOC());

		if (InitLocations(C, -1)) goto emaketopogdal;

		cvec_free(fns);
		free(word);
		return;
emaketopogdal:
		ssdp_print_error_messages();
		ssdp_free_topogrid(C->Tx, C->nTx);
		free(C->Tx); C->Tx=NULL; C->nTx=0;
		C->grid_init=0;
		ssdp_reset_errors();
eTx:
efnlist:
efnspush:
		cvec_free(fns);
ecvec:
epars:
		free(word);
eword:
		if (ifconfig)
				Warning("Error: config_topogdal failed!\n");
		else
				Warning("Error: add_topogdal failed!\n");
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
		configadd_topogrid(in, 1);
}

/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG add_topogrid AddTOPOGrid "C=<out-config> z=<in-array> Nx=<int-value> Ny=<int-value> x1=<float-value> y1=<float-value> x2=<float-value> y2=<float-value>"
DESCRIPTION Add topogrid the existing topogrid. Same arguments as in config_topogrid
ARGUMENT all same as in config_topogrid
END_DESCRIPTION
*/
void AddTOPOGrid (char *in)
{
		configadd_topogrid(in, 0);
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
		configadd_topogdal(in, 1);
}

/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG add_topogdal AddTOPOGDAL "C=<out-config> lat1=<float-value> lon1=<float-value> lat2=<float-value> lon2=<float-value> step=<float-value> f0=<file-str> f1=<file-str> .. fN=<file-str> [flist=<file-str>] [epsg=<int-value>]"
DESCRIPTION Add topography with GDAL to existing topogrid. Arguments same as in config_topogdal
END_DESCRIPTION
*/
void AddTOPOGDAL (char *in)
{
		configadd_topogdal(in, 0);
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


int InitLocations(simulation_config *C, int chunkid)
{
		// case when called for sim_ routines
		if ((-1 == C->chunked) && (chunkid > -1)) return 0;

		// case when called from config_location, postpone tracing for
		// later
		if ((C->chunked > 0) && (-1 == chunkid)) return 0;

		FreeConfigMask(C); // make sure we are clear to allocate new memory
		int err = 0;
		if ((!C->topo_init) && (!C->grid_init))
		{
				Warning("Warning: no topological data "
						"available, omitting horizon\n");
				err |= InitConfigMaskNoH(C, chunkid);
		}

		if ((!C->grid_init) && (C->topo_init))
				err |= InitConfigMask(C, chunkid);
		if (C->grid_init)
				err |= InitConfigGridMask(C, chunkid);

		if (ssdp_error_state || err) goto elocs;
		return 0;
elocs:
		ssdp_print_error_messages();
		FreeConfigLocation(C);
		C->loc_init=0;
		ssdp_reset_errors();

		return -1;
}


/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_locations ConfigLoc "C=<out-config> x=<in-array> y=<in-array> z=<in-array> azimuth=<in-array> zenith=<in-array> [albedo=<float-value>] [xydelta=<float-value>] [zdelta=<float-value] [approx_n=<int-value>] [approx_type=<str>] [chunked=<int-value>]"
DESCRIPTION Setup the topography. Load the x, y, and z data of the unstructured topography mesh into the configuration data.
ARGUMENT x x coordinates
ARGUMENT y y coordinates
ARGUMENT z z coordinates
ARGUMENT azimuth azimuth angle of tilted surface
ARGUMENT zenith zenith angle of tilted surface
ARGUMENT albedo optionally provide an albedo value between 0-1
ARGUMENT xydelta,zdelta the coordinates within xydelta in xy plane and zdelta within z direction are considered the same (default: 0.05)
ARGUMENT approx_n optional, if positive determine number of raster points used for computing the horizon. For sample points are used polar Sobol 2-d set (s_1, s_2), where pixel location is computed using F^{-1}(s_1)*exp(1i*2*pi*s_2), where F^{-1} is provided inverse cumulative distribution function, see `horizon_sample_dstr` (default: 10000)
ARGUMENT approx_type type of sampling used. Either "precise", "sobol", "iid", "rays16", "rays32", "rays64", and "rays128" (default: "precise")
ARGUMENT chunked optional, if positive the locations are processed in chunks of given size. This allows to do simulations with more locations with less memory (default: -1)
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigLoc(char *in)
{
		simulation_config *C;
		double xydelta, zdelta;
		array *x, *y, *z, *az, *ze;
		int N, i, approx_n, chunked=-1;
		char *word;
		enum SampleType stype = PRECISE;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (FetchArray(in, "x", word, &x)) goto eargs;
		if (FetchArray(in, "y", word, &y)) goto eargs;
		if (FetchArray(in, "z", word, &z)) goto eargs;
		if (FetchArray(in, "azimuth", word, &az)) goto eargs;
		if (FetchArray(in, "zenith", word, &ze)) goto eargs;
		if (FetchOptFloat(in, "albedo", word, &(C->albedo))) C->albedo=0.0;
		if (C->albedo < 0 || C->albedo > 1)
				Warning("Warning: albedo=%f lies outside [0,1]\n", C->albedo);
		if (FetchOptFloat(in, "xydelta", word, &(xydelta))) xydelta=0.05;
		if (FetchOptFloat(in, "zdelta", word, &(zdelta))) zdelta=0.05;
		if (FetchOptInt(in, "approx_n", word, &approx_n)) approx_n = 10000;
		if (GetOption(in, "approx_type", word)) {
				if (strncmp(word,"precise",10)==0) stype=PRECISE;
				else if (strncmp(word,"sobol",10)==0) stype=SOBOL;
				else if (strncmp(word,"iid",10)==0) stype=IID;
				else if (strncmp(word,"rays16",10)==0) stype=RAYS16;
				else if (strncmp(word,"rays32",10)==0) stype=RAYS32;
				else if (strncmp(word,"rays64",10)==0) stype=RAYS64;
				else if (strncmp(word,"rays128",10)==0) stype=RAYS128;
				else Warning("Warning: unsupported approx_type=%s, "
							 "continue with precise\n", word);
		}
		if (FetchOptInt(in, "chunked", word, &chunked)) chunked = -1;

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
		TIC();
		C->xyorder = ssdp_init_xyorder(C->Tx, C->x, C->y, C->Nl);
		printf("Sorted xy locations according in %g s\n", TOC());
		if (NULL==C->xyorder) goto exyorder;

		if (NULL==C->hcache)
				if (NULL==(C->hcache=ssdp_rtreecache_init(xydelta, zdelta)))
						goto ehcache;
		C->hcache->xydelta = xydelta / 2.0;
		C->hcache->zdelta = zdelta / 2.0;
		C->approx_n = approx_n;
		C->approx_stype = stype;
		C->chunked = chunked;

		if (InitLocations(C, -1) < 0) goto elocs;

		free(word);
		return;
elocs:
ehcache:
		free(C->xyorder); C->xyorder = NULL;
exyorder:
		free(C->o); C->o = NULL;
eo:
		free(C->z); C->z = NULL;
ez:
		free(C->y); C->y = NULL;
ey:
		free(C->x); C->x = NULL;
ex:
eargs:
		free(word);
eword:
		Warning("config_locations: failed!\n");
		return;
}

/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG horizon_sample_distr HorizonSampleDistr "C=<in-config> [q=<float-array>] [phi=<float-array>] [rasterid=<int-value>]"
DESCRIPTION Set the radial distribution for the approximate horizon sampling. The distribution is set using an array q of quantiles. The distribution can be set per raster. Can only be used with topogrid topogrpaphies.
ARGUMENT C simulation configuration with topogrid
ARGUMENT q an array of quantiles for radial distribution. Must be a monotonically increasing vector, where q[0]: 0-quantile, q[len(q)-1]: 1-quantule.
ARGUMENT phi an array of quantiles for azimuthal distribution. Must be a monotonically increasing vector and must generate distribution in the interval [0,1]
ARGUMENT rasterid optional id of the raster, where the distribution is set (default: 0)
OUTPUT C distribution law of horizon sample set
END_DESCRIPTION
*/
void HorizonSampleDistr(char *in)
{
		int i, r=0;
		char *word;
		simulation_config *C;
		array *q, *phi;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
        if (FetchOptArray(in, "q", word, &q)) q = NULL;
        if (FetchOptArray(in, "phi", word, &phi)) phi = NULL;
		if (FetchOptInt(in, "rasterid", word, &r)) r=0;

		if (C->grid_init==0) {
				Warning("ERROR: simulation config does not contain a topogrid\n");
				goto eargs;
		}
		if (r >= C->nTx) {
				Warning("ERROR: only %d topogrids configured\n", C->nTx);
				goto eargs;

		}

		if (!q && !phi) {
				Warning("Error: neither q nor phi quantiles are provded!\n");
				goto eargs;
		}

		if (q)
				for (i=1; i < q->N; ++i)
						if (q->D[i-1] > q->D[i]) {
								Warning("Error: provided quantiles in q are not monotone\n");
								goto eargs;
						}
		if (phi)
				for (i=1; i < phi->N; ++i)
						if (phi->D[i-1] > phi->D[i]) {
								Warning("Error: provided quantiles in phi are not monotone\n");
								goto eargs;
						}

		TIC();
		if (q && !phi && ssdp_topogrid_approxlaw(C->Tx+r, q->D, q->N, NULL, 0)) goto esetlaw;
		if (!q && phi && ssdp_topogrid_approxlaw(C->Tx+r, NULL, 0, phi->D, phi->N)) goto esetlaw;
		if (q && phi && ssdp_topogrid_approxlaw(C->Tx+r, q->D, q->N, phi->D, phi->N)) goto esetlaw;
		printf("Sampled new horizon sample in %g s\n", TOC());
		ssdp_rtreecache_reset(&(C->hcache));
		printf("Reset horizon cache, since there new sampling distribution\n");

		free(word);
		return;
esetlaw:
eargs:
		free(word);
eword:
		Warning("Error: horizon_sample_distr failed!\n");
}


/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG export_horizon_sampleset ExportHorizSample "C=<in-config> z=<out-array> nx=<out-1dim-array> ny=<out-1dim-array> [rasterid=<int-value>]"
DESCRIPTION Export topography sample set, which is used for horizon computations. Only works when config_locations has been called.
OUTPUT z,nx,ny binary mask array and its dimensions
END_DESCRIPTION
*/
void ExportHorizSample(char *in)
{
		int i, r=0;
		char *word, *oz, *onx, *ony;
		simulation_config *C;
		array z, nx, ny;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(oz=malloc((strlen(in)+1)*sizeof(*oz)))) goto eoz;
		if (NULL==(onx=malloc((strlen(in)+1)*sizeof(*onx)))) goto eonx;
		if (NULL==(ony=malloc((strlen(in)+1)*sizeof(*ony)))) goto eony;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (!GetArg(in, "z", oz)) goto eargs;
		if (!GetArg(in, "nx", onx)) goto eargs;
		if (!GetArg(in, "ny", ony)) goto eargs;
		if (FetchOptInt(in, "rasterid", word, &r)) r=0;

		if (C->grid_init==0) {
				Warning("ERROR: simulation config does not contain a topogrid\n");
				goto eargs;
		}
		if (r >= C->nTx) {
				Warning("ERROR: only %d topogrids configured\n", C->nTx);
				goto eargs;

		}
		if (NULL == C->Tx[r].horizon_sample) {
				Warning("ERROR: approximate horizon sampling set is not computed yet\n");
				goto eargs;
		}

		int k, Nx = C->Tx[r].Nx, Ny = C->Tx[r].Ny;
		z.N = (2*Nx-1)*(2*Ny-1);
		if (NULL==(z.D=calloc(z.N, sizeof(*(z.D))))) goto ez;
		nx.N = 1;
		if (NULL==(nx.D=malloc(nx.N*sizeof(*(nx.D))))) goto enx;
		nx.D[0] = 2*C->Tx[r].Nx-1;
		ny.N = 1;
		if (NULL==(ny.D=malloc(ny.N*sizeof(*(ny.D))))) goto eny;
		ny.D[0] = 2*C->Tx[r].Ny-1;

		for (i=0; i < C->Tx[r].horizon_nsample_eff; ++i) {
				k = (Nx - 1 + C->Tx[r].horizon_sample[i].x)*(2*Ny-1) + Ny - 1 + C->Tx[r].horizon_sample[i].y;
				z.D[k] = (double) 1.0;
		}

		if(AddArray(ony, ny)) goto eo;
		if(AddArray(onx, nx)) goto eo;
		if(AddArray(oz, z)) goto eo;

		free(word);
		return;
eo:
		free(ny.D);
eny:
		free(nx.D);
enx:
		free(z.D);
ez:
eargs:
		free(ony);
eony:
		free(onx);
eonx:
		free(oz);
eoz:
		free(word);
eword:
		Warning("export_horizon_sampleset: failed!\n");
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

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (FetchFloat(in, "dx", word, &dx)) goto eargs;
		if (FetchFloat(in, "dy", word, &dy)) goto eargs;
		if (FetchFloat(in, "dz", word, &dz)) goto eargs;
		if (FetchInt(in, "N1", word, &N1)) goto eargs;
		if (FetchInt(in, "N2", word, &N2)) goto eargs;

		if ((N1<3)||(N2<3)) {
				Warning("N1 and N2 must be larger than 3 in rand_topology\n");
				goto eargs;
		}
		if ((dx<1e-10)||(dy<1e-10)||(dz<1e-10)) {
				Warning("Please provide valid positive ranges for the random topology\n");
				goto eargs;
		}

		if (C->topo_init) ssdp_free_topology(&C->T);
		C->topo_init=1;

		printf("Configuring topology with %d points\n", N2);
		C->T=ssdp_make_rand_topology(dx, dy, dz, N1, N2);
		if (ssdp_error_state) goto emktopo;
		if (InitLocations(C, -1) < 0) goto einitloc;

		free(word);
		return;
einitloc:
emktopo:
		ssdp_free_topology(&C->T);
		C->topo_init=0;
		if (ssdp_error_state) {
				ssdp_print_error_messages();
				ssdp_reset_errors();
		}
eargs:
		free(word);
eword:
		printf("Error: rand_topology failed\n");
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

		int npc = 1;
#ifdef OPENMP
		npc = omp_get_max_threads();
#endif

		struct epsg **pc = calloc(npc, sizeof(*pc));
		if (NULL==pc) goto epc;
		for (i=0; i < npc; ++i) {
				pc[i] = epsg_init_epsg(es, ed);
		}
		if (npc != i) {
				Warning("failed to init epsg context"
						"\t epsg_src = %d"
						"\t epsg_dst = %d", es, ed);
				goto eepsginit;
		}

		printf("Converting coordinates from epsg:%d to epsg:%d\n", es, ed);
		TIC();
#pragma omp parallel private(i)
		{
				int thread = 0;
				struct point p;

#pragma omp for schedule(runtime)
				for (i=0; i<x->N;++i) {
#ifdef OPENMP
						thread=omp_get_thread_num();
#endif

						p.x = x->D[i];
						p.y = y->D[i];
						convert_point(pc[thread],&p,1);
						x->D[i] = p.x;
						y->D[i] = p.y;
				}
		}
		printf("Converted %d coordinates in %g s\n", x->N, TOC());

		free(word);
		for (i=0; i < npc; ++i)
				epsg_free(pc[i]);
		free(pc);
		return;
eepsginit:
		for (i=0; i < npc; ++i)
				if (pc[i]) epsg_free(pc[i]);
		free(pc);
epc:
epars:
		free(word);
eword:
		printf("Error: convert_epsg failed!\n");
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
		if (iz && (NULL==(noz=malloc((strlen(in)+1)*sizeof(*noz))))) goto enoz;
		if (iazi && (NULL==(noazi=malloc((strlen(in)+1)*sizeof(*noazi)))))
				goto enoazi;
		if (izen && (NULL==(nozen=malloc((strlen(in)+1)*sizeof(*nozen)))))
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
		if (iz && (NULL==(oz.D=malloc(lat->N*ix->N*sizeof(*oz.D))))) goto eoz;
		if (iazi && (NULL==(oazi.D=malloc(lat->N*ix->N*sizeof(*oazi.D)))))
				goto eoazi;
		if (izen && (NULL==(ozen.D=malloc(lat->N*ix->N*sizeof(*ozen.D)))))
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

		if (izen && AddArray(nozen, ozen)) {free(nozen); nozen=NULL; free(ozen.D); ozen.D=NULL;}
		if (iazi && AddArray(noazi, oazi)) {free(noazi); noazi=NULL; free(oazi.D); oazi.D=NULL;}
		if (iz && AddArray(noz, oz)) {free(noz); noz=NULL; free(oz.D); oz.D=NULL;}
		if (AddArray(noy, oy)) {free(noy); noy=NULL; free(oy.D); oy.D=NULL;}
		if (AddArray(nox, ox)) {free(nox); nox=NULL; free(ox.D); ox.D=NULL;}

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
enozen:
		free(noazi);
enoazi:
		free(noz);
enoz:
		free(noy);
enoy:
		free(nox);
enox:
		free(word);
epars:
eword:
		Warning("Error: place_body failed!\n");
		return;
}
