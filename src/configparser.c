#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
DESCRIPTION Setup the Angle of Incidence model to model angular dependent relection. The model can be one of "none" (no angularly dependent reflection), "front-cover" (simple refractive index), "anti-relect" (two layer front), and "user" (tabular data).
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

void InitConfigMask(simulation_config *C)
{
	int i;
	if ((C->sky_init)&&(C->topo_init)&&(C->loc_init))
	{
		double dt;
#ifdef OPENMP
		int totals=0, pc;
#endif		
		int pco=0;
		FreeConfigMask(C); // make sure we are clear to allocate new memory
		C->L=calloc(C->Nl,sizeof(location));
		if (C->L==NULL)
		{
			Warning("Could not allocate memory for locations\n");
			return;
		}
		printf("Tracing %d locations\n", C->Nl);
		TIC();
#pragma omp parallel private(i) shared(C,totals)
		{
#pragma omp for schedule(runtime)
			for (i=0;i<C->Nl;i++)
			{ 
				C->L[i]=ssdp_setup_location(&(C->S), &(C->T), C->albedo, C->o[i], C->x[i],C->y[i],C->z[i], &(C->M));
#ifdef OPENMP
				totals++;
				if (omp_get_thread_num()==0)
				{
					pc=100*(totals+1)/C->Nl;
					pco=ProgressBar(pc, pco, ProgressLen, ProgressTics);
				}
#else
				pco=ProgressBar((100*(i+1))/C->Nl, pco, ProgressLen, ProgressTics);
#endif
			}
		}
		pco=ProgressBar(100, pco, ProgressLen, ProgressTics);
		dt=TOC();
		if (!ssdp_error_state)
			printf("%d locations traced in %g s (%e s/horizons)\n", C->Nl, dt, dt/((double)C->Nl));
		else
			FreeConfigMask(C); // make sure we are clear to allocate new memory
			
	}
}

void InitConfigGridMask(simulation_config *C)
{
	int i;
	if ((C->sky_init)&&(C->grid_init)&&(C->loc_init))
	{
		double dt;
#ifdef OPENMP
		int totals=0, pc;
#endif		
		int pco=0;
		FreeConfigMask(C); // make sure we are clear to allocate new memory
		C->L=malloc(C->Nl*sizeof(location));
		if (C->L==NULL)
		{
			Warning("Could not allocate memory for locations\n");
			return;
		}
		printf("Tracing %d locations\n", C->Nl);
		TIC();
#pragma omp parallel private(i) shared(C,totals)
		{
#pragma omp for schedule(runtime)
			for (i=0;i<C->Nl;i++)
			{ 
				C->L[i]=ssdp_setup_grid_location(&(C->S), &(C->Tx), C->albedo, C->o[i], C->x[i],C->y[i],C->z[i], &(C->M));
#ifdef OPENMP
				totals++;
				if (omp_get_thread_num()==0)
				{
					pc=100*(totals+1)/C->Nl;
					pco=ProgressBar(pc, pco, ProgressLen, ProgressTics);
				}
#else
				pco=ProgressBar((100*(i+1))/C->Nl, pco, ProgressLen, ProgressTics);
#endif
			}
		}
		pco=ProgressBar(100, pco, ProgressLen, ProgressTics);
		dt=TOC();
		if (!ssdp_error_state)
			printf("%d locations traced in %g s (%e s/horizons)\n", C->Nl, dt, dt/((double)C->Nl));
		else
			FreeConfigMask(C); // make sure we are clear to allocate new memory
			
	}
}

void InitConfigMaskNoH(simulation_config *C) // same as above but without horizon calculation
{
	int i;
	if ((C->sky_init)&&(C->loc_init))
	{
		double dt;
		int pco=0;
		FreeConfigMask(C); // make sure we are clear to allocate new memory
		C->L=calloc(C->Nl,sizeof(location));
		if (C->L==NULL)
		{
			Warning("Could not allocate memory for locations\n");
			return;
		}
		printf("Tracing %d locations\n", C->Nl);
		TIC();
#pragma omp parallel private(i) shared(C)
		{
#ifdef OPENMP
			int nt=omp_get_num_threads();
#endif
#pragma omp for schedule(runtime)
			for (i=0;i<C->Nl;i++)
			{ 
				C->L[i]=ssdp_setup_location(&(C->S), NULL, C->albedo, C->o[i], C->x[i],C->y[i],C->z[i], &(C->M));
#ifdef OPENMP
				if (omp_get_thread_num()==0)
				{
					int pc;
					pc=100*(nt*i+1)/C->Nl;
					if (pc>100)
						pc=100;
					pco=ProgressBar(pc, pco, ProgressLen, ProgressTics);
				}
#else
				pco=ProgressBar((100*(i+1))/C->Nl, pco, ProgressLen, ProgressTics);
#endif
			}
		}
		pco=ProgressBar(100, pco, ProgressLen, ProgressTics);
		dt=TOC();
		if (!ssdp_error_state)
			printf("%d locations traced in %g s (%e s/horizons)\n", C->Nl, dt, dt/((double)C->Nl));
		else
			FreeConfigMask(C); // make sure we are clear to allocate new memory
			
	}
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
DESCRIPTION Setup the topogrid. Load the z data (column major, from the south-west corner to the north-east corner).
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
PARSEFLAG config_locations ConfigLoc "C=<out-config> x=<in-array> y=<in-array> z=<in-array> azimuth=<in-array> zenith=<in-array> [type=<topology/topogrid> [albedo=<float-value>]"
DESCRIPTION Setup the topography. Load the x, y, and z data of the unstructured topography mesh into the configuration data.
ARGUMENT x x coordinates
ARGUMENT y y coordinates
ARGUMENT z z coordinates
ARGUMENT type Optional argument to select the unstructured topology mesh or the structured topogrid. Valid values are \"topology\" or \"topogrid\" (default topology unless only a topogrid is defined)
ARGUMENT azimuth azimuth angle of tilted surface
ARGUMENT zenith zenith angle of tilted surface
ARGUMENT albedo optionally provide an albedo value between 0-1
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigLoc (char *in)
{
	simulation_config *C;
	array *x, *y, *z, *az, *ze;
	int i;
	char type='t';
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
	if (FetchArray(in, "azimuth", word, &az))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "zenith", word, &ze))
	{
		free(word);
		return;
	}
	if (GetOption(in, "type", word))
	{
		if (strncmp(word,"topology",10)==0)
		{
			if (C->topo_init==0)
			{ 
				Warning("No topology available\n");
				free(word);
				return;
			}
		}
		if (strncmp(word,"topogrid",10)==0)
		{
			if (C->grid_init==0)
			{ 
				Warning("No topogrid available\n");
				free(word);
				return;
			}
			type='g';
		}
	}
	if (C->topo_init==0)
	{
		type='g';
		if (C->grid_init==0)
			Warning("No topology or topogrid available\n");
	}
	if (FetchOptFloat(in, "albedo", word, &(C->albedo)))
		C->albedo=0.0;

	free(word);
	if (x->N!=y->N)
	{
		Warning("x- y- arrays must be of equal length\n"); 
		return;
	}
	if ((z->N!=x->N)&&(z->N!=1))
	{
		Warning("z- array must be either of equal length as x- and y- or of length 1\n"); 
		return;
	}
	
	if (az->N!=ze->N)
	{
		Warning("azimuth and zenith arrays must be of equal length\n"); 
		return;
	}
	if ((az->N!=x->N)&&(az->N!=1))
	{
		Warning("azimuth and zenith arrays must either have length 1 or be of equal length ans z,y,z\n"); 
		return;
	}
	
	FreeConfigLocation(C);
	printf("Configuring locations with %d points\n", x->N);
	
	C->Nl=x->N;
	C->x=malloc(C->Nl*sizeof(double));
	C->y=malloc(C->Nl*sizeof(double));
	C->z=malloc(C->Nl*sizeof(double));
	C->o=malloc(C->Nl*sizeof(sky_pos));
	if ((z->N==1)&&(az->N==1)) // one z, one orientation
	{
		for (i=0;i<C->Nl;i++)
		{
			C->x[i]=x->D[i];
			C->y[i]=y->D[i];
			C->z[i]=z->D[0];
			C->o[i].a=az->D[0];
			C->o[i].z=ze->D[0];
			
		}
	}
	else if (z->N==1) // many orientations, one z
	{
		for (i=0;i<C->Nl;i++)
		{
			C->x[i]=x->D[i];
			C->y[i]=y->D[i];
			C->z[i]=z->D[0];
			C->o[i].a=az->D[i];
			C->o[i].z=ze->D[i];
		}
	}
	else if (az->N==1) // many z's one orientation
	{
		for (i=0;i<C->Nl;i++)
		{
			C->x[i]=x->D[i];
			C->y[i]=y->D[i];
			C->z[i]=z->D[i];
			C->o[i].a=az->D[0];
			C->o[i].z=ze->D[0];
		}
	}
	else // just many
	{
		for (i=0;i<C->Nl;i++)
		{
			C->x[i]=x->D[i];
			C->y[i]=y->D[i];
			C->z[i]=z->D[i];
			C->o[i].a=az->D[i];
			C->o[i].z=ze->D[i];
		}
	}
	C->loc_init=1;
	
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
	}
	
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
DESCRIPTION Rotate template coordinates and place it at a given location. The given coordinates are rotated wrt to the origin in the provided coordinates. The units of the input coordinates are interpreted in the units of the selected epsg coordinate system. The output coordinates are given in WGS84 (epsg:4326) system.
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

        if (FetchOptFloat(in, "azi", word, &azi))
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
