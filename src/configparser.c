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
PARSEFLAG config_coord ConfigCoord "C=<out-config> lat=<float> lon=<float>"
DESCRIPTION Setup the coordinate in the configuration variable.
ARGUMENT lat latitude (in radians)
ARGUMENT lon longitude (in radians)
OUTPUT C configuration variable
END_DESCRIPTION
*/
void ConfigCoord (char *in)
{
	simulation_config *C;
	array *l;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "lon", word, &l))
	{
		free(word);
		return;
	}
	if (l->N!=1)
	{
		Warning("config_coord expects a scalar value (array length 1) as longitude\n");
		free(word);
		return;
	}
	C->lon=l->D[0];
	printf("set longitude to %e degrees (%e rad)\n", rad2deg(C->lon), C->lon);
	if (FetchArray(in, "lat", word, &l))
	{
		free(word);
		return;
	}
	free(word);
	if (l->N!=1)
	{
		Warning("config_coord expects a scalar value (array length 1) as latitude\n");
		return;
	}
	C->lat=l->D[0];
	printf("set latitude to %e degrees (%e rad)\n", rad2deg(C->lat), C->lat);
}

/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_aoi ConfigAOI "C=<out-config> model=<none/front-cover/anti-reflect/user> [nf=<in-float> [nar=<in-float>]] [file=<in-file>]"
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
		for (i=0;C->Nl;i++)
			ssdp_free_location(C->L+i);
	}
	free(C->L);
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
#pragma omp parallel private(i) shared(C)
		{
#ifdef OPENMP
			int nt=omp_get_num_threads();
#endif
#pragma omp for 
			for (i=0;i<C->Nl;i++)
			{ 
				C->L[i]=ssdp_setup_location(&(C->S), &(C->T), C->albedo, C->o[i], C->x[i],C->y[i],C->z[i], &(C->M));
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

void InitConfigGridMask(simulation_config *C)
{
	int i;
	if ((C->sky_init)&&(C->grid_init)&&(C->loc_init))
	{
		double dt;
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
#pragma omp parallel private(i) shared(C)
		{
#ifdef OPENMP
			int nt=omp_get_num_threads();
#endif
#pragma omp for 
			for (i=0;i<C->Nl;i++)
			{ 
				C->L[i]=ssdp_setup_grid_location(&(C->S), &(C->Tx), C->albedo, C->o[i], C->x[i],C->y[i],C->z[i], &(C->M));
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
void InitConfigMaskNoH(simulation_config *C) // same as above but without horizon calculation
{
	int i;
	if ((C->sky_init)&&(C->loc_init))
	{
		double dt;
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
		for (i=0;i<C->Nl;i++)
		{
			C->L[i]=ssdp_setup_location(&(C->S), NULL, C->albedo, C->o[i], C->x[i],C->y[i],C->z[i], &(C->M));
			if (ssdp_error_state)
				break;
			pco=ProgressBar((100*(i+1))/C->Nl, pco, ProgressLen, ProgressTics);
		}
		dt=TOC();
		if (!ssdp_error_state)
			printf("%d locations traced in %g s (%g s/horizons)\n", C->Nl, dt, dt/((double)C->Nl));
		else
		{
			ssdp_print_error_messages();
			FreeConfigMask(C); // make sure we are clear to allocate new memory
		}
	}
}


/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG config_sky ConfigSKY "C=<out-config> N=<int>"
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
PARSEFLAG config_locations ConfigLoc "C=<out-config> x=<in-array> y=<in-array> z=<in-array> azimuth=<in-array> zenith=<in-array> [type=<topology/topogrid> [albedo=<in-float>]"
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
		{ 
			Warning("No topology or topogrid available\n");
			free(word);
			return;
		}
	}
	if (GetOption(in, "albedo", word))
	{
		C->albedo=atof(word);
	}
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
	
	if (type=='t')
		InitConfigMask(C);
	else
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
