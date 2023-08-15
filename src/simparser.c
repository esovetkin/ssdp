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


#define AT(arr, i) (arr->D[i % arr->N])


static int check_simconfig(simulation_config *C, int ok_notopo, int ok_nosky)
{
		if (ok_nosky && (!C->sky_init)) {
				Warning("Error: simulation config has no sky initialised\n");
				return -1;
		}

		if (!C->loc_init)
		{
				Warning("Error: simulation config has no locations initialised\n");
				return -2;
		}

		if ((!ok_notopo) && (!C->topo_init) && (!C->grid_init)) {
				Warning("Error: no topological data available\n");
				return -3;
		}

		if ((!C->topo_init) && (!C->grid_init))
		{
				Warning("Warning: no topological data available, omitting horizon\n");
				InitConfigMaskNoH(C);
		}

		if (ssdp_error_state) {
				ssdp_print_error_messages();
				ssdp_reset_errors();
				return -4;
		}

		return 0;
}


static int check_shapes(int n, array* arrs[])
{
		int i, N = -1;

		for (i=0; i < n; ++i)
				N = N > arrs[i]->N ? N : arrs[i]->N;

		for (i=0; i < n; ++i) {
				if (N != arrs[i]->N && 1 != arrs[i]->N) {
						Warning("Error: invalid length of input arrays\n");
						return -1;
				}
		}

		return N;
}


static int default_array(array **arr, double value)
{
		if (NULL==(*arr=malloc(sizeof(**arr))))
				goto earr;

		if (NULL==((*arr)->D = malloc(sizeof(*((*arr)->D)))))
				goto eD;

		(*arr)->N = 1;
		(*arr)->D[0] = value;
		return 1;
eD:
		free((*arr));
earr:
		return -1;
}


static void free_default_array(array **arr, int iffree)
{
		if (0 == iffree) return;

		free((*arr)->D);
		(*arr)->D = NULL;
		free((*arr));
		(*arr) = NULL;
}

/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG sim_static SimStatic "C=<in-config> t=<in-array> [p=<in-array>] [T=<in-array>] GHI=<in-array> DHI=<in-array> POA=<out-array>"
DESCRIPTION Computes the POA irradiance time series for all configured locations using the Perez All Weather Sky Model. 
ARGUMENT C Simulation config variable
ARGUMENT t unix time array.
ARGUMENT p air pressure in mb (default 1010)
ARGUMENT T air temperature in C (default 10)
ARGUMENT GHI global horizontal irradiance as a function of time
ARGUMENT DHI diffuse horizontal irradiance as a function of time
OUTPUT POA plane of array irradiance as a function of time and location (with n time values and m locations the array contains n times m values)
END_DESCRIPTION
*/
void SimStatic(char *in)
{
	int i, j; // loop through space and time
	char *word;
	simulation_config *C;
	array *t, *GH, *DH, out;
	array *p, *T, pp, TT;
	double tsky=0, tpoa=0;
	int pco=0;
	word=malloc((strlen(in)+1)*sizeof(char));
	pp.N=0;
	pp.D=NULL;
	TT.N=0;
	TT.D=NULL;
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (check_simconfig(C, 1, 0)) {
			free(word);
			return;
	}
	if (FetchArray(in, "t", word, &t))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "GHI", word, &GH))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "DHI", word, &DH))
	{
		free(word);
		return;
	}	
	if (FetchArray(in, "p", word, &p))
	{
		pp.N=1;
		pp.D=malloc(sizeof(double));
		pp.D[0]=1010.0;
		p=&pp;
	}	
	if (FetchArray(in, "T", word, &T))
	{
		TT.N=1;
		TT.D=malloc(sizeof(double));
		TT.D[0]=10.0;
		T=&TT;
	}	
	if ((p->N>1)&&(p->N!=t->N))
	{
		Warning("Length of p array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}	
	if ((T->N>1)&&(T->N!=t->N))
	{
		Warning("Length of T array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}
		
	if ((t->N!=GH->N)||(t->N!=DH->N))
	{
		Warning("Length of t-, GHI-, and DHI-arrays do not match\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}
	// fetch name of output var
	if (!GetArg(in, "POA", word))
	{
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}
	out.D=malloc(t->N*C->Nl*sizeof(double));
	if (out.D==NULL)
	{
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}	
	out.N=t->N*C->Nl;	
	printf("doing static simulation of %d locations at %d instances\n",C->Nl, t->N);
	
	for (j=0;j<t->N;j++)
	{
		TIC();
		ssdp_make_perez_all_weather_sky_coordinate(&(C->S), (time_t) t->D[j], C->lon, C->lat, C->E, p->D[j%p->N], T->D[j%T->N], GH->D[j], DH->D[j]);
		tsky+=TOC();
		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
			for (i=0;i<C->Nl;i++)
				out.D[j*C->Nl+i]=ssdp_total_poa(&(C->S), C->o[i], &(C->M), C->L+i);
		}
		pco=ProgressBar((100*(j+1))/(t->N), pco, ProgressLen, ProgressTics);
		tpoa+=TOC();
	}
	if (pp.D)
		free(pp.D);
	if (TT.D)
		free(TT.D);
	printf("Computed %d skies in %g s (%g s/sky)\n", t->N, tsky, tsky/((double)t->N));
	printf("Computed %d POA Irradiances in %g s (%g s/POA)\n", t->N*C->Nl, tpoa, tpoa/((double)(t->N*C->Nl)));
	printf("Creating array %s\n",word);
	if(AddArray(word, out))
	{
		free(word); // failed to make array
		free(out.D);
	}	
}

/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG sim_static_pos SimStaticPos "C=<in-config> sun_azimuth=<in-array> sun_zenith=<in-array> [dayofyear=<in-array>] GHI=<in-array> DHI=<in-array> POA=<out-array>"
DESCRIPTION Computes the POA irradiance time series for all configured locations using the Perez All Weather Sky Model. Here time and location is ignored and sky is computed based on the given sun position and dayofyear.
ARGUMENT C simulation config variable
ARGUMENT sun_azimuth,sun_zenith sun position in the sky (radians)
ARGUMENT dayofyear day of the year. len(dayofyear)=len(sun_azimuth) or 1 (default 0)
ARGUMENT GHI, DHI global and diffuse horizontal irradiance. len(GHI)=len(DHI)=len(sun_azimuth)=len(sun_zenith)
OUTPUT POA plane of array irradiance len(POA)=len(GHI)*len(locations)
END_DESCRIPTION
*/
void SimStaticPos(char *in)
{
		int i, j, N, pco = 0, fdoy = 0;
		simulation_config *C;
		char *word, *npoa;
		array *sazi, *szen, *doy, *GH, *DH, poa;
		double tsky = 0, tpoa = 0;
		sky_pos sun;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(npoa=malloc((strlen(in)+1)*sizeof(*npoa)))) goto enpoa;

		if (!GetArg(in, "POA", npoa)) goto eargs;
		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (check_simconfig(C, 1, 0)) goto eargs;
		if (FetchArray(in, "sun_azimuth", word, &sazi)) goto eargs;
		if (FetchArray(in, "sun_zenith", word, &szen)) goto eargs;
		if (FetchArray(in, "GHI", word, &GH)) goto eargs;
		if (FetchArray(in, "DHI", word, &DH)) goto eargs;
		if (FetchOptArray(in, "dayofyear", word, &doy))
				if ((fdoy=default_array(&doy, 0)) < 0) goto eargdoy;

		N = check_shapes(5,(array*[]){sazi, szen, GH, DH, doy});
		if (N < 0) goto eshapes;

		if (NULL==(poa.D=malloc(N*C->Nl*sizeof(*poa.D)))) goto epoa;
		poa.N = N*C->Nl;

		for (j=0; j < N; ++j) {
				TIC();
				sun.z = AT(szen,j);
				sun.a = AT(sazi,j);
				ssdp_make_perez_all_weather_sky
						(&(C->S), sun, AT(GH,j), AT(DH,j), AT(doy,j));
				tsky += TOC(); TIC();

#pragma omp parallel private (i) shared(C)
				{
#pragma omp for schedule(runtime)
						for (i=0; i < C->Nl; ++i)
								poa.D[j*C->Nl + i] = ssdp_total_poa
										(&(C->S), C->o[i], &(C->M), C->L+i);
				}
				pco=ProgressBar((100*(j+1))/N, pco, ProgressLen, ProgressTics);
				tpoa+=TOC();
		}

		if (AddArray(npoa, poa)) goto epoaadd;

		free_default_array(&doy, fdoy);
		free(word);
		return;
epoaadd:
		free(poa.D);
epoa:
eshapes:
		free_default_array(&doy, fdoy);
eargdoy:
eargs:
		free(npoa);
enpoa:
		free(word);
eword:
		Warning("Error: sim_static_pos failed\n");
		return;
}

/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG sim_static_integral SimStaticInt "C=<in-config> t=<in-array> [p=<in-array>] [T=<in-array>] GHI=<in-array> DHI=<in-array> POA=<out-array>"
DESCRIPTION Computes the integrated POA irradiance for all configured locations using the Perez All Weather Sky Model.
ARGUMENT C Simulation config variable
ARGUMENT t unix time array.
ARGUMENT p air pressure in mb (default 1010)
ARGUMENT T air temperature in C (default 10)
ARGUMENT GHI global horizontal irradiance as a function of time
ARGUMENT DHI diffuse horizontal irradiance as a function of time
OUTPUT POA plane of array irradiance as a function of time and location (with n time values and m locations the array contains n times m values)
END_DESCRIPTION
*/
void SimStaticInt(char *in)
{
	int i; // loop through space
	char *word;
	simulation_config *C;
	array *t, *GH, *DH, out;
	array *p, *T, pp, TT;
	double tsky, tpoa;
	word=malloc((strlen(in)+1)*sizeof(char));
	pp.N=0;
	pp.D=NULL;
	TT.N=0;
	TT.D=NULL;
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (check_simconfig(C, 1, 0)) {
			free(word);
			return;
	}
	if (FetchArray(in, "t", word, &t))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "GHI", word, &GH))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "DHI", word, &DH))
	{
		free(word);
		return;
	}
	if ((t->N!=GH->N)||(t->N!=DH->N))
	{
		Warning("Length of t-, GHI-, and DHI-arrays do not match\n");
		free(word);
		return;
	}
		
	if (FetchArray(in, "p", word, &p))
	{
		pp.N=t->N;
		pp.D=malloc(t->N*sizeof(double));
		for (i=0;i<t->N;i++)
			pp.D[i]=1010.0;
		p=&pp;
	}	
	if (FetchArray(in, "T", word, &T))
	{
		TT.N=t->N;
		TT.D=malloc(t->N*sizeof(double));
		for (i=0;i<t->N;i++)
			TT.D[i]=10.0;
		T=&TT;
	}	
	if ((p->N>1)&&(p->N!=t->N))
	{
		Warning("Length of p array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}	
	if ((T->N>1)&&(T->N!=t->N))
	{
		Warning("Length of T array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}	
	if ((p->N==1)&&(t->N>1))
	{
		pp.N=t->N;
		pp.D=malloc(t->N*sizeof(double));
		for (i=0;i<t->N;i++)
			pp.D[i]=1010.0;
		p=&pp;
	}
	if ((T->N==1)&&(t->N>1))
	{
		TT.N=t->N;
		TT.D=malloc(t->N*sizeof(double));
		for (i=0;i<t->N;i++)
			TT.D[i]=10.0;
		T=&TT;
	}
	// fetch name of output var
	if (!GetArg(in, "POA", word))
	{
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}
	out.D=calloc(C->Nl, sizeof(double));
	if (out.D==NULL)
	{
		free(word);
		return;
	}	
	out.N=C->Nl;	
	printf("doing static integral simulation of %d locations at %d instances\n",C->Nl, t->N);
	
	TIC();
	// TODO make p and T arrays with different lenths, adapt the integration routine and epand 
	// arrays to t length?
	ssdp_make_perez_cumulative_sky_coordinate(&(C->S),t->D, C->lon, C->lat, C->E, p->D, T->D, GH->D, DH->D, t->N);
	tsky=TOC();
	printf("Integrated %d skies in %g s (%g s/sky)\n", t->N, tsky, tsky/((double)t->N));
	
	TIC();
#pragma omp parallel private(i) shared(C)
	{
#pragma omp for schedule(runtime)
		for (i=0;i<C->Nl;i++)
			out.D[i]+=ssdp_total_poa(&(C->S), C->o[i], &(C->M), C->L+i);
	}
	tpoa=TOC();
	if (pp.D)
		free(pp.D);
	if (TT.D)
		free(TT.D);
	printf("Computed %d POA Irradiances in %g s (%g s/POA)\n", C->Nl, tpoa, tpoa/((double)(C->Nl)));
	printf("Creating array %s\n",word);
	if(AddArray(word, out))
	{
		free(word); // failed to make array
		free(out.D);
	}	
}
/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG sim_route SimRoute "C=<in-config> t=<in-array> [p=<in-array>] [T=<in-array>] GHI=<in-array> DHI=<in-array> POA=<out-array>"
DESCRIPTION Computes the POA irradiance along a route along all configured locations using the Perez All Weather Sky Model.
ARGUMENT C Simulation config variable
ARGUMENT t unix time array. As t progresses from t0-tn we move through locations l0-lm
ARGUMENT p air pressure in mb (default 1010)
ARGUMENT T air temperature in C (default 10)
ARGUMENT GHI global horizontal irradiance as a function of time
ARGUMENT DHI diffuse horizontal irradiance as a function of time
OUTPUT POA plane of array irradiance as a function of time
END_DESCRIPTION
*/
void SimRoute(char *in)
{
	int i, j; // loop through space and time
	char *word;
	simulation_config *C;
	array *t, *GH, *DH, out;
	array *p, *T, pp, TT;
	clock_t tsky0, tpoa0;
	clock_t tsky=0, tpoa=0;
	double ttsky, ttpoa;
	int pco=0;
	word=malloc((strlen(in)+1)*sizeof(char));
	pp.N=0;
	pp.D=NULL;
	TT.N=0;
	TT.D=NULL;
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (check_simconfig(C, 0, 0)) {
			free(word);
			return;
	}
	if (FetchArray(in, "t", word, &t))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "GHI", word, &GH))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "DHI", word, &DH))
	{
		free(word);
		return;
	}	
	
	if (FetchArray(in, "p", word, &p))
	{
		pp.N=1;
		pp.D=malloc(sizeof(double));
		pp.D[0]=1010.0;
		p=&pp;
	}	
	if (FetchArray(in, "T", word, &T))
	{
		TT.N=1;
		TT.D=malloc(sizeof(double));
		TT.D[0]=10.0;
		T=&TT;
	}	
	if ((p->N>1)&&(p->N!=t->N))
	{
		Warning("Length of p array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}	
	if ((T->N>1)&&(T->N!=t->N))
	{
		Warning("Length of T array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}	
	if ((t->N!=GH->N)||(t->N!=DH->N))
	{
		Warning("Length of t-, GHI-, and DHI-arrays do not match\n");
		free(word);
		return;
	}
	if (t->N<C->Nl)
		Warning("Warning: time array contains less points than there are waypoints\n");
	if (t->N>C->Nl)
		Warning("Warning: time array contains more points than there are waypoints\n");
	// fetch name of output var
	if (!GetArg(in, "POA", word))
	{
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}
	out.D=malloc(t->N*sizeof(double));
	if (out.D==NULL)
	{
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}	
	out.N=t->N;	
	printf("doing route simulation along %d locations at %d instances\n",C->Nl, t->N);
	for (j=0;j<t->N;j++)
	{
		// compute sky at evert time instance
		tsky0=clock();
		ssdp_make_perez_all_weather_sky_coordinate(&(C->S), (time_t) t->D[j], C->lon, C->lat, C->E, p->D[j%p->N], T->D[j%T->N],GH->D[j], DH->D[j]);
		tsky+=clock()-tsky0;
		tpoa0=clock();
		if (t->N>1)
			i=(j*(C->Nl-1))/(t->N-1);
		else
			i=0;
		out.D[j]=ssdp_total_poa(&(C->S), C->o[i], &(C->M), C->L+i);
		pco=ProgressBar((100*(j+1))/t->N, pco, ProgressLen, ProgressTics);
		tpoa+=clock()-tpoa0;
	}
	if (pp.D)
		free(pp.D);
	if (TT.D)
		free(TT.D);
	ttsky=(double)tsky/CLOCKS_PER_SEC;
	ttpoa=(double)tpoa/CLOCKS_PER_SEC;
	printf("Computed %d skies in %g s (%g s/sky)\n", t->N, ttsky, ttsky/((double)t->N));
	printf("Computed %d POA Irradiances in %g s (%g s/POA)\n", t->N, ttpoa, ttpoa/((double)(t->N)));
	printf("Creating array %s\n",word);
	if(AddArray(word, out))
	{
		free(word); // failed to make array
		free(out.D);
	}	
}



/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG sim_static_uniform SimStaticUniform "C=<in-config> t=<in-array> [p=<in-array>] [T=<in-array>] GHI=<in-array> DHI=<in-array> POA=<out-array>"
DESCRIPTION Computes the POA irradiance time series for all configured locations using a uniform sky. 
ARGUMENT C Simulation config variable
ARGUMENT t unix time array.
ARGUMENT p air pressure in mb (default 1010)
ARGUMENT T air temperature in C (default 10)
ARGUMENT GHI global horizontal irradiance as a function of time
ARGUMENT DHI diffuse horizontal irradiance as a function of time
OUTPUT POA plane of array irradiance as a function of time and location (with n time values and m locations the array contains n times m values)
END_DESCRIPTION
*/
void SimStaticUniform(char *in)
{
	int i, j; // loop through space and time
	char *word;
	simulation_config *C;
	array *t, *GH, *DH, out;
	array *p, *T, pp, TT;
	double tsky=0, tpoa=0;
	int pco=0;
	word=malloc((strlen(in)+1)*sizeof(char));
	pp.N=0;
	pp.D=NULL;
	TT.N=0;
	TT.D=NULL;
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (check_simconfig(C, 1, 0)) {
			free(word);
			return;
	}
	if (FetchArray(in, "t", word, &t))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "GHI", word, &GH))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "DHI", word, &DH))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "p", word, &p))
	{
		pp.N=1;
		pp.D=malloc(sizeof(double));
		pp.D[0]=1010.0;
		p=&pp;
	}	
	if (FetchArray(in, "T", word, &T))
	{
		TT.N=1;
		TT.D=malloc(sizeof(double));
		TT.D[0]=10.0;
		T=&TT;
	}	
	if ((p->N>1)&&(p->N!=t->N))
	{
		Warning("Length of p array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}	
	if ((T->N>1)&&(T->N!=t->N))
	{
		Warning("Length of T array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}		
	if ((t->N!=GH->N)||(t->N!=DH->N))
	{
		Warning("Length of t-, GHI-, and DHI-arrays do not match\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}
	// fetch name of output var
	if (!GetArg(in, "POA", word))
	{
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}
	out.D=malloc(t->N*C->Nl*sizeof(double));
	if (out.D==NULL)
	{
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}	
	out.N=t->N*C->Nl;	
	printf("doing static simulation of %d locations at %d instances\n",C->Nl, t->N);
	
	
	for (j=0;j<t->N;j++)
	{
		TIC();
		// only update sun position, we do not need to compute the diffuse sky, it is uniform and boring
		ssdp_make_skysunonly_coordinate(&(C->S), (time_t) t->D[j], C->lon, C->lat, C->E, p->D[j%p->N], T->D[j%T->N], GH->D[j], DH->D[j]);
		tsky+=TOC();
		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
			for (i=0;i<C->Nl;i++)
			{
				out.D[j*C->Nl+i]=ssdp_direct_poa(&(C->S), C->o[i], &(C->M), C->L+i);
				out.D[j*C->Nl+i]+=C->L[i].difftrans*DH->D[j]; // diffuse part is directly proportional to diffuse horizontal
			}
		}
		pco=ProgressBar((100*(j+1))/(t->N), pco, ProgressLen, ProgressTics);
		tpoa+=TOC();
	}
	if (pp.D)
		free(pp.D);
	if (TT.D)
		free(TT.D);
	printf("Computed %d skies in %g s (%g s/sky)\n", t->N, tsky, tsky/((double)t->N));
	printf("Computed %d POA Irradiances in %g s (%g s/POA)\n", t->N*C->Nl, tpoa, tpoa/((double)(t->N*C->Nl)));
	printf("Creating array %s\n",word);
	if(AddArray(word, out))
	{
		free(word); // failed to make array
		free(out.D);
	}	
}


/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG solpos SolarPos "t=<in-array> lon=<in-array> lat=<in-array> [E=<in-array>] [p=<in-array>] [T=<in-array>] azimuth=<out-array> zenith=<out-array>"
DESCRIPTION Computes the solar position using freespa [2]. Computation includes atmopheric refraction effects, (true solar position may be obtained by setting p=0).
ARGUMENT t unix time array
ARGUMENT lon longitude array (in degrees)
ARGUMENT lat latitude array (in degrees)
ARGUMENT E elevation (m) float (default 0)
ARGUMENT p pressure in mb array (default 1010)
ARGUMENT T temperature in C array (default 10)
OUTPUT azimuth sun azimuth (radians)
OUTPUT zenith sun zenith (radians)
END_DESCRIPTION
*/
void SolarPos(char *in)
{
		int i, N, fE = 0, fp = 0, fT = 0;
		char *word, *nazi, *nzen;
		sky_pos s;
		array *t, *lon, *lat, *E, *p, *T, azi, zen;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(nazi=malloc((strlen(in)+1)*sizeof(*nazi)))) goto enazi;
		if (NULL==(nzen=malloc((strlen(in)+1)*sizeof(*nzen)))) goto enzen;

		if (!GetArg(in, "azimuth", nazi)) goto earg;
		if (!GetArg(in, "zenith", nzen)) goto earg;
		if (FetchArray(in, "t", word, &t)) goto earg;
		if (FetchArray(in, "lon", word, &lon)) goto earg;
		if (FetchArray(in, "lat", word, &lat)) goto earg;

		if (FetchOptArray(in, "E", word, &E))
				if ((fE=default_array(&E, 1)) < 0) goto eargE;
		if (FetchOptArray(in, "p", word, &p))
				if ((fp=default_array(&p, 1010.0)) < 0) goto eargp;
		if (FetchOptArray(in, "T", word, &T))
				if ((fT=default_array(&T, 10.0)) < 0) goto eargT;

		N = check_shapes(6,(array*[]){t, lon, lat, E, p, T});
		if (N < 0) goto eshapes;

		if (NULL==(azi.D=malloc(N*sizeof(*azi.D)))) goto eazi;
		if (NULL==(zen.D=malloc(N*sizeof(*zen.D)))) goto ezen;
		azi.N = N; zen.N = N;

		printf("Computing the solar position at %d instances\n", N);
		TIC();
# pragma omp parallel private (i)
		{
#pragma omp for schedule(runtime)
				for (i=0; i<N; ++i) {
						s=ssdp_sunpos((time_t)AT(t,i),
									  deg2rad(AT(lat,i)), deg2rad(AT(lon,i)),
									  AT(E,i), AT(p,i), AT(T,i));
						azi.D[i]=s.a;
						zen.D[i]=s.z;
				}
		}
		printf("solpos: computed %d positions in %g s\n", N, TOC());

		if(AddArray(nzen, zen)) goto ezenadd;
		if(AddArray(nazi, azi)) goto ezen;

		free_default_array(&E, fE);
		free_default_array(&p, fp);
		free_default_array(&T, fT);
		free(word);
		return;
ezenadd:
		free(nzen); nzen = NULL;
		free(zen.D);
ezen:
		free(nazi); nazi = NULL;
		free(azi.D);
eazi:
eshapes:
		free_default_array(&T, fT);
eargT:
		free_default_array(&p, fp);
eargp:
		free_default_array(&E, fE);
eargE:
earg:
		free(nzen);
enzen:
		free(nazi);
enazi:
		free(word);
eword:
		Warning("Error: solpos failed!\n");
		return;
}

/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG suntimes SolarTimes "t=<in-array> lon=<float-value> lat=<float-value> [E=<float-value>] [p=<in-array>] [T=<in-array>] sunrise=<out-array> sunset=<out-array> transit=<out-array>"
DESCRIPTION Computes sunrise, transit and sunset times using freespa [2]. Computation includes atmopheric refraction effects and the effect of observer elevation. Note that the solar times are computed w.r.t. the solar transit on the specified UTC date. The sunrise and sunset times signify the period with sunlight. This means that in case of a polar night sunset and sunset coincide at the solar transit (when the sun's altitude is highest). For a midnight sun, sunrise and sunset are 24 hours apart centered around the solar transit with the highest altitude.
ARGUMENT t unix time array
ARGUMENT lon longitude float (in degrees)
ARGUMENT lat latitude float (in degrees)
ARGUMENT E elevation (m) float (default 0)
ARGUMENT p pressure in mb array (default 1010)
ARGUMENT T temperature in C array (default 10)
OUTPUT sunrise (unix time)
OUTPUT transit (unix time)
OUTPUT sunset  (unix time)
END_DESCRIPTION
*/
void SolarTimes(char *in)
{
	int i, r;
	char *word;
	sky_pos s;
	array *t, sunrise, sunset, transit, *p, *T;
	array pp, TT;
	double lon, lat, E;
	time_t t1, t2, t3;
	pp.N=0;
	pp.D=NULL;
	TT.N=0;
	TT.D=NULL;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchArray(in, "t", word, &t))
	{
		free(word);
		return;
	}
	
	if (FetchFloat(in, "lon", word, &lon))
	{
		free(word);
		return;
	}
	lon=deg2rad(lon);
	if (FetchFloat(in, "lat", word, &lat))
	{
		free(word);
		return;
	}
	lat=deg2rad(lat);
	
	if (FetchFloat(in, "E", word, &E))
		E=0;
		
	if (FetchArray(in, "p", word, &p))
	{
		// make a scalar array
		pp.N=1;
		pp.D=malloc(sizeof(double));
		pp.D[0]=1010.0;
		p=&pp;
	}
	if (FetchArray(in, "T", word, &T))
	{
		// make a scalar array
		TT.N=1;
		TT.D=malloc(sizeof(double));
		TT.D[0]=10.0;
		T=&TT;
	}
	
	if ((p->N>1)&&(p->N!=t->N))
	{
		Warning("The p-array must either have length one or be equal in length as the t-array\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}	
	if ((T->N>1)&&(T->N!=t->N))
	{
		Warning("The T-array must either have length one or be equal in length as the t-array\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}	
	
	sunrise.D=malloc(t->N*sizeof(double));
	if (sunrise.D==NULL)
	{
		Warning("Could not allocate memory for the sun rise time\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}		
	sunrise.N=t->N;
	sunset.D=malloc(t->N*sizeof(double));
	if (sunset.D==NULL)
	{
		Warning("Could not allocate memory for the sun set time\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}		
	sunset.N=t->N;
	transit.D=malloc(t->N*sizeof(double));
	if (transit.D==NULL)
	{
		Warning("Could not allocate memory for the sun transit time\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}		
	transit.N=t->N;
	
	
	printf("Computing solar times for %d instances\n",t->N);
	for (i=0;i<t->N;i++)
	{
		r=ssdp_suntimes((time_t)t->D[i], lat, lon, E, p->D[i%p->N], T->D[i%T->N],&t1,&t2,&t3);
		if (r&1)
			sunrise.D[i]=NAN; // is this wise as a signal there was no sunrise that day?
		else
			sunrise.D[i]=(double)t1;
		if (r&2)
			transit.D[i]=NAN; // this is probably a full blown error, how can there be no transit time?
		else
			transit.D[i]=(double)t2;
		if (r&4)
			sunset.D[i]=NAN;
		else
			sunset.D[i]=(double)t3;
	}
	if (pp.D)
		free(pp.D);
	if (TT.D)
		free(TT.D);
	
	if (!GetArg(in, "sunrise", word))
	{
		free(word);
		free(sunrise.D);
	}
	else
	{
		printf("Creating array %s\n",word);
		if(AddArray(word, sunrise))
		{
			free(word); // failed to make array
			free(sunrise.D);
			free(sunset.D);
			free(transit.D);
			return;
		}
		word=malloc((strlen(in)+1)*sizeof(char));
	}
	if (!GetArg(in, "transit", word))
	{
		free(word);
		free(transit.D);
	}
	else
	{
		printf("Creating array %s\n",word);
		if(AddArray(word, transit))
		{
			free(word); // failed to make array
			free(sunset.D);
			free(transit.D);
			return;
		}
		word=malloc((strlen(in)+1)*sizeof(char));
	}
	if (!GetArg(in, "sunset", word))
	{
		free(word);
		free(sunset.D);
	}
	else
	{
		printf("Creating array %s\n",word);
		if(AddArray(word, sunset))
		{
			free(word); // failed to make array
			free(sunset.D);
			return;
		}
	}	
}

/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG export_sky ExportSky "C=<in-config> t=<in-array> [p=<in-array>] [T=<in-array>] GHI=<in-array> DHI=<in-array> index=<int-value> file=<in-string>"
DESCRIPTION Exports a 3D polar plot of a sky according to the Perez All Weather Sky model. It only exports one location. You need to specify the index of the location (starting at index 0).
ARGUMENT C config-variable
ARGUMENT t single value unix time
ARGUMENT p pressure in mb array (default 1010)
ARGUMENT T temperature in C array (default 10)
ARGUMENT GHI single value global horizontal irradiance
ARGUMENT DHI single value diffuse horizontal irradiance
ARGUMENT index index of the location
ARGUMENT file filename
OUTPUT file A file with sky intensities. Organized in 4 columns: x y z I[W/sr]
END_DESCRIPTION
*/
void ExportSky(char *in)
{
	int j; 
	char *word;
	simulation_config *C;
	array *t, *GH, *DH;
	array *p, *T, pp, TT;
	word=malloc((strlen(in)+1)*sizeof(char));
	pp.N=0;
	pp.D=NULL;
	TT.N=0;
	TT.D=NULL;
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (check_simconfig(C, 1, 0)) {
			free(word);
			return;
	}
	if (FetchArray(in, "t", word, &t))
	{
		free(word);
		return;
	}
	if (t->N!=1)
	{
		Warning("Length of t array must be 1\n");
		free(word);
		return;
	}	
	if (FetchArray(in, "p", word, &p))
	{
		pp.N=1;
		pp.D=malloc(sizeof(double));
		pp.D[0]=1010.0;
		p=&pp;
	}	
	if (FetchArray(in, "T", word, &T))
	{
		TT.N=1;
		TT.D=malloc(sizeof(double));
		TT.D[0]=10.0;
		T=&TT;
	}	
	if ((p->N>1)&&(p->N!=t->N))
	{
		Warning("Length of p array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}	
	if ((T->N>1)&&(T->N!=t->N))
	{
		Warning("Length of T array must be 1 or equal to the length of array t\n");
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		free(word);
		return;
	}		
	if ((t->N!=GH->N)||(t->N!=DH->N))
	{
		Warning("Length of t-, GHI-, and DHI-arrays do not match\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}
	if (FetchArray(in, "GHI", word, &GH))
	{
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}
	if (GH->N!=1)
	{
		Warning("Length of GHI array must be 1\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}
	if (FetchArray(in, "DHI", word, &DH))
	{
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}
	if (DH->N!=1)
	{
		Warning("Length of DHI array must be 1\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}	
	if (FetchInt(in, "index", word, &j))
	{
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}		
	if (j<0)
	{
		Warning("index must be larger or equal to 0\n");
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}	
	if (j>=C->Nl)
	{
		Warning("index must be smaller than the number of locations (%d)\n", C->Nl);
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}	
	
	if (!GetArg(in, "file", word))
	{
		free(word);
		if (pp.D)
			free(pp.D);
		if (TT.D)
			free(TT.D);
		return;
	}
	// compute sky
	printf("Writing the sky to file %s as seen from location %d\n",word, j);
	ssdp_make_perez_all_weather_sky_coordinate(&(C->S), (time_t) t->D[0], C->lon, C->lat, C->E, p->D[0], T->D[0], GH->D[0], DH->D[0]);
	WriteDome4D(word, &(C->S), C->L+j);
	free(word);
}


/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG export_horizon ExportHorizon "C=<in-config> index=<int-value> file=<in-string>"
DESCRIPTION Exports the horizon at a certain location.
ARGUMENT C config-variable
ARGUMENT index index of the location
ARGUMENT file filename
OUTPUT file A file with the horizon. Organized in 2 columns: azimuth [rad] zenioth[rad]
END_DESCRIPTION
*/
void ExportHorizon(char *in)
{
	int j; 
	char *word;
	simulation_config *C;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (check_simconfig(C, 1, 1)) {
			free(word);
			return;
	}
	if (FetchInt(in, "index", word, &j))
	{
		free(word);
		return;
	}		
	if (j<0)
	{
		Warning("index must be larger or equal to 0\n");
		free(word);
		return;
	}	
	if (j>=C->Nl)
	{
		Warning("index must be smaller than the number of locations (%d)\n", C->Nl);
		free(word);
		return;
	}	
	
	if (!GetArg(in, "file", word))
	{
		free(word);
		return;
	}
	// export horizon
	printf("Writing the horizon to file %s as seen from location %d\n",word, j); 
	WriteHorizon(word, C->L+j);
	free(word);
}
