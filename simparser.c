#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
/* local includes */
#include "libssdp.h"
#include "io.h"
#include "util.h"
#include "variables.h"
#include "parser.h"
#include "parserutils.h"



// PARSEFLAG sim_static SimStatic "C=<config-variable> t=<array-variable> GHI=<array-variable> DHI=<array-variable> POA=<out-array>"
void SimStatic(char *in)
{
	int i, j; // loop through space and time
	char *word;
	simulation_config *C;
	array *t, *GH, *DH, out;
	clock_t tsky0, tpoa0;
	clock_t tsky=0, tpoa=0;
	double ttsky, ttpoa;
	int pco=0;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}		
	if (!C->sky_init)
	{		
		Warning("Simulation config has no sky initialized\n");
		free(word);
		return;
	}	
	if (!C->topo_init) 
	{	
		Warning("No topological data available, omitting horizon\n");
		InitConfigMaskNoH(C);
		if (ssdp_error_state)
		{
			ssdp_print_error_messages();
			ssdp_reset_errors();
			free(word);
			return;
		}
	}
	if (!C->loc_init) 
	{	
		Warning("Simulation config has no locations initialized\n");
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
	// fetch name of output var
	if (!GetArg(in, "POA", word))
	{
		free(word);
		return;
	}
	out.D=malloc(t->N*C->Nl*sizeof(double));
	if (out.D==NULL)
	{
		free(word);
		return;
	}	
	out.N=t->N*C->Nl;	
	for (j=0;j<t->N;j++)
	{
		// compute sky at evert time instance
		tsky0=clock();
		ssdp_make_perez_all_weather_sky_coordinate(&(C->S), (time_t) t->D[j], C->lon, C->lat, GH->D[j], DH->D[j]);
		tsky+=clock()-tsky0;
		tpoa0=clock();
		for (i=0;i<C->Nl;i++)
			out.D[j*C->Nl+i]=ssdp_total_poa(&(C->S), C->o[i], &(C->M), C->L+i);
		pco=ProgressBar((100*((j+1)*C->Nl))/(t->N*C->Nl), pco, ProgressLen, ProgressTics);
		tpoa+=clock()-tpoa0;
	}
	ttsky=(double)tsky/CLOCKS_PER_SEC;
	ttpoa=(double)tpoa/CLOCKS_PER_SEC;
	printf("Computed %d skies in %g s (%g s/sky)\n", t->N, ttsky, ttsky/((double)t->N));
	printf("Computed %d POA Irradiances in %g s (%g s/POA)\n", t->N*C->Nl, ttpoa, ttpoa/((double)(t->N*C->Nl)));
	if(AddArray(word, out))
	{
		free(word); // failed to make array
		free(out.D);
	}	
}
// PARSEFLAG sim_route SimRoute "C=<config-variable> t=<array-variable> GHI=<array-variable> DHI=<array-variable> POA=<out-array>"
void SimRoute(char *in)
{
	int i, j; // loop through space and time
	char *word;
	simulation_config *C;
	array *t, *GH, *DH, out;
	clock_t tsky0, tpoa0;
	clock_t tsky=0, tpoa=0;
	double ttsky, ttpoa;
	int pco=0;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}		
	if (!C->sky_init)
	{		
		Warning("Simulation config has no sky initialized\n");
		free(word);
		return;
	}	
	if (!C->topo_init) 
	{	
		Warning("No topological data available\n");
		free(word);
		return;
	}
	if (!C->loc_init) 
	{	
		Warning("Simulation config has no locations initialized\n");
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
	if (t->N<C->Nl)
		Warning("Warning: time array contains less points than there are waypoints\n");
	if (t->N>C->Nl)
		Warning("Warning: time array contains more points than there are waypoints\n");
	// fetch name of output var
	if (!GetArg(in, "POA", word))
	{
		free(word);
		return;
	}
	out.D=malloc(t->N*sizeof(double));
	if (out.D==NULL)
	{
		free(word);
		return;
	}	
	out.N=t->N;	
	for (j=0;j<t->N;j++)
	{
		// compute sky at evert time instance
		tsky0=clock();
		ssdp_make_perez_all_weather_sky_coordinate(&(C->S), (time_t) t->D[j], C->lon, C->lat, GH->D[j], DH->D[j]);
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
	ttsky=(double)tsky/CLOCKS_PER_SEC;
	ttpoa=(double)tpoa/CLOCKS_PER_SEC;
	printf("Computed %d skies in %g s (%g s/sky)\n", t->N, ttsky, ttsky/((double)t->N));
	printf("Computed %d POA Irradiances in %g s (%g s/POA)\n", t->N, ttpoa, ttpoa/((double)(t->N)));
	if(AddArray(word, out))
	{
		free(word); // failed to make array
		free(out.D);
	}	
}

// PARSEFLAG solpos SolarPos "t=<array-variable> lon=<longitude> lat=<latitude> azimuth=<output-array> zenith=<output-array>"
void SolarPos(char *in)
{
	int i;
	char *word;
	sky_pos s;
	array *t, azi, zen;
	double lon, lat;
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
	lon=degr2rad(lon);
	if (FetchFloat(in, "lat", word, &lat))
	{
		free(word);
		return;
	}
	lat=degr2rad(lat);
	azi.D=malloc(t->N*sizeof(double));
	if (azi.D==NULL)
	{
		free(word);
		return;
	}	
	azi.N=t->N;
	zen.D=malloc(t->N*sizeof(double));
	if (zen.D==NULL)
	{
		free(word);
		return;
	}	
	zen.N=t->N;
	for (i=0;i<t->N;i++)
	{
		s=ssdp_sunpos((time_t)t->D[i], lat, lon);
		azi.D[i]=s.a;
		zen.D[i]=s.z;
	}
	if (!GetArg(in, "azimuth", word))
	{
		free(word);
		return;
	}
	if(AddArray(word, azi))
	{
		free(word); // failed to make array
		free(azi.D);
	}	
	word=malloc((strlen(in)+1)*sizeof(char));
	if (!GetArg(in, "zenith", word))
	{
		free(word);
		return;
	}
	if(AddArray(word, zen))
	{
		free(word); // failed to make array
		free(zen.D);
	}	
}
