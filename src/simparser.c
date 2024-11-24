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
#include "h5io.h"


static int check_simconfig(simulation_config *C, int ok_notopo, int ok_nosky)
{
		if ((!ok_nosky) && (!C->sky_init)) {
				Warning("Error: simulation config has no sky initialised\n");
				return -1;
		}

		if (!C->loc_init)
		{
				Warning("Error: simulation config has no locations initialised\n");
				return -2;
		}

		if ((!ok_notopo) && (!C->topo_init) && (!C->grid_init)) {
				Warning("Error: no topographical data available\n");
				return -3;
		}

		if (ssdp_error_state) {
				ssdp_print_error_messages();
				ssdp_reset_errors();
				return -4;
		}

		return 0;
}


static int poa_total(simulation_config *C, int o, double* out, int chunkid)
{
		if (InitLocations(C, chunkid, 1)) return -1;
		int i;

#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->Nl_eff; ++i)
						out[o+C->Nl_o+i] = ssdp_total_poa(
								C->S, C->o[C->Nl_o+i], &(C->M), C->L+i);
		}

		return 0;
}


static int poa_unif(simulation_config *C, int o, double* out,
					int chunkid, double DHI)
{
		if (InitLocations(C, chunkid, 1)) return -1;
		int i;

#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0; i < C->Nl_eff; ++i) {
						out[o+C->Nl_o+i] = ssdp_direct_poa
								(C->S, C->o[C->Nl_o+i], &(C->M), C->L+i);
// diffuse part is directly proportional to diffuse horizontal
						out[o+C->Nl_o+i] = C->L[i].difftrans * DHI;

				}
		}

		return 0;
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
		int i, j, N, pco=0, fp = 0, fT = 0;
		char *word, *nout;
		simulation_config *C;
		array *t, *GH, *DH, *p, *T, out;

		double tsky=0, tpoa=0;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(nout=malloc((strlen(in)+1)*sizeof(*nout)))) goto enout;

		if (!GetArg(in, "POA", nout)) goto eargs;
		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (check_simconfig(C, 1, 0)) goto eargs;
		if (FetchArray(in, "t", word, &t)) goto eargs;
		if (FetchArray(in, "GHI", word, &GH)) goto eargs;
		if (FetchArray(in, "DHI", word, &DH)) goto eargs;
		if (FetchOptArray(in, "p", word, &p))
				if ((fp=default_array(&p, 1010.0)) < 0) goto eargsp;
		if (FetchOptArray(in, "T", word, &T))
				if ((fT=default_array(&T, 10.0)) < 0) goto eargsT;

		N = check_shapes(5,(array*[]){t, GH, DH, p, T});
		if (N < 0) goto eshapes;

		if (NULL==(out.D=malloc(N*C->Nl*sizeof(*(out.D))))) goto eout;
		out.N=N*C->Nl;
		printf("doing static simulation of %d locations at %d instances\n", C->Nl, N);

		for (j=0; j < N; ++j) {
				TIC();
				ssdp_make_perez_all_weather_sky_coordinate
						(C->S, (time_t) AT(t,j),
						 C->lon, C->lat, C->E,
						 AT(p,j), AT(T,j), AT(GH,j), AT(DH,j));
				tsky+=TOC(); TIC();

				i = -1;
				while (NextChunk(C, &i))
						if (poa_total(C, j*C->Nl, out.D, i))
								goto esim;
				ProgressBar((100*(j+1))/N, &pco, ProgressLen, ProgressTics);
				tpoa+=TOC();
		}

		printf("Computed %d skies in %g s (%g s/sky)\n",
			   N, tsky, tsky/((double)N));
		printf("Computed %d POA Irradiances in %g s (%g s/POA)\n",
			   N*C->Nl, tpoa, tpoa/((double)(N*C->Nl)));

		printf("Creating array %s\n", nout);
		if(AddArray(nout, out)) goto eoutadd;

		free_default_array(&p, fp);
		free_default_array(&T, fT);
		free(word);
		return;
eoutadd:
esim:
		free(out.D);
eout:
eshapes:
		free_default_array(&T, fT);
eargsT:
		free_default_array(&p, fp);
eargsp:
eargs:
		free(nout);
enout:
		free(word);
eword:
		Warning("Error: sim_static failed!\n");
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
						(C->S, sun, AT(GH,j), AT(DH,j), AT(doy,j));
				tsky += TOC(); TIC();

				i = -1;
				while (NextChunk(C, &i))
						if (poa_total(C, j*C->Nl, poa.D, i))
								goto esim;

				ProgressBar((100*(j+1))/N, &pco, ProgressLen, ProgressTics);
				tpoa+=TOC();
		}

		printf("Computed %d skies in %g s (%g s/sky)\n",
			   N, tsky, tsky/((double)N));
		printf("Computed %d POA Irradiances in %g s (%g s/POA)\n",
			   N*C->Nl, tpoa, tpoa/((double)(N*C->Nl)));

		printf("Creating array %s\n", npoa);
		if (AddArray(npoa, poa)) goto epoaadd;

		free_default_array(&doy, fdoy);
		free(word);
		return;
epoaadd:
esim:
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
		int i, N, fp = 0, fT = 0;
		char *word, *nout;
		simulation_config *C;
		array *t, *GH, *DH, *p, *T, out;
		double tsky = 0, tpoa = 0;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(nout=malloc((strlen(in)+1)*sizeof(*nout)))) goto enout;

		if (!GetArg(in, "POA", nout)) goto eargs;
		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (check_simconfig(C, 1, 0)) goto eargs;
		if (FetchArray(in, "t", word, &t)) goto eargs;
		if (FetchArray(in, "GHI", word, &GH)) goto eargs;
		if (FetchArray(in, "DHI", word, &DH)) goto eargs;
		if (FetchOptArray(in, "p", word, &p))
				if ((fp=default_array(&p, 1010.0)) < 0) goto eargsp;
		if (FetchOptArray(in, "T", word, &T))
				if ((fT=default_array(&T, 10.0)) < 0) goto eargsT;

		N = check_shapes(5,(array*[]){t, GH, DH, p, T});
		if (N < 0) goto eshapes;

		if (N != t->N || t->N != GH->N || t->N != DH->N) {
				Warning("Error: t, GHI, DHI must be of the same length!\n");
				goto eshapes;
		}

		if (NULL==(out.D=malloc(C->Nl*sizeof(*(out.D))))) goto eout;
		out.N=C->Nl;

		TIC();
		printf("doing static integral simulation of %d locations at %d instances\n",C->Nl, N);
		ssdp_make_perez_cumulative_sky_coordinate
				(C->S, t->D, C->lon, C->lat, C->E,
				 p->D, p->N, T->D, T->N, GH->D, DH->D, t->N);
		tsky=TOC();
		printf("Integrated %d skies in %g s (%g s/sky)\n", t->N, tsky, tsky/((double)t->N));

		TIC();
		i = -1;
		while (NextChunk(C, &i))
				if (poa_total(C, 0, out.D, i))
						goto esim;
		tpoa=TOC();

		printf("Computed %d POA Irradiances in %g s (%g s/POA)\n", C->Nl, tpoa, tpoa/((double)(C->Nl)));

		printf("Creating array %s\n",nout);
		if(AddArray(nout, out)) goto eoutadd;

		free_default_array(&p, fp);
		free_default_array(&T, fT);
		free(word);
		return;
eoutadd:
esim:
		free(out.D);
eout:
eshapes:
		free_default_array(&T, fT);
eargsT:
		free_default_array(&p, fp);
eargsp:
eargs:
		free(nout);
enout:
		free(word);
eword:
		Warning("Error: sim_static_integral failed!\n");
}


/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG sim_sky SimSky "C=<in-config> POA=<out-array>"
DESCRIPTION Do simulations with the configured sky
ARGUMENT C Simulation config variable
OUTPUT POA plane of array irradiance for each configured location
END_DESCRIPTION
*/
void SimSky(char *in)
{
		int i;
		char *word, *nout;
		simulation_config *C;
		double tpoa = 0;
		array out;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(nout=malloc((strlen(in)+1)*sizeof(*nout)))) goto enout;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (!GetArg(in, "POA", nout)) goto eargs;
		if (check_simconfig(C, 1, 0)) goto eargs;

		if (NULL==(out.D=malloc(C->Nl*sizeof(*(out.D))))) goto eout;
		out.N=C->Nl;

		TIC();
		i = -1;
		while (NextChunk(C, &i))
				if (poa_total(C, 0, out.D, i))
						goto esim;
		tpoa=TOC();

		printf("Computed %d POA Irradiances in %g s (%g s/POA)\n",
			   C->Nl, tpoa, tpoa/((double)(C->Nl)));
		printf("Creating array %s\n",nout);
		if(AddArray(nout, out)) goto eoutadd;

		free(word);
		return;
eoutadd:
esim:
		free(out.D);
eout:
eargs:
		free(nout);
enout:
		free(word);
eword:
		Warning("Error: sim_sky failed!\n");
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
		int i, j, N, pco=0, fp=0, fT=0;
		char *word, *nout;
		simulation_config *C;
		array *t, *GH, *DH, *p, *T, out;
		double dt;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(nout=malloc((strlen(in)+1)*sizeof(*nout)))) goto enout;

		if (!GetArg(in, "POA", nout)) goto eargs;
		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (check_simconfig(C, 1, 0)) goto eargs;
		if (FetchArray(in, "t", word, &t)) goto eargs;
		if (FetchArray(in, "GHI", word, &GH)) goto eargs;
		if (FetchArray(in, "DHI", word, &DH)) goto eargs;
		if (FetchOptArray(in, "p", word, &p))
				if ((fp=default_array(&p,1010.0)) < 0) goto eargsp;
		if (FetchOptArray(in, "T", word, &T))
				if ((fT=default_array(&T,10.0)) < 0) goto eargsT;

		N = check_shapes(5, (array*[]){t, GH, DH, p, T});
		if (N < 0) goto eshapes;
		if (N<C->Nl)
				Warning("Warning: time array contains less points than there are waypoints\n");
		if (N>C->Nl) goto eshapes;

		if (NULL==(out.D=malloc(C->Nl*sizeof(*out.D)))) goto eout;
		out.N=C->Nl;
		printf("doing route simulation along %d locations at %d instances\n",C->Nl, N);

		TIC();
		i = -1;
		while (NextChunk(C, &i)) {
				if (InitLocations(C, i, 1)) goto esim;

#pragma omp parallel private(j)
				{
#pragma omp for schedule(runtime)
						for (j=0; j < C->Nl_eff; ++j) {
#ifdef OPENMP
								int thread=omp_get_thread_num();
#else
								int thread=0;
#endif
								ssdp_make_perez_all_weather_sky_coordinate
										(C->S+thread, (time_t) AT(t,C->Nl_o+j),
										 C->lon, C->lat, C->E,
										 AT(p,C->Nl_o+j), AT(T,C->Nl_o+j),
										 AT(GH,C->Nl_o+j), AT(DH,C->Nl_o+j));

								out.D[C->Nl_o+j]=ssdp_total_poa(
										C->S+thread, C->o[C->Nl_o+j],
										&(C->M), C->L+j);
								ProgressBar((100*(C->Nl_o+j+1))/C->Nl_eff, &pco, ProgressLen, ProgressTics);
						}
				}
		}
		ProgressBar(100, &pco, ProgressLen, ProgressTics);

		dt = TOC();
		printf("Computed %d skies in %g s (%g s/sky)\n", N, dt, dt/((double)N));

		printf("Creating array %s\n", nout);
		if(AddArray(nout, out)) goto eoutadd;

		free_default_array(&p, fp);
		free_default_array(&T, fT);
		free(word);
		return;
eoutadd:
esim:
		free(out.D);
eout:
eshapes:
		free_default_array(&T, fT);
eargsT:
		free_default_array(&p, fp);
eargsp:
eargs:
		free(nout);
enout:
		free(word);
eword:
		Warning("Error: sim_route failed!\n");
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
		int i, j, N, pco=0, fp = 0, fT = 0;
		char *word, *nout;
		simulation_config *C;
		array *t, *GH, *DH, *p, *T, out;
		double tsky=0, tpoa=0;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(nout=malloc((strlen(in)+1)*sizeof(*nout)))) goto enout;

		if (!GetArg(in, "POA", nout)) goto eargs;
		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (check_simconfig(C, 1, 0)) goto eargs;
		if (FetchArray(in, "t", word, &t)) goto eargs;
		if (FetchArray(in, "GHI", word, &GH)) goto eargs;
		if (FetchArray(in, "DHI", word, &DH)) goto eargs;
		if (FetchOptArray(in, "p", word, &p))
				if ((fp=default_array(&p, 1010.0)) < 0) goto eargsp;
		if (FetchOptArray(in, "T", word, &T))
				if ((fT=default_array(&T, 10.0)) < 0) goto eargsT;

		N = check_shapes(5,(array*[]){t, GH, DH, p, T});
		if (N < 0) goto eshapes;

		if (NULL==(out.D=malloc(N*C->Nl*sizeof(*(out.D))))) goto eout;
		out.N=N*C->Nl;
		printf("doing static simulation of %d locations at %d instances\n",C->Nl, N);

		for (j=0; j < N; ++j) {
				TIC();
				// only update sun position, we do not need to compute
				// the diffuse sky, it is uniform and boring
				ssdp_make_skysunonly_coordinate
						(C->S, (time_t) AT(t,j),
						 C->lon, C->lat, C->E,
						 AT(p,j), AT(T,j), AT(GH,j), AT(DH,j));
				tsky+=TOC(); TIC();

				i = -1;
				while (NextChunk(C, &i))
						if (poa_unif(C, j*C->Nl, out.D, i, DH->D[j]))
								goto esim;
				ProgressBar((100*(j+1))/N, &pco, ProgressLen, ProgressTics);
				tpoa+=TOC();
		}

		printf("Computed %d skies in %g s (%g s/sky)\n",
			   N, tsky, tsky/((double)N));
		printf("Computed %d POA Irradiances in %g s (%g s/POA)\n",
			   N*C->Nl, tpoa, tpoa/((double)(N*C->Nl)));

		printf("Creating array %s\n", nout);
		if(AddArray(nout, out)) goto eoutadd;

		free_default_array(&p, fp);
		free_default_array(&T, fT);
		free(word);
		return;
eoutadd:
esim:
		free(out.D);
eout:
eshapes:
		free_default_array(&T, fT);
eargsT:
		free_default_array(&p, fp);
eargsp:
eargs:
		free(nout);
enout:
		free(word);
eword:
		Warning("Error: sim_static_uniform failed!\n");
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
PARSEFLAG suntimes SolarTimes "t=<in-array> lon=<in-array> lat=<in-array> [E=<in-array>] [p=<in-array>] [T=<in-array>] sunrise=<out-array> sunset=<out-array> transit=<out-array>"
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
		int i, r, N, fE = 0, fp = 0, fT = 0;
		char *word, *nsunrise, *nsunset, *ntransit;
		array *t, *p, *T, *E, *lon, *lat, sunrise, sunset, transit;
		time_t t1, t2, t3;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(nsunrise=malloc((strlen(in)+1)*sizeof(*nsunrise)))) goto enrise;
		if (NULL==(nsunset=malloc((strlen(in)+1)*sizeof(*nsunset)))) goto enset;
		if (NULL==(ntransit=malloc((strlen(in)+1)*sizeof(*ntransit)))) goto entran;

		if (!GetArg(in, "sunrise", nsunrise)) goto eargs;
		if (!GetArg(in, "sunset", nsunset)) goto eargs;
		if (!GetArg(in, "transit", ntransit)) goto eargs;
		if (FetchArray(in, "t", word, &t)) goto eargs;
		if (FetchArray(in, "lon", word, &lon)) goto eargs;
		if (FetchArray(in, "lat", word, &lat)) goto eargs;
		if (FetchOptArray(in, "E", word, &E))
				if ((fE=default_array(&E, 0.0)) < 0) goto eargsE;
		if (FetchOptArray(in, "p", word, &p))
				if ((fp=default_array(&p, 1010.0)) < 0) goto eargsp;
		if (FetchOptArray(in, "T", word, &T))
				if ((fT=default_array(&T, 10.0)) < 0) goto eargsT;

		N = check_shapes(6,(array*[]){t,lon,lat,E,p,T});
		if (N < 0) goto eshapes;

		if (NULL==(sunrise.D=malloc(N*sizeof(*(sunrise.D))))) goto erise;
		sunrise.N = N;
		if (NULL==(sunset.D=malloc(N*sizeof(*(sunset.D))))) goto eset;
		sunset.N = N;
		if (NULL==(transit.D=malloc(N*sizeof(*(transit.D))))) goto etran;
		transit.N = N;

		printf("Computing solar times for %d instances\n", N);
		for (i=0; i < N; ++i) {
				r = ssdp_suntimes(
						(time_t)AT(t,i),
						deg2rad(AT(lat,i)), deg2rad(AT(lon,i)), AT(E, i),
						AT(p,i), AT(T,i), &t1, &t2, &t3);
				if (r&1)
						// TODO is this wise as a signal there was no sunrise that day?
						sunrise.D[i]=NAN;
				else
						sunrise.D[i]=(double)t1;
				if (r&2)
						// TODO this is probably a full blown error, how can there be no transit time?
						transit.D[i]=NAN;
				else
						transit.D[i]=(double)t2;
				if (r&4)
						sunset.D[i]=NAN;
				else
						sunset.D[i]=(double)t3;
		}

		printf("Creating array %s\n", nsunrise);
		if(AddArray(ntransit, transit)) goto eaddtran;
		if(AddArray(nsunset, sunset)) goto etran;
		if(AddArray(nsunrise, sunrise)) goto eset;

		free_default_array(&E, fE);
		free_default_array(&p, fp);
		free_default_array(&T, fT);
		free(word);
		return;
eaddtran:
		free(transit.D);
		free(ntransit); ntransit=NULL;
etran:
		free(sunset.D);
		free(nsunset); nsunset=NULL;
eset:
		free(sunrise.D);
		free(nsunrise); nsunrise=NULL;
erise:
eshapes:
		free_default_array(&T, fT);
eargsT:
		free_default_array(&p, fp);
eargsp:
		free_default_array(&E, fE);
eargsE:
eargs:
		free(ntransit);
entran:
		free(nsunset);
enset:
		free(nsunrise);
enrise:
		free(word);
eword:
		Warning("Error: suntimes failed!\n");
}


/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG export_sky ExportSky "C=<in-config> t=<float-value> [p=<float-value>] [T=<float-value>] GHI=<float-value> DHI=<float-value> index=<int-value> file=<in-string>"
DESCRIPTION Exports a 3D polar plot of a sky according to the Perez All Weather Sky model. It only exports one location. You need to specify the index of the location (starting at index 0).
ARGUMENT C config-variable
ARGUMENT t unix time
ARGUMENT p pressure in mb array (default 1010)
ARGUMENT T temperature in C array (default 10)
ARGUMENT GHI global horizontal irradiance
ARGUMENT DHI diffuse horizontal irradiance
ARGUMENT index index of the location
ARGUMENT file filename
OUTPUT file A file with sky intensities. Organized in 4 columns: x y z I[W/sr]
END_DESCRIPTION
*/
void ExportSky(char *in)
{
		int j;
		char *word, *file;
		simulation_config *C;
		double p, T, t, GH, DH;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(file=malloc((strlen(in)+1)*sizeof(*file)))) goto efile;

		if (!GetArg(in, "file", file)) goto eargs;
		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (check_simconfig(C, 1, 0)) goto eargs;
		if (FetchFloat(in, "t", word, &t)) goto eargs;
		if (FetchFloat(in, "GHI", word, &GH)) goto eargs;
		if (FetchFloat(in, "DHI", word, &DH)) goto eargs;
		if (FetchInt(in, "index", word, &j)) goto eargs;
		if (FetchOptFloat(in, "p", word, &p)) p = 1010.0;
		if (FetchOptFloat(in, "T", word, &T)) T = 10.0;

		if (j<0) {
				Warning("Error: index must be larger or equal to 0\n");
				goto ej;
		}
		if (j>=C->Nl) {
				Warning("Error: index must be smaller than the number of locations (%d)\n", C->Nl);
				goto ej;
		}

		// compute sky
		printf("Writing the sky to file %s as seen from location %d\n", word, j);
		ssdp_make_perez_all_weather_sky_coordinate(
				C->S, (time_t) t, C->lon, C->lat, C->E, p, T, GH, DH);
		WriteDome4D(file, C->S, C->L+j);

		free(file);
		free(word);
		return;
ej:
eargs:
		free(file);
efile:
		free(word);
eword:
		Warning("Error: export_sky failed!\n");
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
		char *word, *file;
		simulation_config *C;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(file=malloc((strlen(in)+1)*sizeof(*file)))) goto efile;

		if (!GetArg(in, "file", file)) goto eargs;
		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (check_simconfig(C, 1, 1)) goto eargs;
		if (FetchInt(in, "index", word, &j)) goto eargs;
		if (j<0) {
				Warning("Error: index must be larger or equal to 0\n");
				goto eargs;
		}
		if (j>=C->Nl) {
				Warning("Error: index must be smaller than the number of locations (%d)\n", C->Nl);
				goto eargs;
		}

		// export horizon
		printf("Writing the horizon to file %s as seen from location %d\n", file, j);
		WriteHorizon(file, C->L+j);

		free(file);
		free(word);
		return;
eargs:
		free(file);
efile:
		free(word);
eword:
		Warning("Error: export_horizon failed!\n");
}

/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG compute_sky ComputeSky "C=<config_file> t=<in-array> GHI=<in-array> DHI=<in-array> [p=<in-array>] [T=<in-array>]"
DESCRIPTION Compute sky. Currently only implemented integral sky
ARGUMENT C Simulation configuration
ARGUMENT t unix time array.
ARGUMENT GHI global horizontal irradiance as a function of time
ARGUMENT DHI diffuse horizontal irradiance as a function of time
ARGUMENT p air pressure in mb (default: 1010)
ARGUMENT T air temperature in C (default: 10)
ARGUMENT type type of sky to compute (default: integral)
OUTPUT C the computed sky is saved in the simulation configuration
END_DESCRIPTION
*/
void ComputeSky(char *in)
{
		int N, fp = 0, fT = 0;
		char *word;
		simulation_config *C;
		array *t, *GH, *DH, *p, *T;
		double tsky=0;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		// if (check_simconfig(C, 1, 0)) goto eargs;
		if (FetchArray(in, "t", word, &t)) goto eargs;
		if (FetchArray(in, "GHI", word, &GH)) goto eargs;
		if (FetchArray(in, "DHI", word, &DH)) goto eargs;
		if (FetchOptArray(in, "p", word, &p))
				if ((fp=default_array(&p, 1010.0)) < 0) goto eargsp;
		if (FetchOptArray(in, "T", word, &T))
				if ((fT=default_array(&T, 10.0)) < 0) goto eargsT;

		N = check_shapes(5,(array*[]){t, GH, DH, p, T});
		if (N < 0) goto eshapes;

		if (N != t->N || t->N != GH->N || t->N != DH->N) {
				Warning("Error: t, GHI, DHI must be of the same length!\n");
				goto eshapes;
		}

		printf("Computing perez integrated sky...\n");
		TIC();
		ssdp_make_perez_cumulative_sky_coordinate
				(C->S, t->D, C->lon, C->lat, C->E,
				 p->D, p->N, T->D, T->N, GH->D, DH->D, t->N);
		tsky=TOC();
		printf("Integrated %d skies in %g s (%g s/sky)\n", t->N, tsky, tsky/((double)t->N));

		free(word);
		return;
eshapes:
		free_default_array(&T, fT);
eargsT:
		free_default_array(&p, fp);
eargsp:
eargs:
		free(word);
eword:
		Warning("Error: compute_sky failed!\n");
}


/*
BEGIN_DESCRIPTION
SECTION Simulation
PARSEFLAG write_sky WriteSky "C=<config_file> file=<in-string> [dataset=<in-string>]"
DESCRIPTION save currently configured sky in h5 file
ARGUMENT C Simulation config variable
ARGUMENT file name of the output file
ARGUMENT dataset optional name of dataset (default "skydata")
OUTPUT file the sky is save uncompressed in a structured h5file
END_DESCRIPTION
*/
void WriteSky(char *in)
{
		char *fn, *word, *dst;
		simulation_config *C;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(fn=malloc((strlen(in)+1)*sizeof(*fn)))) goto efn;
		if (NULL==(dst=malloc((strlen(in)+1)*sizeof(*dst)))) goto edst;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (!GetArg(in, "file", fn)) goto eargs;
		if (!GetOption(in, "dataset", dst)) snprintf(dst, strlen(in), "skydata");

		TIC();
		if (ssdp_write_sky(C->S, fn, dst)) {
				printf("Error: failed to write %s to %s\n", dst, fn);
				goto ewrite;
		}
		printf("Wrote sky to %s::%s in %g s\n", fn, dst, TOC());

		free(dst);
		free(fn);
		free(word);
		return;
ewrite:
eargs:
		free(dst);
edst:
		free(fn);
efn:
		free(word);
eword:
		printf("Error: write_sky failed!\n");
		return;
}


/*
BEGIN_DESCRIPTION
SECTION Simulation Configuration
PARSEFLAG write_horizon WriteHoriz "C=<in-config> file=<file-str> [dataset=<str>] [gzip=<int-value>] [chunksize=<int-value>] [cachemb=<int-value>] [cacheslots=<int-value>]"
DESCRIPTION Write computed horizon for provided locations to a HDF5 file. For each location horizon is stored as one array (see write_h5), setting chunksize might be beneficial for compression ratios.
ARGUMENT C simulation configuration with configured locations
ARGUMENT file name of the output file
ARGUMENT dataset optional name of dataset (default "horizon")
ARGUMENT gzip level of compression to use (default "0")
ARGUMENT chunksize optional number of horizons kept in one chunk (default "10000"). If 0 the array is written contigiously and no compression can be used.
ARGUMENT cachemb,cacheslots optional size of h5 cache (default: cachemb=64, cacheslots=12421).
OUTPUT file output filename
END_DESCRIPTION
*/
void WriteHoriz(char *in)
{
		double **data = NULL;
		int i = 0, chunkid = 0, narr = 0, arrlen = 0;
		char *word = NULL;
		simulation_config *C = NULL;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (FetchConfig(in, "C", word, &C)) goto eC;

		if (NULL == (data=malloc(C->Nl*sizeof(*data)))) goto edata;
		narr = C->Nl;

		chunkid = -1;
		while (NextChunk(C, &chunkid)) {
				if (InitLocations(C, chunkid, 0)) goto einitlocs;
				if (!arrlen) arrlen = C->L->H->N;

				for (i=0; i < C->Nl_eff; ++i)
						data[C->Nl_o + i] = (C->L + i)->H->zen;
		}

		if (!GetArg(in, "file", word)) goto eargs;
		printf("Writing horizons to %s\n", word);
		struct h5io* io = h5io_init(word, 0);
		if (NULL == io) goto eio;

		if (!GetOption(in, "dataset", word)) snprintf(word, strlen(in), "horizon");
		h5io_setdataset(io, word);
		if (h5io_isin(io)) {
				Warning("Error: dataset %s exists!\n", io->dataset);
				goto eexist;
		}

		h5io_setdtype(io, "float64");
		if (FetchOptInt(in, "gzip", word, &io->compression)) io->compression = 0;
		if (FetchOptInt(in, "chunksize", word, &io->chunkarr)) io->chunkarr = 10000;
		if (FetchOptInt(in, "cachemb", word, &io->cachemb)) io->cachemb = 64;
		if (FetchOptInt(in, "cacheslots", word, &io->cacheslots)) io->cacheslots = 12421;

		TIC();
		if (h5io_write(io, (void**)data, "float64", arrlen, narr)) goto ewrite;
		char cmmnt[1024];
		strncpy(cmmnt, "SSDP data: write_horizon ", 1024-1);
		strncat(cmmnt, in, 1024-1);
		h5io_comment(io, cmmnt);
		printf("Wrote %s in %g s\n", io->dataset, TOC());

		h5io_free(io);
		free(data);
		free(word);
		return;
ewrite:
eexist:
		h5io_free(io);
eio:
eargs:
einitlocs:
		free(data);
edata:
eC:
		free(word);
eword:
		Warning("write_horizon: failed!\n");
		return;
}
