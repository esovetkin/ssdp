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
				Warning("Error: no topographical data available\n");
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
						(&(C->S), (time_t) AT(t,j),
						 C->lon, C->lat, C->E,
						 AT(p,j), AT(T,j), AT(GH,j), AT(DH,j));
				tsky+=TOC(); TIC();

#pragma omp parallel private(i) shared(C)
				{
#pragma omp for schedule(runtime)
						for (i=0; i < C->Nl; ++i)
								out.D[j*C->Nl+i]=ssdp_total_poa
										(&(C->S), C->o[i], &(C->M), C->L+i);
				}
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
						(&(C->S), sun, AT(GH,j), AT(DH,j), AT(doy,j));
				tsky += TOC(); TIC();

#pragma omp parallel private (i) shared(C)
				{
#pragma omp for schedule(runtime)
						for (i=0; i < C->Nl; ++i)
								poa.D[j*C->Nl + i] = ssdp_total_poa
										(&(C->S), C->o[i], &(C->M), C->L+i);
				}
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

		if (NULL==(out.D=malloc(N*C->Nl*sizeof(*(out.D))))) goto eout;
		out.N=C->Nl;

		TIC();
		printf("doing static integral simulation of %d locations at %d instances\n",C->Nl, N);
		ssdp_make_perez_cumulative_sky_coordinate
				(&(C->S),t->D, C->lon, C->lat, C->E,
				 p->D, p->N, T->D, T->N, GH->D, DH->D, t->N);
		tsky=TOC();
		printf("Integrated %d skies in %g s (%g s/sky)\n", t->N, tsky, tsky/((double)t->N));

		TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
				for (i=0;i<C->Nl;i++)
						out.D[i]+=ssdp_total_poa
								(&(C->S), C->o[i], &(C->M), C->L+i);
	}
		tpoa=TOC();

		printf("Computed %d POA Irradiances in %g s (%g s/POA)\n", C->Nl, tpoa, tpoa/((double)(C->Nl)));

		printf("Creating array %s\n",nout);
		if(AddArray(nout, out)) goto eoutadd;

		free_default_array(&p, fp);
		free_default_array(&T, fT);
		free(word);
		return;
eoutadd:
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
		int j, N, pco=0, fp=0, fT=0;
		char *word, *nout;
		simulation_config *C;
		array *t, *GH, *DH, *p, *T, out;
		double tsky = 0, tpoa = 0;

		if (NULL==(word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
		if (NULL==(nout=malloc((strlen(in)+1)*sizeof(*nout)))) goto enout;

		if (FetchConfig(in, "C", word, &C)) goto eargs;
		if (check_simconfig(C, 0, 0)) goto eargs;
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
		if (N>C->Nl)
				Warning("Warning: time array contains more points than there are waypoints\n");

		if (NULL==(out.D=malloc(N*sizeof(*out.D)))) goto eout;
		out.N=N;
		printf("doing route simulation along %d locations at %d instances\n",C->Nl, N);
		for (j=0; j < N; ++j) {
				TIC();
				ssdp_make_perez_all_weather_sky_coordinate
						(&(C->S), (time_t) AT(t,j),
						 C->lon, C->lat, C->E,
						 AT(p,j), AT(T,j), AT(GH,j), AT(DH,j));
				tsky += TOC(); TIC();
				out.D[j]=ssdp_total_poa(&(C->S), C->o[j], &(C->M), C->L+j);
				ProgressBar((100*(j+1))/N, &pco, ProgressLen, ProgressTics);
				tpoa+=TOC();
		}

		printf("Computed %d skies in %g s (%g s/sky)\n",
			   N, tsky, tsky/((double)N));
		printf("Computed %d POA Irradiances in %g s (%g s/POA)\n",
			   N, tpoa, tpoa/((double)N));

		printf("Creating array %s\n", nout);
		if(AddArray(nout, out)) goto eoutadd;

		free_default_array(&p, fp);
		free_default_array(&T, fT);
		free(word);
		return;
eoutadd:
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
				// only update sun position, we do not need to compute the diffuse sky, it is uniform and boring
				ssdp_make_skysunonly_coordinate
						(&(C->S), (time_t) AT(t,j),
						 C->lon, C->lat, C->E,
						 AT(p,j), AT(T,j), AT(GH,j), AT(DH,j));
				tsky+=TOC(); TIC();
#pragma omp parallel private(i) shared(C)
		{
#pragma omp for schedule(runtime)
			for (i=0;i < C->Nl; ++i) {
					out.D[j*C->Nl+i]=ssdp_direct_poa
							(&(C->S), C->o[i], &(C->M), C->L+i);
					out.D[j*C->Nl+i]+=C->L[i].difftrans*DH->D[j]; // diffuse part is directly proportional to diffuse horizontal
			}
		}
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
		sky_pos s;
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
				&(C->S), (time_t) t, C->lon, C->lat, C->E, p, T, GH, DH);
		WriteDome4D(file, &(C->S), C->L+j);

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
