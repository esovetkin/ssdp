#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "error.h"
#include "vector.h"
#include "freespa/freespa.h"

/* wrapper for the freespa library */

sky_pos sunpos(time_t t, double lat, double lon, double e, double p, double T)
{
	sky_pos s;
	sol_pos S;
	struct tm ut, *tp;
	tp=gmjtime_r(&t, &ut);
	if (tp)
	{	
		
		S=SPA(tp, NULL, 0, lon, lat, e, p, T);
		if (S.E)
		{
			// ERRORFLAG FREESPAERR  "Error: freespa failed" 
			AddErr(FREESPAERR); 
			s.z=0;         
			s.a=0;
		}
		else
		{
			s.z=S.az;         
			s.a=S.aa;
		}
	}   
	else
	{	 
		// ERRORFLAG GMTIMENULL  "Error: gmtime returned NULL" 
		AddErr(GMTIMENULL); 
		s.z=0;         
		s.a=0;
	}
	return s;
}

int suntimes(time_t t, double lat, double lon, double e, double p, double T, time_t *sunrise, time_t *transit, time_t *sunset)
{
	sol_pos S;
	struct tm ut, *tp;
	struct tm Sr, Ss, St;
	int r;
	tp=gmjtime_r(&t, &ut);
	if (tp)
	{	
		r=SunTimes(ut, NULL, 0, lon, lat, e, p, T,&Sr, &St, &Ss);
		if (sunrise)
			(*sunrise)=mkgmjtime(&Sr);
		if (transit)
			(*transit)=mkgmjtime(&St);
		if (sunset)
			(*sunset)=mkgmjtime(&Ss);
		return r; // -1: polar night, 0: normal, 1 midnight sun
	}   
	else
	{	 
		// ERRORFLAG GMTIMENULL  "Error: gmtime returned NULL" 
		AddErr(GMTIMENULL); 
		(*sunrise)=0;         
		(*transit)=0;         
		(*sunset)=0; 
	}
	return 2;// return vlalue >1 indicates an error
}

