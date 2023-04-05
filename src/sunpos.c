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
		
		S=SPA(tp, NULL, 0, lon, lat, e);
		S=ApSolposBennet(S,NULL, e, p, T);
		if (S.E)
		{
			// ERRORFLAG FREESPAERR  "Error: freespa failed" 
			AddErr(FREESPAERR); 
			s.z=0;         
			s.a=0;
		}
		else
		{
			s.z=S.z;         
			s.a=S.a;
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
	solar_day D;
	struct tm ut, *tp;
	int r=0;
	tp=gmjtime_r(&t, &ut);
	if (tp)
	{	
		// only compute sunrise and sunset
		SDMASK=(_FREESPA_SUNRISE|_FREESPA_SUNSET);
		D=SolarDay(tp, NULL, 0, lon, lat, e, NULL, p, T, &ApSolposBennet);
						
		if (D.status[3]==0)
			(*sunrise)=D.t[3];
		else
		{
			r|=1;
			(*sunrise)=0;       
		}
		if (D.status[1]==0)
			(*transit)=D.t[1];
		else
		{
			r|=2;   
			(*transit)=0;      
		}
		if (D.status[4]==0)
			(*sunset)=D.t[4];
		else
		{
			r|=4;       
			(*sunset)=0; 
		}
		return r;
	}   
	else
	{	 
		// ERRORFLAG GMTIMENULL  "Error: gmtime returned NULL" 
		AddErr(GMTIMENULL); 
		(*sunrise)=0;         
		(*transit)=0;         
		(*sunset)=0; 
	}
	return 8;// return vlalue >0 indicates an error
}

