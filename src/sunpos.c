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
	tp=gmtime_r(&t, &ut);
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

