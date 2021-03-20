#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libssdp.h"
#include "util.h"


void WriteDome3D(char *fn, sky_grid sky, int sun, int mask)
{
	FILE *f;
	int i, si=-1;
	if ((f=fopen(fn,"w"))==NULL)
		Fatal("Cannot open %s for writing\n", fn);
	fprintf(f,"# Polar Plot azimuth and zenith\n");
	fprintf(f,"# x=pi * sin(z) cos(a), y=pi * sin(z) sin(a)\n");
	fprintf(f,"# i.e. distance from origin is the zenith in radians\n");
	fprintf(f,"# 3rd Column (in W/sr)\n");
	if (mask)
		fprintf(f,"# Horizon included\n");
	else
		fprintf(f,"# Horizon ignored\n");
	if (sun)
	{
		si=ssdp_find_skypatch(sky,sky.sp);
		fprintf(f,"# Direct light included\n");
	}
	else
		fprintf(f,"# Diffuse only\n");
	for (i=0;i<sky.N;i++)
	{
		if ((sky.P[i].mask)&&(mask))
			fprintf(f,"%e %e 0\n",M_PI*cos(sky.P[i].p.a)*sin(sky.P[i].p.z)/2, M_PI*sin(sky.P[i].p.a)*sin(sky.P[i].p.z)/2); // W/sr
		else
		{
			if (i==si)
				fprintf(f,"%e %e %e\n",M_PI*cos(sky.P[i].p.a)*sin(sky.P[i].p.z)/2, M_PI*sin(sky.P[i].p.a)*sin(sky.P[i].p.z)/2, sky.N*(sky.P[i].I+sky.sI)/2/M_PI); // W/sr
			else
				fprintf(f,"%e %e %e\n",M_PI*cos(sky.P[i].p.a)*sin(sky.P[i].p.z)/2, M_PI*sin(sky.P[i].p.a)*sin(sky.P[i].p.z)/2, sky.N*sky.P[i].I/2/M_PI); // W/sr
		}
	}	
	fclose(f);
}
void WriteDome4D(char *fn, sky_grid sky, int sun, int mask)
{
	FILE *f;
	int i, si=-1;
	if ((f=fopen(fn,"w"))==NULL)
		Fatal("Cannot open %s for writing\n", fn);
	fprintf(f,"# 4D Dome plot\n");
	if (mask)
		fprintf(f,"# Horizon included\n");
	else
		fprintf(f,"# Horizon ignored\n");
	if (sun)
	{
		si=ssdp_find_skypatch(sky,sky.sp);
		fprintf(f,"# Direct light included\n");
	}
	else
		fprintf(f,"# Diffuse only\n");
	fprintf(f,"# normalized distances\n");
	fprintf(f,"# x[-]\ty[-]\tz\tI [W/sr]\n");
	for (i=0;i<sky.N;i++)
	{
		if ((sky.P[i].mask)&&(mask))
			fprintf(f,"%e %e %e 0\n",cos(sky.P[i].p.a)*sin(sky.P[i].p.z), sin(sky.P[i].p.a)*sin(sky.P[i].p.z), cos(sky.P[i].p.z)); // W/sr
		else
		{
			if (i==si)
				fprintf(f,"%e %e %e %e\n",cos(sky.P[i].p.a)*sin(sky.P[i].p.z), sin(sky.P[i].p.a)*sin(sky.P[i].p.z), cos(sky.P[i].p.z), sky.N*(sky.P[i].I+sky.sI)/2/M_PI); // W/sr
			else
				fprintf(f,"%e %e %e %e\n",cos(sky.P[i].p.a)*sin(sky.P[i].p.z), sin(sky.P[i].p.a)*sin(sky.P[i].p.z), cos(sky.P[i].p.z), sky.N*sky.P[i].I/2/M_PI); // W/sr
		}
	}	
	fclose(f);
	
}

