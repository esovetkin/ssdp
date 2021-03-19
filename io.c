#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "sky_dome.h"
#include "sky_model.h"
#include "project.h"
#include "ground.h"
#include "util.h"


#define MAXSTRLEN 1028
// very basic topology reading routine
// if we will use a commandline version within PV-GRIP
// we better think of a binary file format to speed thins up and keep memory usage low
topology LoadTopo(char *fn)
{
	topology res;
	char c, *line;
	int k, Na=50;
	FILE *f;
	int ddef=0;
	if ((f=fopen(fn,"rb"))==NULL)
		Fatal("Cannot open %s for reading\n", fn);
		
	line=malloc(MAXSTRLEN*sizeof(char));
    fgets(line, MAXSTRLEN-1, f);
	res.N=0;
	res.points=malloc(Na*sizeof(vec));
	while(feof(f)==0)
	{
		
    	k=sscanf(line, " %c", &c);
		if((k==1)&&(c!='#'))
		{
			k=sscanf(line, " %le %le %le", &(res.points[res.N].x), &(res.points[res.N].y), &(res.points[res.N].z));
			if(k==3)
			{
				res.N++;
				if (Na-1==res.N)
				{
					Na+=50;
					res.points=realloc(res.points, Na*sizeof(vec));			
				}
			}
		}
		else if (res.N==0)
		{
			k=sscanf(line, "# d=%le", &(res.d));
			if (k==1)
				ddef=1;
				
		}
    	fgets(line, MAXSTRLEN-1, f);
	}
	free(line);
	fclose(f);
	if (!ddef)
	{
		Warning("Warning: no point diameter d specified, setting to 1\n");
		res.d=1;
	}
	return res;
}
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
		si=FindPatch(sky,sky.sp);
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
		si=FindPatch(sky,sky.sp);
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

