#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libssdp.h"

#define MAXSTRLEN 1028

topology LoadTopo(char *fn)
{
	char c, *line;
	int k, Na=50, N;
	FILE *f;
	double *x, *y, *z;
	topology T;	
	if ((f=fopen(fn,"rb"))==NULL)
	{
		fprintf(stderr,"Cannot open %s for reading\n", fn);
		exit(1);
	}
	line=malloc(MAXSTRLEN*sizeof(char));
    fgets(line, MAXSTRLEN-1, f);
	N=0;
	x=malloc(Na*sizeof(double));
	y=malloc(Na*sizeof(double));
	z=malloc(Na*sizeof(double));
	while(feof(f)==0)
	{
		
    	k=sscanf(line, " %c", &c);
		if((k==1)&&(c!='#'))
		{
			k=sscanf(line, " %le %le %le", x+N, y+N, z+N);
			if(k==3)
			{
				N++;
				if (Na-1==N)
				{
					Na+=50;
					x=realloc(x, Na*sizeof(double));
					y=realloc(y, Na*sizeof(double));
					z=realloc(z, Na*sizeof(double));	
				}
			}
		}
    	fgets(line, MAXSTRLEN-1, f);
	}
	free(line);
	fclose(f);
	T=ssdp_make_topology(x, y, z, N);
	free(x);
	free(y);
	free(z);
	return T;
}

void WriteTopo(char *fn, topology T)
{
	FILE *f;
	int i;
	if ((f=fopen(fn,"w"))==NULL)
	{
		fprintf(stderr,"Cannot open %s for writing\n", fn);
		exit(1);
	}
	fprintf(f,"# topology data\n");
	fprintf(f,"# x y z)\n");
	for (i=0;i<T.N;i++)
		fprintf(f, "%e %e %e\n", T.x[i], T.y[i], T.z[i]);
	fclose(f);
}	
void WriteTriangles(char *fn, topology T)
{
	FILE *f;
	int i;
	if ((f=fopen(fn,"w"))==NULL)
	{
		fprintf(stderr,"Cannot open %s for writing\n", fn);
		exit(1);
	}
	fprintf(f,"# topology triangle data\n");
	fprintf(f,"# x y z)\n");
	for (i=0;i<T.Nt;i++)
	{
		fprintf(f, "%e %e %d\n", T.x[T.T[i].i], T.y[T.T[i].i], i);
		fprintf(f, "%e %e %d\n", T.x[T.T[i].j], T.y[T.T[i].j], i);
		fprintf(f, "%e %e %d\n", T.x[T.T[i].k], T.y[T.T[i].k], i);
		fprintf(f, "%e %e %d\n\n", T.x[T.T[i].i], T.y[T.T[i].i], i);
	}
	fclose(f);
}	
	
void WriteDome3D(char *fn, sky_grid sky, int sun, int mask)
{
	FILE *f;
	int i, si=-1;
	if ((f=fopen(fn,"w"))==NULL)
	{
		fprintf(stderr,"Cannot open %s for writing\n", fn);
		exit(1);
	}
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
	{
		fprintf(stderr,"Cannot open %s for writing\n", fn);
		exit(1);
	}
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

void RasterPOA(char *fn, sky_grid sky, topology T, double albedo, double dz, double a, double tilt, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	FILE *f;
	int i,j;
	double x, y;
	VERB verbstate;
	verbstate=ssdp_verbosity;
	ssdp_verbosity=QUIET;
	if ((f=fopen(fn,"w"))==NULL)
	{
		fprintf(stderr,"Cannot open %s for writing\n", fn);
		exit(1);
	}
	fprintf(f,"# raster plot of topology\n");
	for (i=0;i<Nx;i++)
	{
		x=x1+(x2-x1)*((double)i+0.5)/((double)Nx);
		for (j=0;j<Ny;j++)
		{
			y=y1+(y2-y1)*((double)j+0.5)/((double)Ny);
			ssdp_mask_horizon_z_to_ground(&sky,T,x,y,dz, NULL);
			fprintf(f,"%e %e %e\n", x, y, ssdp_total_poa(sky,albedo,tilt, a,1));
			ssdp_unmask_horizon(&sky);
		}
	}	
	ssdp_verbosity=verbstate;
	fclose(f);
}

void RasterTopology(char *fn, topology T, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	FILE *f;
	int i,j;
	double x, y;
	if ((f=fopen(fn,"w"))==NULL)
	{
		fprintf(stderr,"Cannot open %s for writing\n", fn);
		exit(1);
	}
	fprintf(f,"# raster plot of topology\n");
	for (i=0;i<Nx;i++)
	{
		x=x1+(x2-x1)*((double)i+0.5)/((double)Nx);
		for (j=0;j<Ny;j++)
		{
			y=y1+(y2-y1)*((double)j+0.5)/((double)Ny);
			fprintf(f,"%e %e %e\n", x, y, ssdp_sample_topology(x,y,T, NULL));
		}
	}	
	fclose(f);
}




