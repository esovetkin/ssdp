#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "libssdp.h"
#include "util.h"

#define MAXSTRLEN 1028

topology LoadTopo(char *fn)
{
	char c, *line;
	int k, Na=50, N;
	FILE *f;
	double *x, *y, *z;
	topology T;	
	if ((f=fopen(fn,"rb"))==NULL)
		Fatal("Cannot open %s for reading\n", fn); // should return a null topography so we can recover
		
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

void WriteTopo(char *fn, topology *T)
{
	FILE *f;
	int i;
	if ((f=fopen(fn,"w"))==NULL)
		Fatal("Cannot open %s for writing\n", fn); // should not kill myself over this
	fprintf(f,"# topology data\n");
	fprintf(f,"# x y z)\n");
	for (i=0;i<T->N;i++)
		fprintf(f, "%e %e %e\n", T->x[i], T->y[i], T->z[i]);
	fclose(f);
}	
void WriteTriangles(char *fn, topology *T)
{
	FILE *f;
	int i;
	if ((f=fopen(fn,"w"))==NULL)
		Fatal("Cannot open %s for writing\n", fn);
	fprintf(f,"# topology triangle data\n");
	fprintf(f,"# x y z index)\n");
	for (i=0;i<T->Nt;i++)
	{
		fprintf(f, "%e %e %e %d\n", T->x[T->T[i].i], T->y[T->T[i].i], T->z[T->T[i].i], i);
		fprintf(f, "%e %e %e %d\n", T->x[T->T[i].j], T->y[T->T[i].j], T->z[T->T[i].j], i);
		fprintf(f, "%e %e %e %d\n", T->x[T->T[i].k], T->y[T->T[i].k], T->z[T->T[i].k], i);
		fprintf(f, "%e %e %e %d\n\n", T->x[T->T[i].i], T->y[T->T[i].i], T->z[T->T[i].i], i);
	}
	fclose(f);
}
/* sun solid angle:
 * fe '(appi(64)*(1391400.0/2)^2)/(149.597870700e6^2)'
 * diameter sun:   1391400    km
 * distance sun: 149597870.7  km
 * solid angle : ~ A/d^2 = (pi * 1391400 ^2 / 4) / 149597870.7^2 = 6.7942739713694071218e-05
 * Note that the actual distance sun-earth varies about 3% so the result is not to be taken as accurate
 */
#define SUNSR 6.794e-5	
void WriteDome4D(char *fn, sky_grid *sky, location *l)
{
	FILE *f;
	int i;
	if ((f=fopen(fn,"w"))==NULL)
	{
		fprintf(stderr,"Cannot open %s for writing\n", fn);
		exit(1);
	}
	fprintf(f,"# 4D Dome plot\n");
	fprintf(f,"# x[-]\ty[-]\tz\tI [W/sr]\n");
	if (!ssdp_below_horizon(l, sky->sp))
		fprintf(f,"#sun: %e %e %e %e\n",sin(sky->sp.a)*sin(sky->sp.z), cos(sky->sp.a)*sin(sky->sp.z), cos(sky->sp.z), sky->N*sky->sI/SUNSR); // W/sr
	else
		fprintf(f,"#sun: %e %e %e 0\n",sin(sky->sp.a)*sin(sky->sp.z), cos(sky->sp.a)*sin(sky->sp.z), cos(sky->sp.z)); // W/sr
	
	for (i=0;i<sky->N;i++)
		if (!ssdp_below_horizon(l, sky->P[i].p))
			fprintf(f,"%e %e %e %e\n",sin(sky->P[i].p.a)*sin(sky->P[i].p.z), cos(sky->P[i].p.a)*sin(sky->P[i].p.z), cos(sky->P[i].p.z), sky->N*sky->P[i].I/2/M_PI); // W/sr
		else
			fprintf(f,"%e %e %e 0\n",sin(sky->P[i].p.a)*sin(sky->P[i].p.z), cos(sky->P[i].p.a)*sin(sky->P[i].p.z), cos(sky->P[i].p.z)); // W/sr
	fclose(f);	
}

double **ReadArrays(char *fn, int Narr, int *N)
{
	char c, *line, *p, *q;
	int i, k, Na=50, n;
	FILE *f;
	double **data;
	if ((f=fopen(fn,"rb"))==NULL)
		Fatal("Cannot open %s for reading\n", fn);
		
	line=malloc(MAXSTRLEN*sizeof(char));
    fgets(line, MAXSTRLEN-1, f);
	n=0;
	data=malloc(Narr*sizeof(double *));
	for (i=0;i<Narr;i++)
		data[i]=malloc(Na*sizeof(double));
		
	while(feof(f)==0)
	{
		
    	k=sscanf(line, " %c", &c);
		if((k==1)&&(c!='#'))
		{
			p=line;
			while ((isblank(*p))&&*p)
				p++;
			for (i=0;i<Narr;i++)
			{
				q=p;
				while ((!isblank(*p))&&*p)
					p++;
				c=*p;
				*p='\0';
				data[i][n]=atof(q);
				*p=c;
				while ((isblank(*p))&&*p)
					p++;
				if (!*p)
					break;
			}
			if (i>=Narr-1)
				n++;
					
			if (Na-1==n)
			{
				Na+=50;
				for (i=0;i<Narr;i++)
					data[i]=realloc(data[i], Na*sizeof(double));
			}
			
		}
    	fgets(line, MAXSTRLEN-1, f);
	}
	free(line);
	fclose(f);
	(*N)=n;
	return data;	
}


void WriteArrays(char *fn, double **data, int Narr, int N)
{
	FILE *f;
	int i,j;
	if ((f=fopen(fn,"w"))==NULL)
		Fatal("Cannot open %s for writing\n", fn);
	for (i=0;i<N;i++)
	{
		for (j=0;j<Narr-1;j++)
			fprintf(f, "%12e\t",data[j][i]);
		fprintf(f, "%12e\n",data[Narr-1][i]);			
	}	
	fclose(f);	
}




