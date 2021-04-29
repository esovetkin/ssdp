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


// PARSEFLAG sample_topo SampleTopography "C=<config-variable>  x=<x-array-variable> y=<y-array-variable> z=<z-output-array> azimuth=<azimuth-output-array> zenith=<zenith-output-array>"
void SampleTopography(char *in)
{
	int i;
	char *word;
	simulation_config *C;
	array *x, *y, z, azi, zen;
	sky_pos sn;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}	
	if (!C->topo_init) // TODO: in this case just omit the horizon and compute only one location?
	{	
		Warning("Simulation config has no topology initialized\n");
		free(word);
		return;
	}	
	if (FetchArray(in, "x", word, &x))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "y", word, &y))
	{
		free(word);
		return;
	}
	if (x->N!=y->N)
	{
		Warning("Length of x- and y-arrays do not match\n");
		free(word);
		return;
	}
	z.D=malloc(x->N*sizeof(double));
	z.N=x->N;
	azi.D=malloc(x->N*sizeof(double));
	azi.N=x->N;
	zen.D=malloc(x->N*sizeof(double));
	zen.N=x->N;
	for (i=0;i<x->N;i++)
	{
		z.D[i]=ssdp_sample_topology(x->D[i], y->D[i], &(C->T),&sn);
		azi.D[i]=sn.a;
		zen.D[i]=sn.z;
	}
	
	if (GetOption(in, "z", word))
	{
		printf("Creating array %s\n",word);
		if(AddArray(word, z))
		{
			Warning("Failed to create array %s\n",word);
			free(z.D);	
		}
		else
			word=malloc((strlen(in)+1)*sizeof(char));
	}
	else
		free(z.D);
		
	if (GetOption(in, "azimuth", word))
	{
		printf("Creating array %s\n",word);
		if(AddArray(word, azi))
		{
			Warning("Failed to create array %s\n",word);
			free(azi.D);	
		}
		else
			word=malloc((strlen(in)+1)*sizeof(char));
	}
	else
		free(azi.D);
		
	if (GetOption(in, "zenith", word))
	{
		printf("Creating array %s\n",word);
		if(AddArray(word, zen))
		{
			Warning("Failed to create array %s\n",word);
			free(zen.D);	
		}
		else
			word=NULL;
	}
	else
		free(zen.D);
	if (word)
		free(word);
}
// PARSEFLAG offset_topo OffsetTopography "C=<config-variable>  o=<offset-value> x=<x-array-variable> y=<y-array-variable> xoff=<offset-x-output-array> yoff=<offset-x-output-array> zoff=<offset-z-output-array>"
void OffsetTopography(char *in)
{
	int i;
	char *word;
	simulation_config *C;
	array *x, *y, zoff, yoff, xoff;
	double o;
	sky_pos sn;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}	
	if (!C->topo_init)
	{	
		Warning("Simulation config has no topology initialized\n");
		free(word);
		return;
	}	
	if (FetchFloat(in, "o", word, &o))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "x", word, &x))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "y", word, &y))
	{
		free(word);
		return;
	}
	if (x->N!=y->N)
	{
		Warning("Length of x- and y-arrays do not match\n");
		free(word);
		return;
	}
	zoff.D=malloc(x->N*sizeof(double));
	zoff.N=x->N;
	yoff.D=malloc(x->N*sizeof(double));
	yoff.N=x->N;
	xoff.D=malloc(x->N*sizeof(double));
	xoff.N=x->N;
	for (i=0;i<x->N;i++)
	{
		zoff.D[i]=ssdp_sample_topology(x->D[i], y->D[i], &(C->T),&sn)+o*cos(sn.z);
		xoff.D[i]=x->D[i]+o*sin(sn.z)*cos(sn.a);
		yoff.D[i]=y->D[i]+o*sin(sn.z)*sin(sn.a);
	}
	
	
	if (!GetArg(in, "zoff", word))
	{
		free(word);
		return;
	}
	printf("Creating array %s\n",word);
	if(AddArray(word, zoff))
	{
		Warning("Failed to create array %s\n",word);
		free(zoff.D);	
	}
	else
		word=malloc((strlen(in)+1)*sizeof(char));
		
	if (!GetArg(in, "xoff", word))
	{
		free(word);
		return;
	}	
	printf("Creating array %s\n",word);
	if(AddArray(word, xoff))
	{
		Warning("Failed to create array %s\n",word);
		free(xoff.D);	
	}
	else
		word=malloc((strlen(in)+1)*sizeof(char));
		
	if (!GetArg(in, "yoff", word))
	{
		free(word);
		return;
	}
	printf("Creating array %s\n",word);
	if(AddArray(word, yoff))
	{
		Warning("Failed to create array %s\n",word);
		free(yoff.D);	
	}
	else
		word=NULL;
	
	if (word)
		free(word);
}
// PARSEFLAG rotate_POAto_surface RotatePOA "poa_a=<azimuth-value>  poa_z=<zenith-value> surf_a=<azimuth-array-variable> surf_z=<zenith-array-variable> out_a=<azimuth-output-array> out_z=<zenith-output-array>"
void RotatePOA(char *in)
{
	int i;
	char *word;
	array *rot_a, *rot_z, azi, zen;
	sky_pos poa, poa0;
	sky_pos sn;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchFloat(in, "poa_a", word, &poa0.a))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "poa_z", word, &poa0.z))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "surf_a", word, &rot_a))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "surf_z", word, &rot_z))
	{
		free(word);
		return;
	}
	if (rot_a->N!=rot_z->N)
	{
		Warning("Length of azimuth and zenoth of the surface normal do not match\n");
		free(word);
		return;
	}
	azi.D=malloc(rot_a->N*sizeof(double));
	azi.N=rot_a->N;
	zen.D=malloc(rot_a->N*sizeof(double));
	zen.N=rot_a->N;
	for (i=0;i<rot_a->N;i++)
	{
		sn.a=rot_a->D[i];
		sn.z=rot_z->D[i];
		ssdp_poa_to_surface_normal(poa0, sn, &poa); // orient module w.r.t ground orientation
		azi.D[i]=poa.a;
		zen.D[i]=poa.z;
	}
	
	
	if (!GetArg(in, "out_a", word))
	{
		free(word);
		return;
	}
	printf("Creating array %s\n",word);
	if(AddArray(word, azi))
	{
		Warning("Failed to create array %s\n",word);
		free(azi.D);	
	}
	else
		word=malloc((strlen(in)+1)*sizeof(char));
		
	if (!GetArg(in, "out_z", word))
	{
		free(word);
		return;
	}	
	printf("Creating array %s\n",word);
	if(AddArray(word, zen))
	{
		Warning("Failed to create array %s\n",word);
		free(zen.D);	
	}
	else
		word=NULL;
	if (word)
		free(word);
}
