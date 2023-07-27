#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
/* local includes */
#include "libssdp.h"
#include "lio.h"
#include "util.h"
#include "variables.h"
#include "parser.h"
#include "parserutil.h"


/*
BEGIN_DESCRIPTION
SECTION Topography
PARSEFLAG sample_topo SampleTopography "C=<in-config>  x=<in-array> y=<in-array> type=<topology/topogrid> z=<out-array> azimuth=<out-array> zenith=<out-array>"
DESCRIPTION Samples a topography to obtain the local height and surface normal
ARGUMENT C config-variable
ARGUMENT x x coordinate
ARGUMENT y y coordinate
ARGUMENT type Optional argument to select the unstructured topology mesh or the structured topogrid. Valid values are \"topology\" or \"topogrid\" (default topology unless only a topogrid is defined)
OUTPUT z z coordinate
OUTPUT azimuth surface normal azimuth
OUTPUT zenith surface normal zenith
END_DESCRIPTION
*/
void SampleTopography(char *in)
{
	int i;
	char *word;
	simulation_config *C;
	array *x, *y, z, azi, zen;
	char type='t';
	sky_pos sn;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
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
	if (GetOption(in, "type", word))
	{
		if (strncmp(word,"topology",10)==0)
		{
			if (C->topo_init==0)
			{ 
				Warning("No topology available\n");
				free(word);
				return;
			}
		}
		if (strncmp(word,"topogrid",10)==0)
		{
			if (C->grid_init==0)
			{ 
				Warning("No topogrid available\n");
				free(word);
				return;
			}
			type='g';
		}
	}
	if (C->topo_init==0)
	{
		type='g';
		if (C->grid_init==0)
		{ 
			Warning("No topology or topogrid available\n");
			free(word);
			return;
		}
	}
	
	z.D=malloc(x->N*sizeof(double));
	z.N=x->N;
	azi.D=malloc(x->N*sizeof(double));
	azi.N=x->N;
	zen.D=malloc(x->N*sizeof(double));
	zen.N=x->N;
	for (i=0;i<x->N;i++)
	{
		if (type=='t')
			z.D[i]=ssdp_sample_topology(x->D[i], y->D[i], &(C->T),&sn);
		else
			z.D[i]=ssdp_sample_topogrid(x->D[i], y->D[i], &(C->Tx),&sn);
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
/*
BEGIN_DESCRIPTION
SECTION Topography
PARSEFLAG offset_topo OffsetTopography "C=<in-config> o=<in-array> x=<in-array> y=<in-array> type=<topology/topogrid> xoff=<out-array> yoff=<out-array> zoff=<out-array>"
DESCRIPTION Computes a topography offset in the direction of the surface normal
ARGUMENT C config-variable
ARGUMENT o offset value
ARGUMENT x x coordinate
ARGUMENT y y coordinate
ARGUMENT type Optional argument to select the unstructured topology mesh or the structured topogrid. Valid values are \"topology\" or \"topogrid\" (default topology unless only a topogrid is defined)
OUTPUT xoff x coordinate
OUTPUT yoff y coordinate
OUTPUT zoff z coordinate
END_DESCRIPTION
*/
void OffsetTopography(char *in)
{
	int i;
	char *word;
	simulation_config *C;
	array *o, *x, *y, zoff, yoff, xoff;
	char type='t';
	sky_pos sn;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}	
	if ((!C->topo_init)&&(!C->grid_init))
	{	
		Warning("Simulation config has no topology or topogrid initialized\n");
		free(word);
		return;
	}	
	if (FetchArray(in, "o", word, &o))
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
	if (x->N%o->N)
		Warning("Length x-array (%d) not divisable by the length offset array (%d)\n", x->N, o->N);
		
	if (GetOption(in, "type", word))
	{
		if (strncmp(word,"topology",10)==0)
		{
			if (C->topo_init==0)
			{ 
				Warning("No topology available\n");
				free(word);
				return;
			}
		}
		if (strncmp(word,"topogrid",10)==0)
		{
			if (C->grid_init==0)
			{ 
				Warning("No topogrid available\n");
				free(word);
				return;
			}
			type='g';
		}
	}
	if (C->topo_init==0)
	{
		type='g';
		if (C->grid_init==0)
		{ 
			Warning("No topology or topogrid available\n");
			free(word);
			return;
		}
	}
		
	zoff.D=malloc(x->N*sizeof(double));
	zoff.N=x->N;
	yoff.D=malloc(x->N*sizeof(double));
	yoff.N=x->N;
	xoff.D=malloc(x->N*sizeof(double));
	xoff.N=x->N;
	printf("computing tolology offset\n");
	for (i=0;i<x->N;i++)
	{
		if (type=='t')
			zoff.D[i]=ssdp_sample_topology(x->D[i], y->D[i], &(C->T),&sn)+o->D[i%o->N]*cos(sn.z);
		else
			zoff.D[i]=ssdp_sample_topogrid(x->D[i], y->D[i], &(C->Tx),&sn)+o->D[i%o->N]*cos(sn.z);
		
		xoff.D[i]=x->D[i]+o->D[i%o->N]*sin(sn.z)*sin(sn.a);
		yoff.D[i]=y->D[i]+o->D[i%o->N]*sin(sn.z)*cos(sn.a);
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


static int checkdims(int a, int b)
{
        if (a != b && 1 != a && 1 != b)
                return -1;

        return a >= b ? a : b;
}


/*
BEGIN_DESCRIPTION
SECTION Topography
PARSEFLAG rotate_POA_to_surface RotatePOA "poa_a=<in-array> poa_z=<in-array> surf_a=<in-array> surf_z=<in-array> out_a=<out-array> out_z=<out-array>"
DESCRIPTION Rotate the plane of array (tilted surface) along the surface normal, e.g. a surface on a vehicle changes its tilt an orientation depending on the sufrace nomal. All input array lengths must be equal or it must be a unit length array.
ARGUMENT poa_a azimuth angle of the tilted surface (for a horizonal surface)
ARGUMENT poa_z zenith angle of the tilted surface (for a horizonal surface)
ARGUMENT surf_a azimuth angle of the surface normal
ARGUMENT surf_z zenith angle of the surface normal
OUTPUT out_a azimuth angle of the tilted surface rotated along the surface normal
OUTPUT out_z zenith angle of the tilted surface rotated along the surface normal
END_DESCRIPTION
*/
void RotatePOA(char *in)
{
        int i, Nsurf, Npoa, N;
        char *word;
        array *rot_a, *rot_z, azi, zen;
        array *azi_in, *zen_in;
        sky_pos poa, poa0;
        sky_pos sn;
        word=malloc((strlen(in)+1)*sizeof(char));

        if (FetchArray(in, "poa_a", word, &azi_in)) goto earg;
        if (FetchArray(in, "poa_z", word, &zen_in)) goto earg;
        if (FetchArray(in, "surf_a", word, &rot_a)) goto earg;
        if (FetchArray(in, "surf_z", word, &rot_z)) goto earg;

        if ((Nsurf = checkdims(rot_a->N, rot_z->N)) < 0) {
                Warning("len(surf_a) != len(surf_z) and neither is unit\n");
                goto earg;
        }
        if ((Npoa = checkdims(azi_in->N, zen_in->N)) < 0) {
                Warning("len(poa_a) != len(poa_z) and neither is unit\n");
                goto earg;
        }
        if ((N = checkdims(Npoa, Nsurf)) < 0) {
                Warning("len(poa) != len(surf) and neither is unit\n");
                goto earg;
        }

        if (NULL == (azi.D=malloc(N*sizeof(*azi.D)))) goto eoutazi;
        azi.N = N;
        if (NULL == (zen.D=malloc(N*sizeof(*zen.D)))) goto eoutzen;
        zen.N = N;

        printf("rotating the POA along the surface normal\n");
        for (i=0; i < N; ++i) {
                sn.a=rot_a->D[i % rot_a->N];
                sn.z=rot_z->D[i % rot_z->N];
                poa0.a=azi_in->D[i % azi_in->N];
                poa0.z=zen_in->D[i % zen_in->N];
                ssdp_poa_to_surface_normal(poa0, sn, &poa); // orient module w.r.t ground orientation
                azi.D[i]=poa.a;
                zen.D[i]=poa.z;
        }

        if (!GetArg(in, "out_a", word)) goto eoutarg;
        printf("Creating array %s\n",word);
        if(AddArray(word, azi)) {
                Warning("Failed to create array %s\n",word);
                goto eoutarg;
        }

        if (NULL == (word=malloc((strlen(in)+1)*sizeof(char)))) goto eoutarg;
        if (!GetArg(in, "out_z", word)) goto eoutarg;

        printf("Creating array %s\n",word);
        if(AddArray(word, zen)) {
                Warning("Failed to create array %s\n",word);
                goto eoutarg;
        }

        return;
eoutarg:
        free(zen.D);
eoutzen:
        free(azi.D);
eoutazi:
earg:
        free(word);
}

/*
BEGIN_DESCRIPTION
SECTION Topography
PARSEFLAG bearing Bearing "x=<in-array> y=<in-array> azimuth=<out-array>" 
DESCRIPTION Compute the bearing from the waypoints of a route.
ARGUMENT x x coordinates of the waypoints
ARGUMENT y y coordinates of the waypoints
OUTPUT azimuth bearing azimuth angle
END_DESCRIPTION
*/
void Bearing(char *in)
{
	char *word;
	array *x, *y, azi;
	int i;
	word=malloc((strlen(in)+1)*sizeof(char));
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
		Warning("Length of x and y arrays do not match\n");
		free(word);
		return;
	}
	azi.D=malloc(x->N*sizeof(double));
	azi.N=x->N;
		
	printf("compute bearing\n");
	for(i=0;i<x->N-2;i++)
		azi.D[i+1]=atan2(x->D[i+2]-x->D[i], y->D[i+2]-y->D[i]);
	azi.D[0]=atan2(x->D[1]-x->D[0], y->D[1]-y->D[0]);
	azi.D[azi.N-1]=atan2(x->D[azi.N-1]-x->D[azi.N-2], y->D[azi.N-1]-y->D[azi.N-2]);
	if (!GetArg(in, "azimuth", word))
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
	return;	
}

/*
BEGIN_DESCRIPTION
SECTION Topography
PARSEFLAG export_topo ExportTopo "C=<in-config> file=<file-str>"
DESCRIPTION Export a topography.
ARGUMENT C Configuration variable
OUTPUT file File with x,y,z columns
END_DESCRIPTION
*/
void ExportTopo(char *in)
{
	char *word;
	simulation_config *C;
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
	if (!GetArg(in, "file", word))
	{
		free(word);
		return;
	}
	printf("writing tolography to file %s\n", word);
	WriteTopo(word, &(C->T));
	return;	
}

/*
BEGIN_DESCRIPTION
SECTION Topography
PARSEFLAG export_triangles ExportTriangles "C=<in-config> file=<file-str>"
DESCRIPTION Export a triangulation.
ARGUMENT C Configuration variable
OUTPUT file File with x,y,z columns and one block per triangle with the triangle vertices in right winding order.
END_DESCRIPTION
*/
void ExportTriangles(char *in)
{
	char *word;
	simulation_config *C;
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
	if (!GetArg(in, "file", word))
	{
		free(word);
		return;
	}
	printf("writing triangulation to file %s\n", word);
	WriteTriangles(word, &(C->T));
	free(word);
	return;	
}
