/*
    Simple Sky-Dome Projector Library
    Copyright (C) 2021  B. E. Pieters, 
    IEK-5 Photovoltaik, Forschunszentrum Juelich

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "vector.h"
#include "sky_dome.h"
#include "delaunay.h"
#include "ground.h"
#include "print.h"
#include "print.h"
#include "error.h"


topology MakeTopology(double *x, double *y, double *z, int N)
{
	topology T;
	topology T0={NULL,NULL,NULL,0,NULL,0,NULL};
	int i;
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--MakeTopology\t\t\t");
	Print(VVERBOSE, "\ncreating topology: %d points\n", N);
	// ERRORFLAG MALLOCFAILTOPOLOGY  "Error memory allocation failed in creating a topology"
	T.x=malloc(N*sizeof(double));
	T.y=malloc(N*sizeof(double));
	T.z=malloc(N*sizeof(double));
	if ((T.x==NULL)||(T.y==NULL)||(T.z==NULL))
	{
		AddErr(MALLOCFAILTOPOLOGY);
		return T0;
	}
	for (i=0;i<N;i++)
	{
		T.x[i]=x[i];
		T.y[i]=y[i];
		T.z[i]=z[i];
	}
	T.N=N;
	T.T=Triangulate(T.x, T.y, T.N, &T.Nt);
	if (ssdp_error_state)
		return T0;
	T.P=InitTree(T.T, T.Nt, T.x, T.y, N);
	Print(VVERBOSE, "\ntriangulation: %d triangles\n", T.Nt);
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n");
	return T;
}

void free_topo (topology *T)
{
	if (T)
	{
		if (T->x)
			free(T->x);
		if (T->y)
			free(T->y);
		if (T->z)
			free(T->z);
		if (T->T)
			free(T->T);
		T->x=NULL;
		T->y=NULL;
		T->z=NULL;	
		T->T=NULL;	
		T->N=0;	
		if (T->P)
		{
			FreeTree(T->P);
			free(T->P);
			T->P=NULL;
		}
	}
}

// should export norm somehow as we may want to rotate the pv panel accordingly
double SampleTopo(double x, double y, topology *T, sky_pos *sn)
{ // dangerous, no check for extrapolation
	int n;
	vec nn, a, b, c, v1, v2;
	double D,z, l;
	n=treesearch(T->T, T->P, T->x, T->y, T->Nt, x, y);
	a.x=T->x[T->T[n].i];
	a.y=T->y[T->T[n].i];
	a.z=T->z[T->T[n].i];
	b.x=T->x[T->T[n].j];
	b.y=T->y[T->T[n].j];
	b.z=T->z[T->T[n].j];
	c.x=T->x[T->T[n].k];
	c.y=T->y[T->T[n].k];
	c.z=T->z[T->T[n].k];
	v1=diff(a,b);
	v2=diff(a,c);
	nn=cross(v2,v1);
	l=norm(nn);
	
	/* the delaunay code should have arranged the vertices such that we do not need this check
	 if (nn.z<0)
		nn=scalevec(nn,-1);
	*/
	if (nn.z<1e-10*l) // avoid div 0
		nn.z=1e-10*l;
	//A x + B y +C z = D
	D=nn.x*a.x+nn.y*a.y+nn.z*a.z;
	z=(D-nn.x*x-nn.y*y)/nn.z;
	if (sn)
		(*sn)=vecdir(nn);
	return z;	
}
topology CreateRandomTopology(double dx, double dy, double dz, double fN, int N)
{
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double *x, *y, *z;
	topology T;
	topology T0={NULL,NULL,NULL,0,NULL,0,NULL};
	int i, n;
	VERB verbstate;
	verbstate=ssdp_verbosity;
	ssdp_verbosity=QUIET;
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--Create Random Topology\t\t\t");
	Print(VVERBOSE, "\ncreating topology: %d points\n", N);
	fN=fabs(fN);
	if (fN>1)
		fN=1/fN;
		
	n=3+(int)round(fN*((double)N));
	
	// ERRORFLAG MALLOCFAILRANDTOPOLOGY  "Error memory allocation failed in generating a random topology"
	x=malloc(n*sizeof(double));
	y=malloc(n*sizeof(double));
	z=malloc(n*sizeof(double));
	if ((x==NULL)||(y==NULL)||(z==NULL))
	{
		AddErr(MALLOCFAILRANDTOPOLOGY);
		return T0;
	}
	srand(time(NULL));
	xmin=-dx/2;
	xmax=dx/2;
	ymin=-dy/2;
	ymax=dy/2;
	zmin=-dz/2;
	zmax=dz/2;
	for (i=0;i<n;i++)
	{
		x[i]=xmin+(xmax-xmin)*((double)rand())/RAND_MAX;
		y[i]=ymin+(ymax-ymin)*((double)rand())/RAND_MAX;
		z[i]=zmin+(zmax-zmin)*((double)rand())/RAND_MAX;
	}
	T=MakeTopology(x, y, z, n);
	if (ssdp_error_state)
	
	x=realloc(x,N*sizeof(double));
	y=realloc(y,N*sizeof(double));
	z=realloc(z,N*sizeof(double));
	if ((x==NULL)||(y==NULL)||(z==NULL))
	{
		AddErr(MALLOCFAILRANDTOPOLOGY);
		return T0;
	}
	for (i=0;i<N;i++)
	{
		x[i]=xmin+(xmax-xmin)*((double)rand())/RAND_MAX;
		y[i]=ymin+(ymax-ymin)*((double)rand())/RAND_MAX;
		z[i]=SampleTopo(x[i], y[i], &T, NULL);
		if (z[i]<zmin)
			z[i]=zmin;
		if (z[i]>zmax)
			z[i]=zmax;
	}
	free_topo (&T);
	T=MakeTopology(x, y, z, N);
	free(x);
	free(y);
	free(z);
	ssdp_verbosity=verbstate;
	Print(VVERBOSE, "\ntriangulation: %d triangles\n", T.Nt);
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n");
	return T;
}



typedef struct horizon {
	int N;
	double *zen;
	double *azi;
} horizon;

horizon InitHorizon(int Nz)
{
	horizon H;
	horizon H0={0,NULL,NULL};
	int i;
	double astep;
	H.N=6*Nz;
	astep=2*M_PI/((double)H.N);
	// ERRORFLAG MALLOCFAILHORIZON  "Error memory allocation failed horizon initialization"
	H.zen=malloc(H.N*sizeof(double));
	H.azi=malloc(H.N*sizeof(double));
	if ((H.zen==NULL)||(H.azi==NULL))
	{
		AddErr(MALLOCFAILHORIZON);
		return H0;
	}
	for (i=0;i<H.N;i++)
	{
		H.zen[i]=M_PI/2;
		H.azi[i]=((double)i)*astep;
	}
	return H;
}
void FreeHorizon(horizon *H)
{
	free(H->zen);
	free(H->azi);
	H->N=0;
	H->zen=NULL;
	H->azi=NULL;
}

void RizeHorizon(horizon *H, double azi1, double azi2, double zen)
{
	int i, j, k;
	double astep=2*M_PI/((double)H->N);
	if (zen>M_PI/2)
		return;
	i=(int)round(azi1/astep);
	j=(int)round(azi2/astep);
	if (i<0)
		i+=H->N;
	if (j<0)
		j+=H->N;
	i=i%H->N;
	j=j%H->N;
	if (abs(j-i)>H->N/2)
	{
		// i must be larger than j
		if (i<j)
		{
			k=i;
			i=j;
			j=k;
		}
	}
	else
	{
		// i must be smaller than j
		if (i>j)
		{
			k=i;
			i=j;
			j=k;
		}
	}
	for (;i%H->N<j;i++)
	{
		k=i%H->N;
		if (H->zen[k]>zen)
			H->zen[k]=zen;
	}	
}

/* MakeHorizon2: old routine kept for reference.
 * This routine is slightly slower as it relies more on atan2 operations than the new one
 * In this code three atan2's are used to find the azimuthal range of an element
 * In the new code the edges of the triangles are analyzed to determine which two corners
 * determine the azimuthal range and we can skip one atan2 (of 4 per element).
 */ 

void MakeHorizon2(horizon *H, topology *T, double xoff, double yoff, double zoff)
{
	int i;
	double d;
	double z, a1, a2, a3, w1,w2,w3;
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--MakeHorizon2\t\t\t");
	Print(VVERBOSE, "\ntopology: %d points\n", T->N);
	Print(VVERBOSE, "horizon: %d points\n", H->N);
	Print(VVERBOSE, "Computing Horizon\n");
	
	for (i=0;i<T->Nt;i++)
	{
		// compute sky position and diameter in radians
		
		d=sqrt((T->T[i].ccx-xoff)*(T->T[i].ccx-xoff)+(T->T[i].ccy-yoff)*(T->T[i].ccy-yoff));
		z=(T->z[T->T[i].i]+T->z[T->T[i].j]+T->z[T->T[i].k])/3;
		
		a1=atan2(T->y[T->T[i].i]-yoff,T->x[T->T[i].i]-xoff);
		a2=atan2(T->y[T->T[i].j]-yoff,T->x[T->T[i].j]-xoff);
		a3=atan2(T->y[T->T[i].k]-yoff,T->x[T->T[i].k]-xoff);
		w1=fabs(adiff(a1, a2));
		w2=fabs(adiff(a2, a3));
		w3=fabs(adiff(a1, a3));		
		
		if ((w1>w2)&&(w1>w3))
			RizeHorizon(H, a1, a2, M_PI/2-atan2(z-zoff,d));
		else if ((w2>w1)&&(w2>w3))
			RizeHorizon(H, a2, a3, M_PI/2-atan2(z-zoff,d));
		else if ((w3>w1)&&(w3>w2))
			RizeHorizon(H, a1, a3, M_PI/2-atan2(z-zoff,d));
	}
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n\n");
	
}

// alternative MakeHorizon code relying on the triangulation producing right handed triangles
// This code computes onl;y two atan2's to get the azimuth range of a triangle instead of three
// speeds up the code maybe 15% or so
double pcross(double ax, double ay, double bx, double by, double cx, double cy) 
{ 
	return (bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
} 

int EdgeVis(double px, double py, double ax, double ay, double bx, double by)
{
	return (pcross(px,py,ax,ay,bx,by)<0);
}


#define AX x[T.i]
#define AY y[T.i]
#define BX x[T.j]
#define BY y[T.j]
#define CX x[T.k]
#define CY y[T.k]
int TriangleAziRange(triangles T, double *x, double *y, double xoff, double yoff, double *a1, double *a2)
{
	// assume right handed triangles
	
	// fprintf(stderr, "Winding %f\n", pcross(AX, AY, BX, BY, CX, CY));
	if (EdgeVis(xoff,yoff,AX,AY,BX,BY)&&(EdgeVis(xoff,yoff,BX,BY,CX,CY))&&(EdgeVis(xoff,yoff,CX,CY,AX,AY))) // right winding with each edge
	{
		(*a1)=-M_PI;
		(*a2)=M_PI;
		return 0; // we are inside the triangle
	}
	if (pcross(xoff,yoff,AX,AY,BX,BY)*(pcross(xoff,yoff,BX,BY,CX,CY))>0) // consitently left or right winding from A via B to C
	{
		// p1 p3
		(*a1)=atan2(AY,AX);
		(*a2)=atan2(CY,CX);
		return 1;
	}
	if (pcross(xoff,yoff,BX,BY,CX,CY)*(pcross(xoff,yoff,CX,CY,AX,AY))>0) // consitently left or right winding from B via C to A
	{
		// p2 p1
		(*a1)=atan2(BY,BX);
		(*a2)=atan2(AY,AX);
		return 1;
	}
	if (pcross(xoff,yoff,CX,CY,AX,AY)*(pcross(xoff,yoff,AX,AY,BX,BY))>0) // consitently left or right winding from C via A to B
	{
		// p2 p3
		(*a1)=atan2(CY,CX);
		(*a2)=atan2(BY,BX);
		return 1;
	}
	// ERRORFLAG TRIANGLEMESS  "Error cannot figure out the azimuth range of a triangle"
	AddErr(TRIANGLEMESS);
	(*a1)=-M_PI;
	(*a2)=M_PI;
	return 0;
}
#undef AX 
#undef AY 
#undef BX 
#undef BY 
#undef CX 
#undef CY 


void MakeHorizon(horizon *H, topology *T, double xoff, double yoff, double zoff)
{
	int i;
	double d;
	double z, a1, a2;
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--MakeHorizon2\t\t\t");
	Print(VVERBOSE, "\ntopology: %d points\n", T->N);
	Print(VVERBOSE, "horizon: %d points\n", H->N);
	Print(VVERBOSE, "Computing Horizon\n");
	
	for (i=0;i<T->Nt;i++)
	{
		// compute sky position and diameter in radians
		
		d=sqrt((T->T[i].ccx-xoff)*(T->T[i].ccx-xoff)+(T->T[i].ccy-yoff)*(T->T[i].ccy-yoff));
		z=(T->z[T->T[i].i]+T->z[T->T[i].j]+T->z[T->T[i].k])/3;
		if (TriangleAziRange(T->T[i], T->x, T->y, xoff, yoff, &a1, &a2)) // the horizon routine is not equipped to put a roof over the PV panel...
			RizeHorizon(H, a1, a2, M_PI/2-atan2(z-zoff,d));
	}
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n\n");
	
}

double MinZentith(horizon *H)
{
	int i;
	double min;
	min=H->zen[0];
	for (i=1;i<H->N;i++)
		if (min>H->zen[i])
			min=H->zen[i];
	return min;
}

int BelowHorizon(horizon *H, sky_pos p)
{
	int i,j;
	double astep=2*M_PI/((double)H->N);
	double a1,a2;
	i=(int)p.a/astep;
	if (i<0)
		i+=H->N;
	i=i%H->N;
	if (p.a<0)
		p.a+=M_PI*2;
	while (adiff(i*astep, p.a)>0)
	{
		i--;
		if (i<0)
			i=H->N-1;
	}
	j=(i+1)%H->N;
	a1=adiff(p.a, ((double)i)*astep);
	a2=adiff(((double)j)*astep, p.a);
	return (p.z>(H->zen[i]*a2+H->zen[j]*a1)/(a1+a2));
}

/* takes a sky, topology anmd a transfer function.
 * creates a new transfer function which multiplies the old one with 0 for elements below the horizon 
 * if ST==NULL it creates a new transfer function with 1's for non masked elements
 */
sky_transfer MaskHorizon(sky_grid *sky, topology *T, sky_transfer *ST, double xoff, double yoff, double zoff) 
{
	int i, nz;
	horizon H;
	sky_transfer M={NULL,0};
	double minz;
	if (ST)
	{
		// ERRORFLAG SKYTRANSSKYMISMATCH  "Error: sky transfer function not compatible with provided sky"
		if (ST->N!=sky->N)
		{
			AddErr(SKYTRANSSKYMISMATCH);
			return M;
		}
	}
	H=InitHorizon(sky->Nz);
	if (ssdp_error_state)
		return M;
	MakeHorizon(&H, T, xoff, yoff, zoff);
	M=InitSkyTransfer(sky->N);
	minz=MinZentith(&H);
	nz=(int) floor(minz/((M_PI/2)/((double)sky->Nz)));
	if (ST)
	{
		for (i=3*(nz-1)*nz+1;i<sky->N;i++)
			if (BelowHorizon(&H, sky->P[i].p))
				M.t[i]=ST->t[i]*0.0;
	}
	else
	{
		for (i=3*(nz-1)*nz+1;i<sky->N;i++)
			if (BelowHorizon(&H, sky->P[i].p))
				M.t[i]=0.0;
	}
	FreeHorizon(&H);
	return M;
}

