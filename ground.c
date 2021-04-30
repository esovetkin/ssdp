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
#include "trace.h"
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
	{
		a=nn;
		// in our coordinate system we go counter clockwise with 0 degrees pointing up (swap x and y)
		nn.x=nn.y;
		nn.y=a.x;
		(*sn)=vecdir(nn);
	}
	return z;	
}
topology CreateRandomTopology(double dx, double dy, double dz, double fN, int N)
{
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double *x, *y, *z;
	topology T;
	topology T0={NULL,NULL,NULL,0,NULL,0,NULL};
	int i, n;
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
	return T;
}


void RizeHorizon(horizon *H, double azi1, double azi2, double zen)
{
	int i, j, k;
	if (zen>M_PI/2)
		return;
	i=(int)(azi1/H->astep);
	j=(int)(azi2/H->astep);
	
	if (i<0)
		i+=H->N;
	if (j<0)
		j+=H->N;
		
	i=i%H->N;
	j=j%H->N;
	
	if (i>j)
	{
		k=i;
		i=j;
		j=k;
	}
	if (j-i>=H->N/2)
	{
		for (k=j;k<H->N+i;k++)
		{
			if (H->zen[k%H->N]>zen)
				H->zen[k%H->N]=zen;
		}
	}
	else
	{	
		for (k=i;k<=j;k++)
		{
			if (H->zen[k]>zen)
				H->zen[k]=zen;
		}
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

int TriangleAziRange(triangles T, double *x, double *y, double dz, double xoff, double yoff, double *a1, double *a2)
{
	// assume right handed triangles
	double dx, dy, l, dz2=dz*dz;
	// fprintf(stderr, "Winding %f\n", pcross(AX, AY, BX, BY, CX, CY));
	if (EdgeVis(xoff,yoff,AX,AY,BX,BY)&&(EdgeVis(xoff,yoff,BX,BY,CX,CY))&&(EdgeVis(xoff,yoff,CX,CY,AX,AY))) // right winding with each edge
	{
		(*a1)=-M_PI;
		(*a2)=M_PI;
		return 0; // we are inside the triangle
	}
	if (pcross(xoff,yoff,AX,AY,BX,BY)*(pcross(xoff,yoff,BX,BY,CX,CY))>0) // consitently left or right winding from A via B to C
	{
		/* the cartesian and angular systems used here orientate to how we generally use 
		 * maps and compasses this means that north is up (i.e. the y-axis points north) 
		 * and is at 0°. The x-axis points to the east which is at 90° (i.e. we go clockwise). 
		 * Thus in the following atan2's we swap the usual positions of the x and y arguments.
		 * Note that the definition of angles is also affected by the sunpos.c code.
		 */ 
		// p1 p3
		dx=AX-xoff;
		dy=AY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a1)=atan2(l*dx,l*dy); // I want north to be 0° and east
		
		dx=CX-xoff;
		dy=CY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a2)=atan2(l*dx,l*dy);
		return 1;
	}
	if (pcross(xoff,yoff,BX,BY,CX,CY)*(pcross(xoff,yoff,CX,CY,AX,AY))>0) // consitently left or right winding from B via C to A
	{
		// p2 p1
		dx=BX-xoff;
		dy=BY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a1)=atan2(l*dx,l*dy);
		
		dx=AX-xoff;
		dy=AY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a2)=atan2(l*dx,l*dy);
		return 1;
	}
	if (pcross(xoff,yoff,CX,CY,AX,AY)*(pcross(xoff,yoff,AX,AY,BX,BY))>0) // consitently left or right winding from C via A to B
	{
		// p2 p3
		dx=CX-xoff;
		dy=CY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a1)=atan2(l*dx,l*dy);
		
		dx=BX-xoff;
		dy=BY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a2)=atan2(l*dx,l*dy);
		
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

// take a topology and rize the horizon accordingly
void ComputeHorizon(horizon *H, topology *T, double xoff, double yoff, double zoff)
// in principle one can later add a second topology to an existing horizon, is this useful? probably not.
{
	int i;
	double d;
	double z, a1, a2;	
	for (i=0;i<T->Nt;i++)
	{
		// compute sky position and diameter in radians
		
		d=sqrt((T->T[i].ccx-xoff)*(T->T[i].ccx-xoff)+(T->T[i].ccy-yoff)*(T->T[i].ccy-yoff));
		z=(T->z[T->T[i].i]+T->z[T->T[i].j]+T->z[T->T[i].k])/3;
		if (TriangleAziRange(T->T[i], T->x, T->y, z-zoff, xoff, yoff, &a1, &a2)) // the horizon routine is not equipped to put a roof over the PV panel...
		{
			//fprintf(stderr, "%e %e %e %e %e\n", d, z-zoff, a1, a2,  atan2(d, z-zoff));
			RizeHorizon(H, a1, a2, atan2(d, z-zoff));
		}
	}	
}

horizon MakeHorizon(sky_grid *sky, topology *T, double xoff, double yoff, double zoff)
{
	horizon H={0,0,NULL};
	H=InitHorizon(sky->Nz);
	if (ssdp_error_state)
		return H;
	ComputeHorizon(&H, T, xoff, yoff, zoff);
	return H;
}


