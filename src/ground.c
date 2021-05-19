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
#include "config.h"
#include "fatan2.h"
// sort triangle list by height

int tcomp(const void *a, const void *b)
{ 
	if (((triangles *)a)->ccz<((triangles *)b)->ccz)
		return -1;
	if (((triangles *)a)->ccz>((triangles *)b)->ccz)
		return 1;
	return 0;
} 
void tsort(triangles *T, int Nt)
{ 
	qsort(T, Nt, sizeof(triangles), tcomp);
}

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
	for (i=0;i<T.Nt;i++)
		T.T[i].ccz=(T.z[T.T[i].i]+T.z[T.T[i].j]+T.z[T.T[i].k])/3;
	tsort(T.T, T.Nt);
	
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
topology CreateRandomTopology(double dx, double dy, double dz, int N1, int N2)
{
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double *x, *y, *z;
	topology T;
	topology T0={NULL,NULL,NULL,0,NULL,0,NULL};
	int i, n;
	if (N1>N2)
	{
		n=N1;
		N1=N2;
		N2=n;
	}
	if (N1<3)
	{
		N1=3;
		N2+=3;
	}
	// ERRORFLAG MALLOCFAILRANDTOPOLOGY  "Error memory allocation failed in generating a random topology"
	x=malloc(N2*sizeof(double));
	y=malloc(N2*sizeof(double));
	z=malloc(N2*sizeof(double));
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
	for (i=0;i<N1;i++)
	{
		x[i]=xmin+(xmax-xmin)*((double)rand())/RAND_MAX;
		y[i]=ymin+(ymax-ymin)*((double)rand())/RAND_MAX;
		z[i]=zmin+(zmax-zmin)*((double)rand())/RAND_MAX;
	}
	T=MakeTopology(x, y, z, N1);
	//if (ssdp_error_state)
	for (i=0;i<N2;i++)
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
	T=MakeTopology(x, y, z, N2);
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
	if (zen<0)
		return;
	i=(int)floor(azi1/H->astep);
	j=(int)floor(azi2/H->astep);
	
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
		j=j+1;
		for (k=j;k<H->N;k++)
		{
			if (H->zen[k]>zen)
				H->zen[k]=zen;
		}
		for (k=0;k<i;k++)
		{
			if (H->zen[k]>zen)
				H->zen[k]=zen;
		}
	}
	else
	{	
		i=i+1;
		for (k=i;k<=j;k++)
		{
			if (H->zen[k]>zen)
				H->zen[k]=zen;
		}
	}
}

// MakeHorizon relies on the triangulation producing right handed triangles
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
		(*a1)=ATAN2(l*dx,l*dy); // I want north to be 0° and east
		
		dx=CX-xoff;
		dy=CY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a2)=ATAN2(l*dx,l*dy);
		return 1;
	}
	if (pcross(xoff,yoff,BX,BY,CX,CY)*(pcross(xoff,yoff,CX,CY,AX,AY))>0) // consitently left or right winding from B via C to A
	{
		// p2 p1
		dx=BX-xoff;
		dy=BY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a1)=ATAN2(l*dx,l*dy);
		
		dx=AX-xoff;
		dy=AY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a2)=ATAN2(l*dx,l*dy);
		return 1;
	}
	if (pcross(xoff,yoff,CX,CY,AX,AY)*(pcross(xoff,yoff,AX,AY,BX,BY))>0) // consitently left or right winding from C via A to B
	{
		// p2 p3
		dx=CX-xoff;
		dy=CY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a1)=ATAN2(l*dx,l*dy);
		
		dx=BX-xoff;
		dy=BY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a2)=ATAN2(l*dx,l*dy);
		
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
void ComputeHorizon(horizon *H, topology *T, double minzen, double xoff, double yoff, double zoff)
// the minzen parameter can save alot of work
// it specifies the minimum zenith angle to consider a triangle for the horizon
// I would set it to 0.5 times the zenith step in the sky. Especially for large topographies it reduces the amount of work considerably as far away triangles are less likely of consequence
{
	int i;
	double d;
	double z, a1, a2;	
	double r;
	r=tan(minzen);// compute threshold height over distance ratio
	i=T->Nt-1;
	while ((i>=0)&&(T->T[i].ccz>zoff))
	{
		// compute sky position and diameter in radians
		d=sqrt((T->T[i].ccx-xoff)*(T->T[i].ccx-xoff)+(T->T[i].ccy-yoff)*(T->T[i].ccy-yoff));
		// need to sort the triangles, then we do not need this test and we do not need to go through all triangles
		if ((T->T[i].ccz-zoff)/d>r) // do not compute anything for triangles below the zenith threshold
			if (TriangleAziRange(T->T[i], T->x, T->y, T->T[i].ccz-zoff, xoff, yoff, &a1, &a2)) // the horizon routine is not equipped to put a roof over the PV panel...
				RizeHorizon(H, a1, a2, d/(T->T[i].ccz-zoff));
		i--;
	}	
}


horizon MakeHorizon(sky_grid *sky, topology *T, double xoff, double yoff, double zoff)
{
	horizon H={0,0,NULL};
	H=InitHorizon(sky->Nz);
	if (ssdp_error_state)
		return H;
	ComputeHorizon(&H, T, M_PI/4.0/((double)sky->Nz), xoff, yoff, zoff);
	AtanHorizon(&H);
	return H;
}


