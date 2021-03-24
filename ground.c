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
#include <math.h>
#include "vector.h"
#include "sky_dome.h"
#include "delaunay.h"
#include "ground.h"
#include "util.h"


topology MakeTopology(double *x, double *y, double *z, int N)
{
	topology T;
	int i;
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--MakeTopology\t\t\t");
	Print(VVERBOSE, "\ncreating topology: %d points\n", N);
	T.x=malloc(N*sizeof(double));
	T.y=malloc(N*sizeof(double));
	T.z=malloc(N*sizeof(double));
	for (i=0;i<N;i++)
	{
		T.x[i]=x[i];
		T.y[i]=y[i];
		T.z[i]=z[i];
	}
	T.N=N;
	T.T=Triangulate(T.x, T.y, T.N, &T.Nt);
	T.P=InitTree(T.T, T.Nt, T.x, T.y, N);
	Print(VVERBOSE, "\ntriangulation: %d triangles\n", T.Nt);
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n");
	return T;
}

void free_topo (topology *T)
{
	free(T->x);
	free(T->y);
	free(T->z);
	free(T->T);
	T->x=NULL;
	T->y=NULL;
	T->z=NULL;	
	T->T=NULL;	
	T->N=0;	
	FreeTree(T->P);
	free(T->P);
	T->P=NULL;
}

// should export norm somehow as we may want to rotate the pv panel accordingly
double SampleTopo(double x, double y, topology T, sky_pos *sn)
{ // dangerous, no check for extrapolation
	int n;
	vec nn, a, b, c, v1, v2;
	double D,z, l;
	n=treesearch(T.T, T.P, T.x, T.y, T.Nt, x, y);
	a.x=T.x[T.T[n].i];
	a.y=T.y[T.T[n].i];
	a.z=T.z[T.T[n].i];
	b.x=T.x[T.T[n].j];
	b.y=T.y[T.T[n].j];
	b.z=T.z[T.T[n].j];
	c.x=T.x[T.T[n].k];
	c.y=T.y[T.T[n].k];
	c.z=T.z[T.T[n].k];
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
	
	x=malloc(n*sizeof(double));
	y=malloc(n*sizeof(double));
	z=malloc(n*sizeof(double));
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
	x=realloc(x,N*sizeof(double));
	y=realloc(y,N*sizeof(double));
	z=realloc(z,N*sizeof(double));
	for (i=0;i<N;i++)
	{
		x[i]=xmin+(xmax-xmin)*((double)rand())/RAND_MAX;
		y[i]=ymin+(ymax-ymin)*((double)rand())/RAND_MAX;
		z[i]=SampleTopo(x[i], y[i], T, NULL);
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

// this routine recursively asks all elements below element index
// which lie less than W/2 away from the azimuth in p
void MarkBelow(int index, sky_grid *sky, sky_pos p, double W)
{
	int i, k;
	
	i=0;
	while ((k=sky->P[index].NL[i])>=0) // check next level patches
	{
		if ((sky->P[k].p.z>sky->P[index].p.z)&&(fabs(adiff(sky->P[k].p.a,p.a))<W/2)&&(!sky->P[k].mask))
		{
			sky->P[k].mask=1; // mask patch
			MarkBelow(k, sky, p, W); // recursively mask patches below
		}
		i++;
	}
}

void UpdateHorizon(sky_grid *sky, sky_pos p, double W)
{
	int i, j, k, g;
	
	if ((fabs(adiff(sky->sp.a,p.a))>W/2)&&(sky->sp.z>p.z))
		sky->smask=1;
	
	i=FindPatch(*sky, p);
	if (fabs(adiff(sky->P[i].p.a,p.a))>W/2)
		return; // shade too small
		
	if (sky->P[i].p.z<p.z) // patch not below p.z
	{
		g=1;
		j=0;
		while (((k=sky->P[i].NL[j])>=0)&&g) // check next level patches
		{
			if ((sky->P[k].p.z>p.z)&&(fabs(adiff(sky->P[k].p.a,p.a))<W/2))
				g=0;
			j++;
		}		
		if (g)
			return; // cannot find a patch under p.z
	}
	
	// patch is found within azimuth range and below p.z 
	sky->P[i].mask=1;
	MarkBelow(i, sky, p, W); // recursively mask all patches below
	
	// previous patch at same zenith is N[0], next patch N[1]	
	j=i;
	while ((fabs(adiff(p.a,sky->P[j].p.a))<W/2)&&(sky->P[j].PI>0))
	{		
		j=sky->P[j].PI;
		if ((fabs(adiff(p.a,sky->P[j].p.a))<W/2)&&(sky->P[j].mask==0))
		{
			sky->P[j].mask=1;
			MarkBelow(j, sky, p, W); // recursively mask all patches below
		}
	}	
	j=i;
	while ((fabs(adiff(p.a,sky->P[j].p.a))<W/2)&&(sky->P[j].NI>0))
	{		
		j=sky->P[j].NI;
		if ((fabs(adiff(p.a,sky->P[j].p.a))<W/2)&&(sky->P[j].mask==0))
		{
			sky->P[j].mask=1;
			MarkBelow(j, sky, p, W); // recursively mask all patches below
		}
	}
	
}

// warning: zoff is absolute and not w.r.t. local ground level
void MakeHorizon(sky_grid *sky, topology T, double xoff, double yoff, double zoff) 
{
	int i;
	sky_pos p;
	double d, W;
	double z, a1, a2, a3, w1,w2,w3;
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--MakeHorizon\t\t\t");
	Print(VVERBOSE, "\ntopology: %d points\n", T.N);
	Print(VVERBOSE, "sky dome: %d patches\n", sky->N);
	Print(VVERBOSE, "Computing Horizon\n");
	
	for (i=0;i<T.Nt;i++)
	{
		// compute sky position and diameter in radians
		
		d=sqrt((T.T[i].ccx-xoff)*(T.T[i].ccx-xoff)+(T.T[i].ccy-yoff)*(T.T[i].ccy-yoff));
		z=(T.z[T.T[i].i]+T.z[T.T[i].j]+T.z[T.T[i].k])/3;
		p.z=M_PI/2-atan2(z-zoff,d);		
		a1=atan2(T.y[T.T[i].i]-yoff,T.x[T.T[i].i]-xoff);
		a2=atan2(T.y[T.T[i].j]-yoff,T.x[T.T[i].j]-xoff);
		a3=atan2(T.y[T.T[i].k]-yoff,T.x[T.T[i].k]-xoff);
		w1=fabs(adiff(a1, a2));
		w2=fabs(adiff(a2, a3));
		w3=fabs(adiff(a1, a3));		
		
		if ((w1>w2)&&(w1>w3))
		{
			W=w1;
			p.a=amean(a1, a2);
		}
		else if ((w2>w1)&&(w2>w3))
		{
			W=w2;
			p.a=amean(a2, a3);
		}
		else if ((w3>w1)&&(w3>w2))
		{
			W=w3;
			p.a=amean(a1, a3);
		}
		if (p.z<M_PI/2)
			UpdateHorizon(sky, p, W);
	}
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n\n");
}

void ClearHorizon(sky_grid *sky) 
{
	int i;
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--ClearHorizon\t\t\t");
	Print(VVERBOSE, "\n");
	for (i=0;i<sky->N;i++)
		sky->P[i].mask=0;
	sky->smask=0;
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n\n");
}


