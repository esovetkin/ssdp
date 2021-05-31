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
#include <math.h>
#include "vector.h"
#include "sky_dome.h"
#include "error.h"
#include "print.h"

// helper routines for indexing hex meshes
int i_sqrt(int x)	// integer sqrt
{
	if(x<=0)
		return 0;
	else
		return (int)sqrt((float)(x));
}
int NNZ(int n)	// number of elements in a mesh given a number of levels (zenith discretizations)
{
	return 3*n*(n+1)+1;
}
int NZN(int n)	// inverse of above, i.e. give it a index and it returns the level (zenith) index
{
	return (1+(i_sqrt(12*(n-1)+9)-3)/6);
}

// compute angles for a given index, i, and mesh with Nz zenith disdcretizations
sky_pos GridPos(int Nz, int i)
{
	sky_pos r;
	int Na;
	int nz;
	if (i<=0)
	{
		r.a=0;
		r.z=0;
		return r;
	}
	if (i>NNZ(Nz))
	{
		r.a=0;
		r.z=M_PI/2.0;
		return r;
	}
	nz=NZN(i);
	Na=nz*6;
	r.z=((double)nz)*M_PI/2.0/((double)Nz);
	r.a=fmod(2*((double)(i-NNZ(nz-1)))*M_PI/(double)Na, 2*M_PI);
	return r;
}

// compute index for a point p in a mesh with Nz zenith discretizations
// is a bit approximate though, you may end up in the patch next to the one you want.
int GridIndex(int Nz, sky_pos p)
{
	int nz,Na,i;
	double zstep, astep;
	if (p.z<=0)
		return 0;
	if (p.z>M_PI/2)
		return NNZ(Nz)-1;
		
	zstep=M_PI/2.0/((double)Nz);
	nz=(int)round(p.z/zstep);
	if (nz==0)
		return 0;
	if (nz>Nz)
		nz=Nz;
	Na=nz*6; 
	i=NNZ(nz-1);	
	astep=2*M_PI/((double)Na);
	if (p.a>=0)
		i+=(int)round(p.a/astep);
	else
		i+=(int)round((p.a+2*M_PI)/astep);
	
	if (i>NNZ(Nz)-1)
		i=NNZ(Nz)-1;
	if (i<0)
		i=0;
	return i;
}


// routines to connect patches
int NextIsoL(int Nz, int index) // same level next
{
	int nz;
	if (index==0)
		return -1;
		
	nz=NZN(index);
	if (nz==NZN(index+1))
		return index+1;
	return index-nz*6+1;// last of this level, must go to first of this level
}

int PrevIsoL(int Nz, int index) // same level previous
{
	int nz;
	if (index==0)
		return -1;
		
	nz=NZN(index);
	if (nz==NZN(index-1))// previous at same level
		return index-1;
	return index+nz*6-1;// last of this level, must go to first of this level
}

void NextL(int Nz, int index, int *N) // next level
{
	int nz, i, aoff, n=0;
	nz=NZN(index); // this level
	i=NNZ(nz-1);   // start index of this level
	aoff=index-i;  // offset at this level
	
	if (index==0) // hexpatch 0 only has next level neighbors
	{
		N[n++]=1;
		N[n++]=2;
		N[n++]=3;
		N[n++]=4;
		N[n++]=5;
		N[n++]=6;
		N[n]=-1;
		return;
	}
	if (nz==Nz) // last level has no next
	{
		N[n]=-1;
		return;
	}
	// from the center of the grid there are 3 axis where the hexpatches are perfectly aligned.
	// this means that there are 6 azimuth angles that always repeat for every level
	// this identifies whether we are on one of the 6 major azimuth angles
	if (aoff%nz==0)
	{
		N[n++]=index+(NNZ(nz)-i)+aoff/nz;
		N[n]=N[n-1]+1;
		n++;
		N[n]=N[n-2]-1;
		n++;
		if (NZN(N[n-1])==nz) // this one is two levels down...
			N[n-1]=NNZ(nz+1)-1;
	}
	else
	{
		N[n++]=index+(NNZ(nz)-i)+aoff/nz;
		N[n]=N[n-1]+1;
		n++;
	}
	N[n]=-1;
}

void PrevL(int Nz, int index, int *N) // previous level
{
	int nz, i, aoff, n=0;
	nz=NZN(index); // this level
	i=NNZ(nz-1);   // start index of this level
	aoff=index-i;  // offset at this level
	
	if (index==0) // hexpatch 0 only has no previous level
	{
		N[n++]=-1;
		return;
	}
	// from the center of the grid there are 3 axis where the hexpatches are perfectly aligned.
	// this means that there are 6 azimuth angles that always repeat for every level
	// this identifies whether we are on one of the 6 major azimuth angles
	if (aoff%nz==0)
	{
		if (nz==1)
			N[n++]=0;
		else
			N[n++]=index-(nz-1)*6-aoff/nz;
		N[n++]=-1;
	}
	else
	{
		N[n++]=index-(nz-1)*6-aoff/nz-1;
		N[n]=N[n-1]+1;
		n++;
		if (NZN(N[n-1])==nz)
			N[n-1]=NNZ(nz-1);
	}
	N[n]=-1;
}


// for debugging meshing some routines to export connections
void Neighbors(int Nz, int index, int *N) // collect all neighbors
{
	int *n;
	n=N;
	(*n)=PrevIsoL(Nz, index);
	if(*n>=0)
		n++;
	
	(*n)=NextIsoL(Nz, index);
	if (*n>=0)
		n++;	
		
	PrevL(Nz, index, n);
	while (*n>=0)
		n++;	
		
	NextL(Nz, index, n);
	while (*n>=0)
		n++;	
}

// used for degugging, see if patches know their neighbors
void PlotConn(int Nz, int N[], int index)
{
	sky_pos p0,p;
	vec v0, v;
	int i;
	p0=GridPos(Nz, index);
	v0=unit(p0);
	i=0;
	while (N[i]>=0)
	{
		p=GridPos(Nz, N[i]);
		v=diff(unit(p), v0);
		// gnuplot compatible vectors (splot '...' u 1:2:3:4:5:6 w vectors)
		printf("%e %e %e %e %e %e\n", v0.x, v0.y, v0.z, v.x, v.y, v.z);
		i++;
	}
}	
void Connectivity(int Nz)
{
	int i;
	int N[7];
	for (i=0;i<NNZ(Nz);i++)
	{
		Neighbors(Nz, i, N);
		PlotConn(Nz,N,i);
	}
}


// findpatch routine
int FindPatch(sky_grid *sky, sky_pos p)
{
	int r;
	int i,j;
	int dmin, d, rmin;
	r=GridIndex(sky->Nz, p);
	rmin=r;
	// result of GridIndex(...) may be a bit approximate
	// look through the neigbors for the best patch
	dmin=(p.a-sky->P[r].p.a)*(p.a-sky->P[r].p.a)+(p.z-sky->P[r].p.z)*(p.z-sky->P[r].p.z);
	i=0;
	while ((j=sky->P[r].NL[i])>=0)
	{
		d=(p.a-sky->P[j].p.a)*(p.a-sky->P[j].p.a)+(p.z-sky->P[j].p.z)*(p.z-sky->P[j].p.z);
		if (d<dmin)
		{
			dmin=d;
			rmin=j;
		}
		i++;
	}
	i=0;
	while ((j=sky->P[r].PL[i])>=0)
	{
		d=(p.a-sky->P[j].p.a)*(p.a-sky->P[j].p.a)+(p.z-sky->P[j].p.z)*(p.z-sky->P[j].p.z);
		if (d<dmin)
		{
			dmin=d;
			rmin=j;
		}
		i++;
	}
	if ((j=sky->P[r].PI)>=0)
	{
		d=(p.a-sky->P[j].p.a)*(p.a-sky->P[j].p.a)+(p.z-sky->P[j].p.z)*(p.z-sky->P[j].p.z);
		if (d<dmin)
		{
			dmin=d;
			rmin=j;
		}
	}
	if ((j=sky->P[r].NI)>=0)
	{
		d=(p.a-sky->P[j].p.a)*(p.a-sky->P[j].p.a)+(p.z-sky->P[j].p.z)*(p.z-sky->P[j].p.z);
		if (d<dmin)
		{
			dmin=d;
			rmin=j;
		}
	}
	return rmin;	
}

// as the number of elements rises quadratically we better set a
// practical limit on the Nz values
// if you ever change the 1 Mio below you should consider modifying 
// ZENEPS in project.c.
#define MAXNZ NZN(1000000)

sky_grid InitSky(int Nz)
{
	sky_grid sky;
	sky_pos sun={0,0};// default sun, straight above
	int i;
	if (Nz>MAXNZ)
	{
		Print(WARNING,"Warning: number of zenith discretizations %d too large, using %d instead\n", Nz, MAXNZ);
		Nz=MAXNZ;
	}
	sky.N=NNZ(Nz);	
	sky.Nz=Nz;
	sky.sp=sun; // The default sun: is strainght above
	sky.sI=0;	// and pitch black
	sky.suni=0;
	// ERRORFLAG MALLOCFAILSKYDOME  "Error memory allocation failed in creating a sky-dome"
	if ((sky.P=malloc((sky.N+1)*sizeof(hexpatch)))==NULL)
	{
		AddErr(MALLOCFAILSKYDOME);
		sky.Nz=0;
		sky.N=0;
		return sky;
	}
	if ((sky.cosz=malloc((sky.N+1)*sizeof(double)))==NULL)
	{
		free(sky.P);
		sky.P=NULL;
		AddErr(MALLOCFAILSKYDOME);
		sky.Nz=0;
		sky.N=0;
		return sky;
	}
	sky.icosz=0;
	for (i=0;i<sky.N;i++)
	{
		sky.P[i].I=0;
		sky.P[i].p=GridPos(Nz, i);
		NextL(Nz, i, sky.P[i].NL);
		PrevL(Nz, i, sky.P[i].PL);
		sky.P[i].NI=NextIsoL(Nz, i);
		sky.P[i].PI=PrevIsoL(Nz, i);
		sky.cosz[i]=cos(sky.P[i].p.z);
		sky.icosz+=sky.cosz[i];
	}
	return sky;
}

void free_sky_grid(sky_grid *sky)
{
	if (sky->P)
		free(sky->P);
	sky->P=NULL;
	if (sky->cosz)
		free(sky->cosz);
	sky->cosz=NULL;	
	sky->Nz=0;
	sky->N=0;
}
