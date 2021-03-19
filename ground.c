#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "sky_dome.h"
#include "ground.h"



void free_topo (topology *T)
{
	free(T->points);
	T->N=0;	
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


void MakeHorizon(sky_grid *sky, topology T, double height) 
{
	int i;
	sky_pos p;
	double d, W;
	for (i=0;i<T.N;i++)
	{
		// compute sky position and diameter in radians
		d=sqrt(T.points[i].x*T.points[i].x+T.points[i].y*T.points[i].y);
		p.z=M_PI/2-fabs(atan2((T.points[i].z-height),d));
		p.a=atan2(T.points[i].y,T.points[i].x);
		W=2*atan(T.d/(2*d));
		if (p.z>M_PI/2)
			p.z=M_PI/2;
		UpdateHorizon(sky, p, W);
	}
}

void ClearHorizon(sky_grid *sky) 
{
	int i;
	for (i=0;i<sky->N;i++)
		sky->P[i].mask=0;
	sky->smask=0;
}

