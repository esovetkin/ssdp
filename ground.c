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
#include "ground.h"


void free_topo (topology *T)
{
	free(T->x);
	free(T->y);
	free(T->z);
	T->x=NULL;
	T->y=NULL;
	T->z=NULL;	
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


void MakeHorizon(sky_grid *sky, topology T, double xoff, double yoff, double zoff) 
{
	int i;
	sky_pos p;
	double d, W;
	for (i=0;i<T.N;i++)
	{
		// compute sky position and diameter in radians
		d=sqrt((T.x[i]-xoff)*(T.x[i]-xoff)+(T.y[i]-yoff)*(T.y[i]-yoff));
		p.z=M_PI/2-fabs(atan2(T.z[i]-zoff,d));
		p.a=atan2(T.y[i]-yoff,T.x[i]-xoff);
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

