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
#include "sky_model.h"

double DiffusePlaneOfArray(sky_grid sky, double tilt, double a, int mask)
{
	sky_pos axis, r;
	double POA=0;
	int i;
	
	axis.a=a+M_PI/2;// axis points perpendicular to a
	axis.z=M_PI/2;
	for (i=0;i<sky.N;i++)
	{
		if ((!sky.P[i].mask)||(!mask))
		{
			r=rrf(sky.P[i].p, axis, -tilt);  
			if ((r.z>=0)&&(r.z<=M_PI/2))                                        
				POA+=sky.P[i].I*cos(r.z);
		}
	}
	
	return POA;
}
double DirectPlaneOfArray(sky_grid sky, double tilt, double a, int mask)
{
	sky_pos axis, r;
	double POA=0;
	if ((!sky.smask)||(!mask))
	{
		r=rrf(sky.sp, axis, -tilt); 
		if ((r.z>=0)&&(r.z<=M_PI/2))       
			POA=sky.sI*cos(r.z);
	}
	return POA;
}

double DiffuseHorizontal(sky_grid sky, int mask)
{
	double POA=0;
	int i;	
	for (i=0;i<sky.N;i++)   
		if ((!sky.P[i].mask)||(!mask))
			POA+=sky.P[i].I*cos(sky.P[i].p.z);    
	return POA;
}

double DirectHorizontal(sky_grid sky, int mask)
{
	double POA=0;
	if ((!sky.smask)||(!mask))
		POA+=sky.sI*cos(sky.sp.z);
	return POA;
}

double POA_Albedo(sky_grid sky, double albedo, double tilt, double a, int mask)
{
	double GHI, POA;
	sky_grid ground;
	sky_pos z={0,0};
	GHI=DiffuseHorizontal(sky, mask);
	GHI+=DirectHorizontal(sky, mask);
	ground=InitSky(sky.Nz); // this should not be necessary but I have to so some math
	UniformSky(&ground, z, albedo*GHI, albedo*GHI);
	POA=DiffusePlaneOfArray(ground, tilt+M_PI, a, mask);
	free_sky_grid(&ground);	
	return POA;
}
