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
double PlaneOfArray(sky_grid sky, double tilt, double a, int mask)
{
	sky_pos axis, r;
	double POA;
	POA=DiffusePlaneOfArray(sky, tilt, a, mask);	
	if ((!sky.smask)||(!mask))
	{
		r=rrf(sky.sp, axis, -tilt); 
		if ((r.z>=0)&&(r.z<=M_PI/2))       
			POA+=sky.sI*cos(r.z);
	}
	return POA;
}
double DiffusteHorizontal(sky_grid sky, int mask)
{
	double POA=0;
	int i;	
	for (i=0;i<sky.N;i++)   
		if ((!sky.P[i].mask)||(!mask))
			POA+=sky.P[i].I*cos(sky.P[i].p.z);    
	return POA;
}
double GlobalHorizontal(sky_grid sky, int mask)
{
	double POA;
	POA=DiffusteHorizontal(sky, mask);
	if ((!sky.smask)||(!mask))
		POA+=sky.sI*cos(sky.sp.z);
	return POA;
}

double POA_Albedo(sky_grid sky, double albedo, double tilt, double a, int mask)
{
	double GHI, POA;
	sky_grid ground;
	sky_pos z={0,0};
	GHI=GlobalHorizontal(sky, mask);
	ground=InitSky(sky.Nz); // this should not be necessary but I have to so some math
	UniformSky(&ground, z, albedo*GHI, albedo*GHI);
	POA=DiffusePlaneOfArray(ground, tilt+M_PI, a, mask);
	free_sky_grid(&ground);	
	return POA;
}
