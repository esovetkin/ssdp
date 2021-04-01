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
#include "project.h"
#include "error.h"
#include "print.h"
#include "pv-aoi.h"

AOI_Model_Data InitAOIModel(AOI_Model model, double ng, double nar , double *theta, double *effT, int N)
{
	AOI_Model_Data M;
	M.M=model;
	M.ng=ng;
	M.nar=nar;
	if (M.M==AOI_USER)
	{
		int i;
		// ERRORFLAG NODATAFORUSERAOIMODEL "Error: no data provided for the user AOI model"
		if ((!theta)||(!effT)||(N==0))
		{
			AddErr(NODATAFORUSERAOIMODEL);
			// fallback none
			M.M=AOI_NONE;
			M.N=0;
			M.theta=NULL;
			M.effT=NULL;
			return M;
		}
		// point to data, user must take care of its own arrays and not free them till we are done
		M.theta=theta;
		M.effT=effT;
		M.N=N;
	}
	else
	{
		M.theta=NULL;
		M.effT=NULL;
		M.N=0;
	}
	return M;
}
	

#define AOIEPS 1e-3
double EffectiveT(AOI_Model_Data *M, double z, double R0)
{
	switch (M->M)
	{
		case AOI_GLASS:
			return Transmission(1.0,M->ng,z)/R0;
		case AOI_GLASS_AR:
			return Transmission_ar(1.0,M->nar,M->ng,z)/R0;
		case AOI_USER:
		{
			int imin=0, imax=M->N-1, i;
			// find aoi in theta array
			if (M->effT[imin]>z)
				return 1.0;
			if (M->theta[imax]<z)
				return M->effT[imax]/M->theta[0];
			i=(imin+imax)/2;
			while ((i!=imin)&&(i!=imax))
			{
				if (M->theta[i]<z)
					imin=i;
				else
					imax=i;
				i=(imin+imax)/2;
			}
			if ((M->theta[imax]-M->theta[imin])>AOIEPS)
				return 	((z-M->theta[imin])*M->effT[imax]+(M->theta[imax]-z)*M->effT[imin])/(M->theta[imax]-M->theta[imin])/M->theta[0];
			else 
				return 	(M->effT[imax]+M->effT[imin])/2.0/M->theta[0];
		}		
		default:
			return 1/R0;
	}
}


double DiffusePlaneOfArray(sky_grid *sky, double tilt, double a, AOI_Model_Data *M, int mask)
{
	sky_pos axis, r;
	double POA=0, R0;
	int i;
	R0=EffectiveT(M, 0, 1);
	
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--DiffusePlaneOfArray\t\t");
	Print(VVERBOSE, "\nTilt:    %f\n", rad2degr(tilt));
	Print(VVERBOSE, "Azimuth: %f\n", rad2degr(a));
	axis.a=a+M_PI/2;// axis points perpendicular to a
	axis.z=M_PI/2;
	for (i=0;i<sky->N;i++)
	{
		if ((!sky->P[i].mask)||(!mask))
		{
			r=rrf(sky->P[i].p, axis, -tilt);  
			if ((r.z>=0)&&(r.z<=M_PI/2))                                        
				POA+=sky->P[i].I*cos(r.z)*EffectiveT(M, r.z, R0);
		}
	}
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n\n");
	
	return POA;
}
double DirectPlaneOfArray(sky_grid *sky, double tilt, double a, AOI_Model_Data *M, int mask)
{
	sky_pos axis, r;
	double POA=0, R0;
	R0=EffectiveT(M, 0, 1);
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--DirectPlaneOfArray\t\t");
	Print(VVERBOSE, "\nTilt:    %f\n", rad2degr(tilt));
	Print(VVERBOSE, "Azimuth: %f\n", rad2degr(a));
	axis.a=a+M_PI/2;// axis points perpendicular to a
	axis.z=M_PI/2;
	if ((!sky->smask)||(!mask))
	{
		r=rrf(sky->sp, axis, -tilt); 
		if ((r.z>=0)&&(r.z<=M_PI/2))       
			POA=sky->sI*cos(r.z)*EffectiveT(M, r.z, R0);
	}
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n\n");
	
	return POA;
}

double DiffuseHorizontal(sky_grid *sky, AOI_Model_Data *M, int mask)
{
	double POA=0, R0;
	int i;	
	R0=EffectiveT(M, 0, 1);
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--DiffuseHorizontal\t\t");
	Print(VVERBOSE, "\n");
	for (i=0;i<sky->N;i++)   
		if ((!sky->P[i].mask)||(!mask))
			POA+=sky->P[i].I*cos(sky->P[i].p.z)*EffectiveT(M, sky->P[i].p.z, R0);  
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n\n");  
	return POA;
}

double DirectHorizontal(sky_grid *sky, AOI_Model_Data *M, int mask)
{
	double POA=0, R0;
	R0=EffectiveT(M, 0, 1);
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--DirectHorizontal\t\t");
	Print(VVERBOSE, "\n");
	if ((!sky->smask)||(!mask))
		POA+=sky->sI*cos(sky->sp.z)*EffectiveT(M, sky->sp.z, R0);
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n\n");  
	return POA;
}

double POA_Albedo(sky_grid *sky, double albedo, double tilt, double a, AOI_Model_Data *M, int mask)
{
	double GHI, POA;
	sky_grid ground;
	sky_pos z={0,0};
	VERB val;
	Print(VVERBOSE, "********************************************************************************\n");
	Print(VERBOSE, "--POA_Albedo\t\t\t");
	Print(VVERBOSE, "\nAlbedo:  %f\n", albedo);
	Print(VVERBOSE, "Tilt:    %f\n", rad2degr(tilt));
	Print(VVERBOSE, "Azimuth: %f\n", rad2degr(a));
	val=ssdp_verbosity;
	ssdp_verbosity=QUIET;
	GHI=DiffuseHorizontal(sky, M, mask);
	GHI+=DirectHorizontal(sky, M, mask);
	ground=InitSky(sky->Nz); // this should not be necessary, so some math
	if (ssdp_error_state)
		return 0;
	UniformSky(&ground, z, albedo*GHI, albedo*GHI);
	POA=DiffusePlaneOfArray(&ground, tilt+M_PI, a, M, mask);
	free_sky_grid(&ground);	
	ssdp_verbosity=val;
	Print(VERBOSE, "Done\n");
	Print(VVERBOSE, "********************************************************************************\n\n");
	return POA;
}


void POA_to_SurfaceNormal(double *tilt, double *a, sky_pos sn)
{
	sky_pos axis, r;
	axis.a=sn.a;
	axis.z=M_PI/2;
	r.a=(*a);
	r.z=(*tilt);
	r=rrf(r, axis, sn.z);
	(*a)=r.a;
	(*tilt)=r.z;
}
	
	
