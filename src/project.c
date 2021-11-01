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
#include "trace.h"
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
			// binary search of aoi in theta array
			if (M->theta[imin]>z)
				return 1.0;
			if (M->theta[imax]<z)
				return M->effT[imax]/M->effT[0];
			i=(imin+imax)/2;
			while ((i!=imin)&&(i!=imax))
			{
				if (M->theta[i]<z)
					imin=i;
				else
					imax=i;
				i=(imin+imax)/2;
			}
			return 	((z-M->theta[imin])*M->effT[imax]+(M->theta[imax]-z)*M->effT[imin])/(M->theta[imax]-M->theta[imin])/M->theta[0];
		}		
		default:
			return 1/R0;
	}
}
// in sky_dome.c the maximum number af sky patches is set to 1 Mio, which constitutes 577 
// zenith steps and thus about 2.7e-3 radians. Thus a rotation in zenith over 
// considerably less than this cannot be adequately resolved, hence 2.7e-4 rad
#define ZENEPS 2.7e-4 
void POA_Sky_Transfer(sky_grid *sky, sky_transfer *T, sky_pos pn, AOI_Model_Data *M)
{
	sky_pos axis, r;
	int i;
	char rot=1;
	double R0;
	R0=EffectiveT(M, 0, 1);
	if (T->N!=sky->N)
	{
		// ERRORFLAG SKYTRANSSKYMISMATCH "Error: sky_transfer and sky_grid data structures size mismatch"
		AddErr(SKYTRANSSKYMISMATCH);
		return;
	}
	if (fabs(pn.z)<ZENEPS)
		rot=0;
	else
	{
		axis.a=pn.a+M_PI/2;// rotation axis points perpendicular in azimuth
		axis.z=M_PI/2;// rotation axis in ground plane
	}
	
	for (i=0;i<sky->N;i++)
	{
		if (rot)
		{
			r=rrf(sky->P[i].p, axis, -pn.z);  
			if ((r.z>=-M_PI/2)&&(r.z<=M_PI/2))                                        
				T->t[i]=T->t[i]*cos(r.z)*EffectiveT(M, r.z, R0);
			else
				T->t[i]=0;
		}
		else              
			T->t[i]=T->t[i]*sky->cosz[i]*EffectiveT(M, sky->P[i].p.z, R0);
	}
}

double DiffusePlaneOfArray(sky_grid *sky, sky_transfer *T)
{
	double POA=0;
	int i;
	for (i=0;i<sky->N;i++)                            
		POA+=sky->P[i].I*T->t[i];
	return POA;
}

double DirectPlaneOfArray(sky_grid *sky, horizon *H, sky_pos pn, AOI_Model_Data *M)
{
	sky_pos axis, r;
	char rot=1;
	double R0;
	if (sky->suni<0) // sun not in sky
		return 0;
	if (BelowHorizon(H, sky->sp))
		return 0;
		
	R0=EffectiveT(M, 0, 1);
	if (fabs(pn.z)<ZENEPS)
		rot=0;
	else
	{
		axis.a=pn.a+M_PI/2;// rotation axis points perpendicular in azimuth
		axis.z=M_PI/2;// rotation axis in ground plane
	}
	
	if (rot)
	{
		r=rrf(sky->sp, axis, -pn.z);    
		if ((r.z>=-M_PI/2)&&(r.z<=M_PI/2))                                 
			return sky->sI*cos(r.z)*EffectiveT(M, r.z, R0);
		else
			return 0;
	}
	else
		return sky->sI*cos(sky->sp.z)*EffectiveT(M, sky->sp.z, R0);
	return 0;
}

double POA_Albedo_Transfer(sky_grid *sky, sky_pos pn, AOI_Model_Data *M)
{
	sky_transfer Th;
	double g=0;
	sky_grid ground;
	int i;
	
	/* compute the sky transfer function to the horizontal plane */ 	
	ground=InitSky(sky->Nz); 
	if (ssdp_error_state)
	{
		free_sky_grid(&ground);
		return 0;
	}
	pn.z+=M_PI;
	pn.z=fmod(pn.z,2*M_PI);
	Th=InitSkyTransfer(sky->N);	
	POA_Sky_Transfer(&ground,&Th, pn, M);	
	free_sky_grid(&ground);
	
	for (i=0;i<Th.N;i++)
		g+=Th.t[i];
	g=g/((double)Th.N);
	return g;
}
// diffuse light albedo transfer : albedo*g*sky->cosz[i];
// direct light albedo transfer : albedo*g*cos(sky->sp.z);


void POA_to_SurfaceNormal(sky_pos *pn, sky_pos sn)
{
	sky_pos axis;
	
	axis.z=M_PI/2;
	axis.a=pn->a-M_PI/2;
	(*pn)=rrf(sn, axis, -pn->z);	// which way to rotate?
}
	
	
