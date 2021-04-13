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
#define ZENEPS 1e-4
sky_transfer POA_Sky_Transfer(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, sky_transfer *ST)
{
	sky_transfer T={NULL,0};
	sky_pos axis, r;
	int i;
	char rot=1;
	double R0;
	R0=EffectiveT(M, 0, 1);
	if (ST)
	{
		if (ST->N!=sky->N)
		{
			AddErr(SKYTRANSSKYMISMATCH);// label defined in ground.c
			return T;
		}
	}
	T=InitSkyTransfer(sky->N);
	if (fabs(pn.z)<ZENEPS)
		rot=0;
	else
	{
		axis.a=pn.a+M_PI/2;// axis points perpendicular to a
		axis.z=M_PI/2;
	}
	
	if (ST)
	{
		for (i=0;i<sky->N;i++)
		{
			if (rot)
				r=rrf(sky->P[i].p, axis, -pn.z);  
			else
				r=sky->P[i].p;
			if ((r.z>=-M_PI/2)&&(r.z<=M_PI/2))                                        
				T.t[i]=ST->t[i]*cos(r.z)*EffectiveT(M, r.z, R0);
			else
				T.t[i]=0;
		}
	}
	else
	{
		for (i=0;i<sky->N;i++)
		{
			if (rot)
				r=rrf(sky->P[i].p, axis, -pn.z);  
			else
				r=sky->P[i].p;
			if ((r.z>=-M_PI/2)&&(r.z<=M_PI/2))  
				T.t[i]=cos(r.z)*EffectiveT(M, r.z, R0);
			else
				T.t[i]=0;
		}
	}
	return T;
}

double DiffusePlaneOfArray(sky_grid *sky, sky_transfer *T)
{
	double POA=0;
	int i;
	for (i=0;i<sky->N;i++)                            
			POA+=sky->P[i].I*T->t[i];
	return POA;
}

double DirectPlaneOfArray(sky_grid *sky, sky_transfer *T)
{
	if (sky->suni>=0)
		return sky->sI*T->t[sky->suni]; // we approximate a bit and use the transfer function of the sun patch for the direct contribution
	return 0;
}

sky_transfer POA_Albedo_Transfer(sky_grid *sky, sky_pos pn, AOI_Model_Data *M, sky_transfer *ST)
{
	sky_transfer T={NULL,0};
	sky_transfer Th;
	double g=0;
	sky_grid ground;
	int i;
	sky_pos n={0,0};
	
	if (ST)
	{
		if (ST->N!=sky->N)
		{
			AddErr(SKYTRANSSKYMISMATCH);// label defined in ground.c
			return T;
		}
	}
	/* compute the sky transfer function to the horizontal plane */ 
	Th=POA_Sky_Transfer(sky, n, M, ST);
	ground=InitSky(sky->Nz); 
	if (ssdp_error_state)
	{
		free_sky_grid(&ground);
		FreeSkyTransfer(&Th);
		return T;
	}
	pn.z+=M_PI;
	pn.z=fmod(pn.z,2*M_PI);
	T=POA_Sky_Transfer(&ground, pn, M, NULL);
	
	free_sky_grid(&ground);
	for (i=0;i<T.N;i++)
		g+=T.t[i];
	g=g/((double)T.N);
		
	for (i=0;i<T.N;i++)	
		T.t[i]=g*Th.t[i];
		
	FreeSkyTransfer(&Th);
	return T;
}

double POA_Albedo(sky_grid *sky, double albedo, sky_transfer *T)
{
	double POA=0;
	int i;
	for (i=0;i<sky->N;i++)                            
		POA+=sky->P[i].I*albedo*T->t[i];
	return POA;
}


void POA_to_SurfaceNormal(sky_pos *pn, sky_pos sn)
{
	sky_pos axis;
	
	axis.z=M_PI/2;
	axis.a=pn->a-M_PI/2;
	(*pn)=rrf(sn, axis, pn->z);	
}
	
	
