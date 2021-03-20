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
#include "util.h"
	
void UniformSky(sky_grid *sky, sky_pos sun, double GHI, double DHI)
{
	double dhi0=0, dir;
	int i;
	Print(VERBOSE, "Creating uniform sky-dome\n");
	Print(VERBOSE, "GHI:           %e\nDHI:           %e\n\n", GHI, DHI);
	Print(VERBOSE, "Solar Zenith:  %f\nSolar Azimuth: %f\n", rad2degr(sun.z), rad2degr(sun.a));
	
	// sum intensity cos(z) product to normalize intensities
	for (i=0;i<sky->N;i++)
		dhi0+=cos(sky->P[i].p.z);
	dhi0=DHI/dhi0;
	for (i=0;i<sky->N;i++)
		sky->P[i].I=dhi0;
	dir=GHI-DHI;
	sky->sp=sun;
	sky->sI=dir/cos(sun.z);
}

double SolarAngle(double z1, double z2, double a1, double a2)
{
	double p=fmod(a1-a2,2*M_PI);
	double f;
	f=cos(z1)*cos(z2)+sin(z1)*sin(z2)*cos(p);
	if (f>1 && f<1.1)
		return 0.0;
	else if (f > 1.1 )
	{
		Warning("Warning: cannot compute gamma (angle between point and sun)");
		return 0.0;
	}
	return acos(f);
}

double 	air_mass(sky_pos sun)
{
	if (sun.z>M_PI/2)
		sun.z=M_PI/2;
	// following Kasten and Young (1989) airmass expression
	// 1/(cos(z)+0.50572*(96.07995-deg(z))^-1.6364
	return 1/(cos(sun.z)+0.50572*pow(96.07995-sun.z*180/M_PI,-1.6364));
}

#define	solconst 1367    /* solar constant W/m^2 */
double ExtraSolPower(double dayofyear)
{
	double day_angle;
	double E0;

	day_angle=2*M_PI*(dayofyear-1)/365;
	E0=1.00011+0.034221*cos(day_angle)+0.00128*sin(day_angle)+0.000719*cos(2*day_angle)+0.000077*sin(2*day_angle);
	return E0*solconst;
}


/* Perez sky's brightness */
double sky_brightness(sky_pos sun, double DHI, double dayofyear)
{
	double delta;	
	delta=DHI*air_mass(sun)/ExtraSolPower(dayofyear);
	if (delta<0.01)
		delta=0.01;
	return delta;
}

/* Perez sky's clearness */
#define kappa 1.041
double sky_clearness(sky_pos sun, double DHI, double GHI)
{
	double eps;	
	eps=((DHI+(GHI-DHI)/cos(sun.z))/DHI+kappa*sun.z*sun.z*sun.z)/(1+kappa*sun.z*sun.z*sun.z);
	return eps;
}

int ClearBin(double eps)
{
	if (eps < 1.065)
		return 0;
	else if (eps < 1.230) 
		return 1;
	else if (eps < 1.500)
		return 2;
	else if (eps < 1.950)
		return 3;
	else if (eps < 2.800)
		return 4;
	else if (eps < 4.500)
		return 5;
	else if (eps < 6.200)
		return 6;
	return 7;
}
/* Perez sky model coefficients */
/* Reference:	Perez, R., R. Seals, and J. Michalsky, 1993. "All- */
/*				Weather Model for Sky Luminance Distribution - */
/*				Preliminary Configuration and Validation," Solar */
/*				Energy 50(3):235-245, Table 1. */

static const double PerezCoeff[8][5][4] =
{
	{   /* clearness bin 0 epsilon: < 1.065 */
		{1.352500e+00, -2.576000e-01, -2.690000e-01, -1.436600e+00},
		{-7.670000e-01, 7.000000e-04, 1.273400e+00, -1.233000e-01},
		{2.800000e+00, 6.004000e-01, 1.237500e+00, 1.000000e+00},
		{1.873400e+00, 6.297000e-01, 9.738000e-01, 2.809000e-01},
		{3.560000e-02, -1.246000e-01, -5.718000e-01, 9.938000e-01},
	},
	{   /* clearness bin 1 epsilon: 1.065 to 1.230 */
		{-1.221900e+00, -7.730000e-01, 1.414800e+00, 1.101600e+00},
		{-2.054000e-01, 3.670000e-02, -3.912800e+00, 9.156000e-01},
		{6.975000e+00, 1.774000e-01, 6.447700e+00, -1.239000e-01},
		{-1.579800e+00, -5.081000e-01, -1.781200e+00, 1.080000e-01},
		{2.624000e-01, 6.720000e-02, -2.190000e-01, -4.285000e-01},
	},
	{   /* clearness bin 2 epsilon: 1.230 to 1.500 */
		{-1.100000e+00, -2.515000e-01, 8.952000e-01, 1.560000e-02},
		{2.782000e-01, -1.812000e-01, -4.500000e+00, 1.176600e+00},
		{2.472190e+01, -1.308120e+01, -3.770000e+01, 3.484380e+01},
		{-5.000000e+00, 1.521800e+00, 3.922900e+00, -2.620400e+00},
		{-1.560000e-02, 1.597000e-01, 4.199000e-01, -5.562000e-01},
	},
	{   /* clearness bin 3 epsilon: 1.500 to 1.950 */
		{-5.484000e-01, -6.654000e-01, -2.672000e-01, 7.117000e-01},
		{7.234000e-01, -6.219000e-01, -5.681200e+00, 2.629700e+00},
		{3.333890e+01, -1.830000e+01, -6.225000e+01, 5.207810e+01},
		{-3.500000e+00, 1.600000e-03, 1.147700e+00, 1.062000e-01},
		{4.659000e-01, -3.296000e-01, -8.760000e-02, -3.290000e-02},
	},
	{   /* clearness bin 4 epsilon: 1.950 to 2.800 */
		{-6.000000e-01, -3.566000e-01, -2.500000e+00, 2.325000e+00},
		{2.937000e-01, 4.960000e-02, -5.681200e+00, 1.841500e+00},
		{2.100000e+01, -4.765600e+00, -2.159060e+01, 7.249200e+00},
		{-3.500000e+00, -1.554000e-01, 1.406200e+00, 3.988000e-01},
		{3.200000e-03, 7.660000e-02, -6.560000e-02, -1.294000e-01},
	},
	{   /* clearness bin 5 epsilon: 2.800 to 4.500 */
		{-1.015600e+00, -3.670000e-01, 1.007800e+00, 1.405100e+00},
		{2.875000e-01, -5.328000e-01, -3.850000e+00, 3.375000e+00},
		{1.400000e+01, -9.999000e-01, -7.140600e+00, 7.546900e+00},
		{-3.400000e+00, -1.078000e-01, -1.075000e+00, 1.570200e+00},
		{-6.720000e-02, 4.016000e-01, 3.017000e-01, -4.844000e-01},
	},
	{   /* clearness bin 6 epsilon: 4.500 to 6.200 */
		{-1.000000e+00, 2.110000e-02, 5.025000e-01, -5.119000e-01},
		{-3.000000e-01, 1.922000e-01, 7.023000e-01, -1.631700e+00},
		{1.900000e+01, -5.000000e+00, 1.243800e+00, -1.909400e+00},
		{-4.000000e+00, 2.500000e-02, 3.844000e-01, 2.656000e-01},
		{1.046800e+00, -3.788000e-01, -2.451700e+00, 1.465600e+00},
	},
	{   /* clearness bin 7 epsilon: > 6.200 */
		{-1.050000e+00, 2.890000e-02, 4.260000e-01, 3.590000e-01},
		{-3.250000e-01, 1.156000e-01, 7.781000e-01, 2.500000e-03},
		{3.106250e+01, -1.450000e+01, -4.611480e+01, 5.537500e+01},
		{-7.231200e+00, 4.050000e-01, 1.335000e+01, 6.234000e-01},
		{1.500000e+00, -6.426000e-01, 1.856400e+00, 5.636000e-01},
	}	
};

#define _X PerezCoeff[clearbin]
void ParamPerez(double eps, double delta, sky_pos sun, double *a, double *b, double *c, double *d, double *e)
{
	int clearbin;
	
	clearbin=ClearBin(eps);
	if ((eps > 1.065)&&(eps < 2.8))
		if (delta < 0.2)
			delta = 0.2;

	if (clearbin != 0)
	{
		/* Calculate parameter a, b, c, d and e (Eqn. 6) */
		(*a)=_X[0][0]+_X[0][1]*sun.z+delta*(_X[0][2]+_X[0][3]*sun.z);
		(*b)=_X[1][0]+_X[1][1]*sun.z+delta*(_X[1][2]+_X[1][3]*sun.z);
		(*c)=_X[2][0]+_X[2][1]*sun.z+delta*(_X[2][2]+_X[2][3]*sun.z);
		(*d)=_X[3][0]+_X[3][1]*sun.z+delta*(_X[3][2]+_X[3][3]*sun.z);
		(*e)=_X[4][0]+_X[4][1]*sun.z+delta*(_X[4][2]+_X[4][3]*sun.z);
	}
	else
	{
		/* Parameters a, b and e (Eqn. 6) */
		(*a)=_X[0][0]+_X[0][1]*sun.z+delta*(_X[0][2]+_X[0][3]*sun.z);
		(*b)=_X[1][0]+_X[1][1]*sun.z+delta*(_X[1][2]+_X[1][3]*sun.z);
		(*e)=_X[4][0]+_X[4][1]*sun.z+delta*(_X[4][2]+_X[4][3]*sun.z);

		/* Parameter c (Eqn. 7) */
		(*c)=exp(pow(delta*(_X[2][0]+_X[2][1]*sun.z),_X[2][2]))-_X[2][3];
		/* Parameter d (Eqn. 8) */
		(*d)=-exp(delta*(_X[3][0]+_X[3][1]*sun.z))+_X[3][2]+delta*_X[3][3];
	}
}
#undef _X

double Fperez(double z,double g, double a, double b, double c, double d, double e)
{
	return (1+a*exp(b/cos(z)))*(1+c*exp(d*g)+e*cos(g)* cos(g));
}

void PerezSky(sky_grid * sky, sky_pos sun, double GHI, double DHI, double dayofyear)
{
	double dhi0=0, dir, g;
	int i;
	double a, b, c, d, e;
	double eps, delta;
	eps=sky_clearness(sun, DHI, GHI);
	delta=sky_brightness(sun, DHI, dayofyear);
	ParamPerez(eps, delta, sun, &a, &b, &c, &d, &e);
	Print(VERBOSE,"Creating Perez all weather sky-dome\n");	
	Print(VERBOSE,"GHI:           %e\nDHI:           %e\nDay:           %f\n", GHI, DHI, dayofyear);
	Print(VERBOSE,"Solar Zenith:  %f\nSolar Azimuth: %f\n", rad2degr(sun.z), rad2degr(sun.a));
	Print(VERBOSE,"******** Derived Parametrization\n");
	Print(VERBOSE,"Clearness:     %e\nDelta:         %e\n", eps, delta);	
	Print(VERBOSE,"a:             %e\nb:             %e\nc:             %e\nd:             %e\ne:             %e\n", a,b,c,d,e);
	Print(VERBOSE,"********\n");
	
	// sum intensity cos(z) product to normalize intensities
	for (i=0;i<sky->N;i++)
	{
		g=SolarAngle(sky->P[i].p.z,sun.z, sky->P[i].p.a, sun.a);   
		sky->P[i].I=Fperez(sky->P[i].p.z,g, a, b, c, d, e);                                                   
		dhi0+=sky->P[i].I*cos(sky->P[i].p.z);
	}
	dhi0=DHI/dhi0;// correction factor
	for (i=0;i<sky->N;i++)
		sky->P[i].I*=dhi0;
	dir=GHI-DHI; // direct contribution
	sky->sp=sun;
	sky->sI=dir/cos(sun.z);
}
