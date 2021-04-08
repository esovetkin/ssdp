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


vec cross(vec a, vec b)
{
	vec r;
	r.x=a.y*b.z-a.z*b.y;
	r.y=a.z*b.x-a.x*b.z;
	r.z=a.x*b.y-a.y*b.x;
	return r;
}

vec sum(vec a, vec b)
{
	vec r;
	r.x=a.x+b.x;
	r.y=a.y+b.y;
	r.z=a.z+b.z;
	return r;
}
vec diff(vec a, vec b)
{
	vec r;
	r.x=a.x-b.x;
	r.y=a.y-b.y;
	r.z=a.z-b.z;
	return r;
}

vec scalevec(vec a, double b)
{
	vec r;
	r.x=a.x*b;
	r.y=a.y*b;
	r.z=a.z*b;
	return r;
}

double dot(vec a, vec b)
{
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

double norm(vec v)
{
	return sqrt(dot(v,v));
}

vec unit(sky_pos a)
{
	vec r;
	r.x=sin(a.z)*cos(a.a);
	r.y=sin(a.z)*sin(a.a);
	r.z=cos(a.z);
	return r;
}
sky_pos vecdir(vec a)
{
	sky_pos r;
	r.a=atan2(a.y,a.x);
	r.z=atan2(sqrt(a.x*a.x+a.y*a.y), a.z);
	return r;
}
// projection routines
sky_pos rrf(sky_pos p, sky_pos axis, double beta)
{
	// Rodrigues' Rotation Formula
	vec a,k,b;
	k=unit(axis);
	a=unit(p);
	b=scalevec(a, cos(beta));
	b=sum(b,scalevec(cross(k,a), sin(beta)));
	b=sum(b,scalevec(k, dot(k,a)*(1-cos(beta))));
	return vecdir(b);	
}

double amean(double a1, double a2)
{
	vec v;
	sky_pos p;
	v.x=(cos(a1)+cos(a2))/2;
	v.y=(sin(a1)+sin(a2))/2;
	v.z=0;
	p=vecdir(v);
	return p.a;
}

double adiff(double a1, double a2)
{
	double d;
	a1=fmod(a1,2*M_PI);
	a2=fmod(a2,2*M_PI);
	
	d=a1-a2;
	if (fabs(d)>M_PI)
	{
		if (d>0)
			return d-2*M_PI;
		else
			return d+2*M_PI;
	}
	return d;
}


double rad2degr(double rad)
{
	return rad*180/M_PI;
}
double degr2rad(double degr)
{
	return degr*M_PI/180;
}
