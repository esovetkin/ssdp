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
#ifndef _VECTOR_H
#define _VECTOR_H

#include <hdf5.h>

// basic 3D vector calculation
//BEGIN_SSDP_EXPORT
typedef struct sky_pos {
	double z,a;
} sky_pos;
//END_SSDP_EXPORT

typedef struct vec {
	double x,y,z;
} vec;
vec cross(vec a, vec b);
vec sum(vec a, vec b);
vec diff(vec a, vec b);
vec scalevec(vec a, double b);
double dot(vec a, vec b);
double norm(vec v);
vec unit(sky_pos a);
sky_pos vecdir(vec a);
double amean(double a1, double a2);
double adiff(double a1, double a2);
sky_pos rrf(sky_pos p, sky_pos axis, double beta);
#define RADPDEG 1.745329251994329e-02
#define DEGPRAD 5.729577951308232e+01
#define rad2deg(r) (DEGPRAD*(r))
#define deg2rad(d) (RADPDEG*(d))

hid_t h5t_sky_pos();

#endif
