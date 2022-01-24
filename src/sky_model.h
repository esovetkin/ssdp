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
#ifndef _SKY_MODEL_H
#define _SKY_MODEL_H
typedef struct {
	char gradation, indicatrix;
	double a,b,c,d,e;
	char *descr;
} CIE_SKY;

// table1 from CIE S 011/E:2003 (ISO 15469:2004(E))
extern CIE_SKY CIE_SKIES[];
// sky type defines
typedef enum CIE_SKY_TYPE {
	CIE_SKY_TYPE_1,
	CIE_SKY_TYPE_2,
	CIE_SKY_TYPE_3,
	CIE_SKY_TYPE_4,
	CIE_SKY_TYPE_5,
	CIE_SKY_TYPE_6,
	CIE_SKY_TYPE_7,
	CIE_SKY_TYPE_8,
	CIE_SKY_TYPE_9,
	CIE_SKY_TYPE_10,
	CIE_SKY_TYPE_11,
	CIE_SKY_TYPE_12,
	CIE_SKY_TYPE_13,
	CIE_SKY_TYPE_14,
	CIE_SKY_TYPE_15
} CIE_SKY_TYPE;

void UniformSky(sky_grid *sky, sky_pos sun, double GHI, double DHI);
void SkySunOnly(sky_grid *sky, sky_pos sun, double GHI, double DHI);
void PerezSky(sky_grid * sky, sky_pos sun, double GHI, double DHI, double dayofyear);
void CumulativePerezSky(sky_grid * sky, sky_pos *sun, double *t, double *GHI, double *DHI, double *dayofyear, int N);
#endif
