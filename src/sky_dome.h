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
#ifndef _SKY_DOME_H
#define _SKY_DOME_H

//BEGIN_SSDP_EXPORT
typedef struct hexpatch {
	double I;
	sky_pos p;
	int NL[7]; // next level neighbors
	int PL[3]; // previous level neighbors
	int NI;    // iso level next
	int PI;    // iso level previous
} hexpatch;

typedef struct sky_grid {
	hexpatch *P;
	// get sun out of the hexpatch dome so we can separate contributions
	// if we want
	sky_pos sp;	// solar position
	double sI;	// solar intensity
	int suni; // index of the sky patch the sun is in
	double *cosz; // per patch cosine of zenith
	double *sa; // per patch solid angle
	double icosz; // integral of all cosz and solid angle factors, used for normalization
	int N;		
	int Nz;
} sky_grid;
//END_SSDP_EXPORT


void Connectivity(int Nz);
sky_grid InitSky(int Nz);
void free_sky_grid(sky_grid *sky);
int FindPatch(sky_grid *sky, sky_pos p);

int NNZ(int n);	// number of elements in a mesh given a number of levels (zenith discretizations)
int NZN(int n);	// inverse of above, i.e. give it a index and it returns the level (zenith) index
/* sun solid angle:
 * fe '(appi(64)*(1391400.0/2)^2)/(149.597870700e6^2)'
 * diameter sun:   1391400    km
 * distance sun: 149597870.7  km
 * solid angle : ~ A/d^2 = (pi * 1391400 ^2 / 4) / 149597870.7^2 = 6.7942739713694071218e-05
 * Note that the actual distance sun-earth varies about 3% so the result is not to be taken as accurate
 */
#define SUNSR 6.794e-5	
#endif
