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
#ifndef _RROJECT_H
#define _RROJECT_H
//BEGIN_SSDP_EXPORT
typedef enum {AOI_NONE, AOI_GLASS, AOI_GLASS_AR, AOI_USER} AOI_Model;
typedef struct {
	double ng, nar;
	double *theta, *effT;
	int N;
	AOI_Model M;
} AOI_Model_Data;
//END_SSDP_EXPORT

AOI_Model_Data InitAOIModel(AOI_Model model, double ng, double nar , double *theta, double *effT, int N);
void POA_Sky_Transfer(sky_grid *sky, sky_transfer *T, sky_pos pn, AOI_Model_Data *M);
double POA_Albedo_Transfer(sky_grid *sky, sky_pos pn, AOI_Model_Data *M);
double DiffusePlaneOfArray(sky_grid *sky, sky_transfer *T);
double DirectPlaneOfArray(sky_grid *sky, horizon *H, sky_pos pn, AOI_Model_Data *M);
void POA_to_SurfaceNormal(sky_pos *pn, sky_pos sn);
#endif
