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
typedef struct topology {
	double *x, *y, *z;  // 3D coordinates
	int N;		 		// number of points
	triangles *T;
	int Nt;		 		// number of triangles
	nodetree *P;
} topology;

topology MakeTopology(double *x, double *y, double *z, int N);
topology CreateRandomTopology(double dx, double dy, double dz, double fN, int N);
void free_topo (topology *T);
void MakeHorizon(sky_grid *sky, topology T, double xoff, double yoff, double zoff, vec *sn);
void ClearHorizon(sky_grid *sky);
