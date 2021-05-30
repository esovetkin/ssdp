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
//BEGIN_SSDP_EXPORT
typedef struct topology {
	double *x, *y, *z;  // 3D coordinates
	int N;		 		// number of points
	triangles *T;
	int Nt;		 		// number of triangles
	nodetree *P;
} topology;
typedef struct topogrid {
	double *z;  // z coordinate in column-mayor format
	int *sort;  // sorted indexing with increasing height
	float *A1, *A2; // pre computed discrete angles
	int Na;		// number of discrete angles in A1 and A2 arrays
	int Nx,Ny;	// number of points
	double x1, y1; // lower left corner
	double x2, y2; // upper right corner
} topogrid;
//END_SSDP_EXPORT


topology MakeTopology(double *x, double *y, double *z, int N);
topogrid MakeTopogrid(double *z, double x1, double y1, double x2, double y2, int Nx, int Ny);

topology CreateRandomTopology(double dx, double dy, double dz, int N1, int N2);
void free_topo (topology *T);
void free_topogrid (topogrid *T);
double SampleTopo(double x, double y, topology *T, sky_pos *sn);
double SampleTopoGrid(double x, double y, topogrid *T, sky_pos *sn);

void ComputeHorizon(horizon *H, topology *T, double minzen, double xoff, double yoff, double zoff);
void ComputeGridHorizon(horizon *H, topogrid *T, double minzen, double xoff, double yoff, double zoff);

horizon MakeHorizon(sky_grid *sky, topology *T, double xoff, double yoff, double zoff);
