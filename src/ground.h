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


enum SampleType {
/** Approximate horizon sampling strategy type
 */
		PRECISE = 0, // precise horizon algorithm. Walk through all
				     // available topography pixels
		SOBOL   = 1, // Sobol low-discrepancy sequence (quasi MC)
		IID     = 2, // draw iid random varaibles (MC)
		/** use fixed number of rays and measure horizon along the
		 * centreline
		 */
		RAYS16  = 16,
		RAYS32  = 32,
		RAYS64  = 64,
		RAYS128 = 128,
};


struct hsample_data {
		int x, y;
		double d;
};


typedef struct topogrid {
		double *z;  // z coordinate in column-major format
		int *sort;  // sorted indexing with increasing height
		double *A1, *A2; // pre computed discrete angles
		int Na;		// number of discrete angles in A1 and A2 arrays
		int Nx,Ny;	// number of points
		double x1, y1; // lower left corner
		double x2, y2; // upper right corner
		double dx, dy; // dx and dy step size
		enum SampleType horizon_stype; // type of sampling
		struct hsample_data* horizon_sample; // locations where topography for horizon is sampled
		int horizon_nsample; // number of points horizon sample points
		int horizon_nsample_eff; // effective number of horizon sample points
		double *horizon_dstr; // sampling distribution (ordered quantiles)
		int horizon_dstrn; // len(horizon_dstr)
		double *horizon_phid; // azimuth sample distribution
		int horizon_nphid; // len(horizon_phid)
} topogrid;
//END_SSDP_EXPORT


topology MakeTopology(double *x, double *y, double *z, int N);
topogrid MakeTopogrid(double *z, double x1, double y1, double x2, double y2, int Nx, int Ny);
topogrid MakeTopoGDAL(double x1, double y1, double x2, double y2, char **fns, int nfns, double step, int epsg);

topology CreateRandomTopology(double dx, double dy, double dz, int N1, int N2);
void free_topo (topology *T);
void free_topogrid (topogrid *T);
double SampleTopo(double x, double y, topology *T, sky_pos *sn);
double SampleTopoGrid(double x, double y, topogrid *T, sky_pos *sn);

int FillMissingTopoGrid(topogrid *T, double na, int maxwalk);
int AddHeightTopoGrid(topogrid *T, double *x, double *y, double *z, int nx, int nz);
int BlurTopoGrid(topogrid *T, int size);

int HorizonSet(topogrid *T, int nsample, enum SampleType type);
int HorizonSetDstr(topogrid *T, double *d, int nd);
int HorizonSetPhi(topogrid *T, double *phi, int nphi);

void ComputeHorizon(horizon *H, topology *T, double minzen, double xoff, double yoff, double zoff);
void ComputeGridHorizon(horizon *H, topogrid *T, double minzen, double xoff, double yoff, double zoff);

horizon MakeHorizon(sky_grid *sky, topology *T, double xoff, double yoff, double zoff);
