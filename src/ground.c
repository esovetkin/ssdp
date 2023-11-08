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
#include <time.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "vector.h"
#include "sky_dome.h"
#include "trace.h"
#include "delaunay.h"
#include "ground.h"
#include "print.h"
#include "print.h"
#include "error.h"
#include "config.h"
#include "fatan2.h"
#include "epsg.h"
#include "iset.h"
#include "edt.h"
#include "minfill.h"
#include "topogdal.h"
#include "sobolseq.h"
#include "filterimage/filter.h"

#ifdef RUNMEMTEST
#include "random_fail_malloc.h"
#define malloc(x) random_fail_malloc(x)
#endif

// sort triangle list by height
int tcomp(const void *a, const void *b)
{ 
	if (((triangles *)a)->ccz<((triangles *)b)->ccz)
		return -1;
	if (((triangles *)a)->ccz>((triangles *)b)->ccz)
		return 1;
	return 0;
} 
void tsort(triangles *T, int Nt)
{ 
	qsort(T, Nt, sizeof(triangles), tcomp);
}

topology MakeTopology(double *x, double *y, double *z, int N)
{
	topology T;
	topology T0={NULL,NULL,NULL,0,NULL,0,NULL};
	int i;
	// ERRORFLAG MALLOCFAILTOPOLOGY  "Error memory allocation failed in creating a topology"
	T.x=malloc(N*sizeof(double));
	T.y=malloc(N*sizeof(double));
	T.z=malloc(N*sizeof(double));
	if ((T.x==NULL)||(T.y==NULL)||(T.z==NULL))
	{
		AddErr(MALLOCFAILTOPOLOGY);
		return T0;
	}
	for (i=0;i<N;i++)
	{
		T.x[i]=x[i];
		T.y[i]=y[i];
		T.z[i]=z[i];
	}
	T.N=N;
	T.T=Triangulate(T.x, T.y, T.N, &T.Nt);
	if (ssdp_error_state)
	{
		free_topo(&T);
		return T0;
	}
	for (i=0;i<T.Nt;i++)
		T.T[i].ccz=(T.z[T.T[i].i]+T.z[T.T[i].j]+T.z[T.T[i].k])/3;
	tsort(T.T, T.Nt); // sorting triangles
	
	T.P=InitTree(T.T, T.Nt, T.x, T.y, N);
	if (ssdp_error_state)
	{
		free_topo(&T);
		return T0;
	}
	return T;
}

// sort grid index by height
// ERRORFLAG MALLOCFAILTOPOGRID  "Error memory allocation failed in creating a topogrid"
#ifdef _WIN32 // WIN32 misses the qsort_r function, have to use a global variable
double *_SortHeight;
int gcomp(const void *a, const void *b)
{
	if (((double *)_SortHeight)[((int *)a)[0]]<((double *)_SortHeight)[((int *)b)[0]])
		return -1;
	if (((double *)_SortHeight)[((int *)a)[0]]>((double *)_SortHeight)[((int *)b)[0]])
		return 1;
	return 0;
} 
int *gsort(double *z, int N)
{ 
	int i;
	int *sort;
	sort=malloc(N*sizeof(int));
	if (sort==NULL)
	{
		AddErr(MALLOCFAILTOPOGRID);
		return NULL;
	}
	_SortHeight=z;
	
	for (i=0;i<N;i++)
		sort[i]=i;
	qsort(sort, N, sizeof(int), gcomp);
	return sort;
}
#else
int gcomp(const void *a, const void *b, void *c)
{
	if (((double *)c)[((int *)a)[0]]<((double *)c)[((int *)b)[0]]) // thats a lot of brackets...
		return -1;
	if (((double *)c)[((int *)a)[0]]>((double *)c)[((int *)b)[0]])
		return 1;
	return 0;
} 
int *gsort(double *z, int N)
{ 
	int i;
	int *sort;
	sort=malloc(N*sizeof(int));
	if (sort==NULL)
	{
		AddErr(MALLOCFAILTOPOGRID);
		return NULL;
	}
	for (i=0;i<N;i++)
		sort[i]=i;
	qsort_r(sort, N, sizeof(int), gcomp, z);
	return sort;
}
#endif /*_WIN32 */

// describes discrete angular ranges for elements in a regular grid
// only computes fro 1/8th of the elements, the reast may be inferred
void MakeAngles(int n, double dx, double dy, double **A1, double **A2, int *N)
{
	int i, j, k;
	double d, w;
	(*N)=(4*n*n+20*n)/8+2;
	(*A1)=malloc((*N)*sizeof(double));
	(*A2)=malloc((*N)*sizeof(double));
	if (((*A1)==NULL)||((*A2)==NULL))
	{
		AddErr(MALLOCFAILTOPOGRID);
		return;
	}
	for (i=1;i<n;i++)
		for (j=0;j<=i;j++)
		{
			k=i-1;
			k=(4*k*k+20*k)/8+j;
			(*A1)[k]=M_PI/2-ATAN((dy*((double)j))/(dx*((double)i))); // base angle, swap x and y for the angle
			d=sqrt(dx*dx*((double)(i*i))+dy*dy*((double)(j*j)));
			w=ATAN(sqrt(dx*dx+dy*dy)/d)/2; // width measured against the diagonal of the element
			(*A2)[k]=(*A1)[k]+w;
			(*A1)[k]-=w;
		}
}

int Arange(int dx, int dy, double *a1, double *a2, double *A1, double *A2)
{
	int k;
	if ((dx==0)&&(dy==0))
	{
		(*a1)=0;
		(*a2)=2*M_PI;
		return 0;
	}
	if ((dx>=0)&&(dy>0)) // note that north (x=0, y>0) is 0 rad and east (x=1, y=0) is pi/2 rad! 
	{
		if (dx>dy)
		{
			// second 1/8th
			k=dx-1;
			k=(4*k*k+20*k)/8+dy;
			(*a1)=A1[k];
			(*a2)=A2[k];
			return 1;
		}
		// first  1/8th
		k=dy-1;
		k=(4*k*k+20*k)/8+dx;
		(*a1)=M_PI/2-A2[k];
		(*a2)=M_PI/2-A1[k];
		return 1;
	}
	if ((dx<0)&&(dy>=0))
	{
		dx=abs(dx);
		if (dx>dy)
		{
			// seventh 1/8th
			k=dx-1;
			k=(4*k*k+20*k)/8+dy;
			
			(*a1)=2*M_PI-A2[k];
			(*a2)=2*M_PI-A1[k];
			return 1;
		}
		// eighth  1/8th
		k=dy-1;
		k=(4*k*k+20*k)/8+dx;
		(*a1)=3*M_PI/2+A1[k];
		(*a2)=3*M_PI/2+A2[k];
		return 1;
	}
	if ((dx<=0)&&(dy<0))
	{
		dx=abs(dx);
		dy=abs(dy);
		if (dx>dy)
		{
			// sixth 1/8th
			k=dx-1;
			k=(4*k*k+20*k)/8+dy;
			(*a1)=M_PI+A1[k];
			(*a2)=M_PI+A2[k];
			return 1;
		}
		// fifth  1/8th
		k=dy-1;
		k=(4*k*k+20*k)/8+dx;
		(*a1)=3*M_PI/2-A2[k];
		(*a2)=3*M_PI/2-A1[k];
		return 1;
	}
	if ((dx>0)&&(dy<=0))
	{
		dy=abs(dy);
		if (dx>dy)
		{
			// third 1/8th
			k=dx-1;
			k=(4*k*k+20*k)/8+dy;
		
			(*a1)=M_PI-A2[k];
			(*a2)=M_PI-A1[k];
			return 1;
		}
		// fourth  1/8th
		k=dy-1;
		k=(4*k*k+20*k)/8+dx;
		(*a1)=M_PI/2+A1[k];
		(*a2)=M_PI/2+A2[k];
		return 1;
	}
	(*a1)=0;
	(*a2)=0;
	return 0;
}


topogrid MakeTopogrid(double *z, double x1, double y1, double x2, double y2, int Nx, int Ny)
{
	topogrid T;
	topogrid T0={NULL,NULL,NULL,NULL,0,0,0,0,0,1,1,-1,1,NULL,NULL};
	int i, N;
	double dx, dy;
	N=Nx*Ny;
	T.z=malloc(N*sizeof(double));
	if (T.z==NULL)
	{
		AddErr(MALLOCFAILTOPOGRID);
		return T0;
	}
	for (i=0;i<N;i++)
		T.z[i]=z[i];
	T.sort=gsort(T.z, N);	
	if (ssdp_error_state)
	{
		free_topogrid(&T);
		return T0;
	}

	if (NULL==(T.horizon_idx=malloc(N*sizeof(*T.horizon_idx)))) {
		free_topogrid(&T);
		return T0;
	}
	for (i=0; i < N; ++i)
		T.horizon_idx[i] = (2*Ny-1)*(T.sort[i]/Ny) + T.sort[i]%Ny;

	dx=(x2-x1)/((double)Nx);
	dy=(y2-y1)/((double)Nx);
	if (Nx>Ny)
		MakeAngles(Nx, dx, dy, &(T.A1), &(T.A2), &(T.Na));
	else
		MakeAngles(Ny, dx, dy, &(T.A1), &(T.A2), &(T.Na));	
	if (ssdp_error_state)
	{
		free_topogrid(&T);
		return T0;
	}
	T.Nx=Nx;
	T.Ny=Ny;
	T.horizon_sample = NULL;
	T.horizon_nsample = -1;
	T.horizon_decay = (double) 1;
	T.x1=x1;
	T.y1=y1;
	T.x2=x2;
	T.y2=y2;
	return T;
}


topogrid MakeTopoGDAL(double x1, double y1, double x2, double y2, char **fns, int nfns, double step, int epsg)
{
        topogrid T={NULL,NULL,NULL,NULL,0,0,0,0,0,1,1,-1,1,NULL,NULL};
        struct coordinates *lcs = box2coordinates(x1,y1,x2,y2,step,epsg);
        if (NULL == lcs) goto maketopogdal_elcs;

        printf("Location box dimensions: %d, %d\n", lcs->nx, lcs->ny);

        struct gdaldata *gd = gdaldata_init((const char **)fns, nfns);
        if (NULL == gd) goto maketopogdal_egd;

        double *z;
        z = topogrid_from_gdal(gd, lcs);
        if (NULL == z) goto maketopogdal_ez;

        T = MakeTopogrid(z, lcs->y1, lcs->x1, lcs->y2, lcs->x2, lcs->ny, lcs->nx);

        free(z);
        gdaldata_free(gd);
        coordinates_free(lcs);
        return T;
maketopogdal_ez:
        gdaldata_free(gd);
maketopogdal_egd:
        coordinates_free(lcs);
maketopogdal_elcs:
        AddErr(MALLOCFAILTOPOGRID);
        return T;
}


void free_topo (topology *T)
{
	if (T)
	{
		if (T->x)
			free(T->x);
		if (T->y)
			free(T->y);
		if (T->z)
			free(T->z);
		if (T->T)
			free(T->T);
		T->x=NULL;
		T->y=NULL;
		T->z=NULL;	
		T->T=NULL;	
		T->N=0;	
		if (T->P)
		{
			FreeTree(T->P);
			free(T->P);
			T->P=NULL;
		}
	}
}
void free_topogrid (topogrid *T)
{
	if (T)
	{
		if (T->z)
		{
			free(T->z);
			T->z=NULL;
		}
		if (T->A1)
		{
			free(T->A1);
			T->A1=NULL;
		}
		if (T->A2)
		{
			free(T->A2);
			T->A2=NULL;
		}
		if (T->sort)
		{
			free(T->sort);
			T->sort=NULL;
		}
		T->Nx=0;
		T->Ny=0;
		T->x1=0;
		T->y1=0;
		T->x2=1;
		T->y1=1;
		T->horizon_nsample = -1;
		T->horizon_decay = (double) 1;
		free(T->horizon_sample); T->horizon_sample=NULL;
		free(T->horizon_idx); T->horizon_idx=NULL;
	}
}

// should export norm somehow as we may want to rotate the pv panel accordingly
double SampleTopo(double x, double y, topology *T, sky_pos *sn)
{ // dangerous, no check for extrapolation
	int n;
	vec nn, a, b, c, v1, v2;
	double D,z, l;
	n=treesearch(T->T, T->P, T->x, T->y, T->Nt, x, y);
	a.x=T->x[T->T[n].i];
	a.y=T->y[T->T[n].i];
	a.z=T->z[T->T[n].i];
	b.x=T->x[T->T[n].j];
	b.y=T->y[T->T[n].j];
	b.z=T->z[T->T[n].j];
	c.x=T->x[T->T[n].k];
	c.y=T->y[T->T[n].k];
	c.z=T->z[T->T[n].k];
	v1=diff(a,b);
	v2=diff(a,c);
	nn=cross(v2,v1);
	l=norm(nn);
	
	/* the delaunay code should have arranged the vertices such that we do not need this check
	 if (nn.z<0)
		nn=scalevec(nn,-1);
	*/
	if (nn.z<1e-10*l) // avoid div 0
		nn.z=1e-10*l;
	//A x + B y +C z = D
	D=nn.x*a.x+nn.y*a.y+nn.z*a.z;
	z=(D-nn.x*x-nn.y*y)/nn.z;
	if (sn)
	{
		a=nn;
		// in our coordinate system we go counter clockwise with 0 degrees pointing up (swap x and y)
		nn.x=nn.y;
		nn.y=a.x;
		(*sn)=vecdir(nn);
	}
	return z;	
}

static inline int YINDEX(int N, int Ny)
{
	return N%Ny;
} 
static inline int XINDEX(int N, int Ny)
{
	return N/Ny;
}
#ifdef INDEX
#undef INDEX
#endif
static inline int INDEX(int x, int y, int Ny)
{
	return x*Ny+y;
} 

static inline int IndexGridX(double x, topogrid *T, int *o)
{
	double dx=(T->x2-T->x1)/((double)T->Nx);
	double nx;
	nx=round((x-T->x1)/dx);
	if (o)
	{
		if (nx*dx>x)
			(*o)=(int)nx-1;
		else
			(*o)=(int)nx+1;
	}
	return (int)nx;
}

static inline int IndexGridY(double y, topogrid *T, int *o)
{
	double dy=(T->y2-T->y1)/((double)T->Ny);
	double ny;
	ny=round((y-T->y1)/dy);
	if (o)
	{
		if (ny*dy>y)
			(*o)=(int)ny -1;
		else
			(*o)=(int)ny +1;
	}
	return (int)ny;
}

static inline double YVALUE(int i,topogrid *T)
{
	double dy=(T->y2-T->y1)/((double)T->Ny);
	return T->y1+((double)i)*dy;
} 
static inline double XVALUE(int i,topogrid *T)
{
	double dx=(T->x2-T->x1)/((double)T->Nx);
	return T->x1+((double)i)*dx;
} 

// should export norm somehow as we may want to rotate the pv panel accordingly
double SampleTopoGrid(double x, double y, topogrid *T, sky_pos *sn)
{ // dangerous, no check for extrapolation
	int i, j, k, l;
	vec nn, a, b, c, v1, v2;
	double D,z, len;
	i=IndexGridX(x, T, &k);
	j=IndexGridY(y, T, &l);
	if (i<=0)
	{
		i=0;
		k=1;
	}
	if (i>=T->Nx-1)
	{
		i=T->Nx-1;
		k=i-1;// look backward for a triangle
	}
		
	if (j<=0)
	{
		j=0;
		l=1;
	}
	if (j>=T->Ny-1)
	{
		j=T->Ny-1;
		l=j-1;// look backward for a triangle
	}
	
	a.x=XVALUE(i,T);
	a.y=YVALUE(j,T);
	a.z=T->z[INDEX(i, j, T->Ny)];
	b.x=XVALUE(k,T);
	b.y=YVALUE(j,T);
	b.z=T->z[INDEX(k, j, T->Ny)];
	c.x=XVALUE(i,T);
	c.y=YVALUE(l,T);
	c.z=T->z[INDEX(i, l, T->Ny)];
	
	v1=diff(a,b);
	v2=diff(a,c);
	nn=cross(v2,v1);
	len=norm(nn);
	
	if (nn.z<0)
		nn=scalevec(nn,-1);
	
	if (nn.z<1e-10*len) // avoid div 0
		nn.z=1e-10*len;
	//A x + B y +C z = D
	D=nn.x*a.x+nn.y*a.y+nn.z*a.z;
	z=(D-nn.x*x-nn.y*y)/nn.z;
	if (sn)
	{
		a=nn;
		// in our coordinate system we go counter clockwise with 0 degrees pointing up (swap x and y)
		nn.x=nn.y;
		nn.y=a.x;
		(*sn)=vecdir(nn);
	}
	return z;	
}


static int _edt_fillmissing(topogrid *T, double na)
{
		struct edt *dt;
		if (NULL==(dt=edt_init(T->z, T->Nx*T->Ny, T->Ny, T->Nx, na)))
				goto edt;

		edt_compute(dt);
		edt_fill(dt);
		edt_free(dt);
		return 0;
edt:
		return -1;
}


static int _minmax_fillmissing(topogrid *T, double na, int maxwalk)
{
		struct minfill *mf = minfill_init(T->z, T->Ny, T->Nx, na);
		if (NULL == mf) goto emf;

		mf->maxwalk = abs(maxwalk);
		mf->direction = (maxwalk > 0) - (maxwalk < 0);
		if (minfill_fill(mf)) goto efill;
		minfill_free(mf);

		return 0;
efill:
		minfill_free(mf);
emf:
		return -1;
}


int FillMissingTopoGrid(topogrid *T, double na, int maxwalk)
{
		if (0 == maxwalk)
				return _edt_fillmissing(T, na);

		return _minmax_fillmissing(T, na, maxwalk);
}


int AddHeightTopoGrid(topogrid *T, double *x, double *y, double *z, int n, int nz)
{
		int i, ix, iy, iz, t;
		struct iset *set;
		if (NULL==(set=iset_init(2*n))) goto eset;

		for (i=0; i < n; ++i) {
				ix = IndexGridX(x[i], T, &t);
				iy = IndexGridY(y[i], T, &t);

				// do not alter height at borders
				if (ix<=0 || ix>=T->Nx-1 ||
					iy<=0 || iy>=T->Ny-1)
						continue;

				iz = INDEX(ix, iy, T->Ny);

				// do not alter height at already altered pixels
				if (iset_isin(set, (unsigned int) iz)) continue;

				T->z[iz] += z[i % nz];
				if (iset_insert(set, (unsigned int) iz)) goto einsert;
		}

		iset_free(set);
		return 0;
einsert:
		iset_free(set);
eset:
		return -1;
}


int BlurTopoGrid(topogrid *T, int size)
{
		double f[1] = {1};
		int dmx[1] = {0}, dmy[1] = {0};
		filterset F = DerivOperatorSet2D(size,size,size,size,0,dmx,dmy,f,1,'p');
		if (ssdp_error_state) goto efilter;

		image Iin, Iout;
		Iin.I = T->z;
		Iin.N = T->Ny;
		Iin.M = T->Nx;
		Iout = ApplyFilter(Iin, 1, 1, F);
		if (ssdp_error_state) goto efilter;

		int i;
		for (i=0; i < T->Nx*T->Ny; ++i)
				T->z[i] = Iout.I[i];

		FreeImage(&Iout);
		FreeFilterSet(&F);
		return 0;
efilter:
		return -1;
		// errors?
}


topology CreateRandomTopology(double dx, double dy, double dz, int N1, int N2)
{
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double *x, *y, *z;
	topology T;
	topology T0={NULL,NULL,NULL,0,NULL,0,NULL};
	int i, n;
	if (N1>N2)
	{
		n=N1;
		N1=N2;
		N2=n;
	}
	if (N1<3)
	{
		N1=3;
		N2+=3;
	}
	// ERRORFLAG MALLOCFAILRANDTOPOLOGY  "Error memory allocation failed in generating a random topology"
	x=malloc(N2*sizeof(double));
	y=malloc(N2*sizeof(double));
	z=malloc(N2*sizeof(double));
	if ((x==NULL)||(y==NULL)||(z==NULL))
	{
		AddErr(MALLOCFAILRANDTOPOLOGY);
		return T0;
	}
	srand(time(NULL));
	xmin=-dx/2;
	xmax=dx/2;
	ymin=-dy/2;
	ymax=dy/2;
	zmin=-dz/2;
	zmax=dz/2;
	for (i=0;i<N1;i++)
	{
		x[i]=xmin+(xmax-xmin)*((double)rand())/RAND_MAX;
		y[i]=ymin+(ymax-ymin)*((double)rand())/RAND_MAX;
		z[i]=zmin+(zmax-zmin)*((double)rand())/RAND_MAX;
	}
	T=MakeTopology(x, y, z, N1);
	//if (ssdp_error_state)
	for (i=0;i<N2;i++)
	{
		x[i]=xmin+(xmax-xmin)*((double)rand())/RAND_MAX;
		y[i]=ymin+(ymax-ymin)*((double)rand())/RAND_MAX;
		z[i]=SampleTopo(x[i], y[i], &T, NULL);
		if (z[i]<zmin)
			z[i]=zmin;
		if (z[i]>zmax)
			z[i]=zmax;
	}
	free_topo (&T);
	T=MakeTopology(x, y, z, N2);
	free(x);
	free(y);
	free(z);
	return T;
}


void RizeHorizon(horizon *H, double azi1, double azi2, double zen)
{
	int i, j, k;
	if (zen<0)
		return;
	i=(int)round(azi1/H->astep);
	j=(int)round(azi2/H->astep);
	
	if (i<0)
		i+=H->N;
	if (j<0)
		j+=H->N;
		
	i=i%H->N;
	j=j%H->N;
	
	if (i>j)
	{
		k=i;
		i=j;
		j=k;
	}
	if (j-i>=H->N/2)
	{
		for (k=j;k<H->N;k++)
		{
			if (H->zen[k]>zen)
				H->zen[k]=zen;
		}
		for (k=0;k<=i;k++)
		{
			if (H->zen[k]>zen)
				H->zen[k]=zen;
		}
	}
	else
	{	
		for (k=i;k<=j;k++)
		{
			if (H->zen[k]>zen)
				H->zen[k]=zen;
		}
	}
}

// MakeHorizon relies on the triangulation producing right handed triangles
double pcross(double ax, double ay, double bx, double by, double cx, double cy) 
{ 
	return (bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
} 

int EdgeVis(double px, double py, double ax, double ay, double bx, double by)
{
	return (pcross(px,py,ax,ay,bx,by)<0);
}

#define AX x[T.i]
#define AY y[T.i]
#define BX x[T.j]
#define BY y[T.j]
#define CX x[T.k]
#define CY y[T.k]
int TriangleAziRange(triangles T, double *x, double *y, double dz, double xoff, double yoff, double *a1, double *a2)
{
	// assume right handed triangles
	double dx, dy, l, dz2=dz*dz;
	// fprintf(stderr, "Winding %f\n", pcross(AX, AY, BX, BY, CX, CY));
	if (EdgeVis(xoff,yoff,AX,AY,BX,BY)&&(EdgeVis(xoff,yoff,BX,BY,CX,CY))&&(EdgeVis(xoff,yoff,CX,CY,AX,AY))) // right winding with each edge
	{
		(*a1)=-M_PI;
		(*a2)=M_PI;
		return 0; // we are inside the triangle
	}
	if (pcross(xoff,yoff,AX,AY,BX,BY)*(pcross(xoff,yoff,BX,BY,CX,CY))>0) // consitently left or right winding from A via B to C
	{
		/* the cartesian and angular systems used here orientate to how we generally use 
		 * maps and compasses this means that north is up (i.e. the y-axis points north) 
		 * and is at 0°. The x-axis points to the east which is at 90° (i.e. we go clockwise). 
		 * Thus in the following atan2's we swap the usual positions of the x and y arguments.
		 * Note that the definition of angles is also affected by the sunpos.c code.
		 */ 
		// p1 p3
		dx=AX-xoff;
		dy=AY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a1)=ATAN2(l*dx,l*dy); // I want north to be 0° and east
		
		dx=CX-xoff;
		dy=CY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a2)=ATAN2(l*dx,l*dy);
		return 1;
	}
	if (pcross(xoff,yoff,BX,BY,CX,CY)*(pcross(xoff,yoff,CX,CY,AX,AY))>0) // consitently left or right winding from B via C to A
	{
		// p2 p1
		dx=BX-xoff;
		dy=BY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a1)=ATAN2(l*dx,l*dy);
		
		dx=AX-xoff;
		dy=AY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a2)=ATAN2(l*dx,l*dy);
		return 1;
	}
	if (pcross(xoff,yoff,CX,CY,AX,AY)*(pcross(xoff,yoff,AX,AY,BX,BY))>0) // consitently left or right winding from C via A to B
	{
		// p2 p3
		dx=CX-xoff;
		dy=CY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a1)=ATAN2(l*dx,l*dy);
		
		dx=BX-xoff;
		dy=BY-yoff;
		l=sqrt((dx*dx+dz2)/(dx*dx+dy*dy));
		(*a2)=ATAN2(l*dx,l*dy);
		
		return 1;
	}
	// probably the point in question coincides with one of the triangle corners
	return 0;
	// ERRORFLAG TRIANGLEMESS  "Error cannot figure out the azimuth range of a triangle"
	//AddErr(TRIANGLEMESS);
	//(*a1)=-M_PI;
	//(*a2)=M_PI;
	return 0;
}

#undef AX 
#undef AY 
#undef BX 
#undef BY 
#undef CX 
#undef CY 
// take a topology and rize the horizon accordingly
void ComputeHorizon(horizon *H, topology *T, double minzen, double xoff, double yoff, double zoff)
// the minzen parameter can save alot of work
// it specifies the minimum zenith angle to consider a triangle for the horizon
// I would set it to 0.5 times the zenith step in the sky. Especially for large topographies it reduces the amount of work considerably as far away triangles are less likely of consequence
{
	int i;
	double d;
	double a1, a2;	
	double r;
	r=tan(minzen);// compute threshold height over distance ratio
	i=T->Nt-1;
	while ((i>=0)&&(T->T[i].ccz>zoff)) // go backwards through the list till the triangles are lower than the projection point
	{
		d=sqrt((T->T[i].ccx-xoff)*(T->T[i].ccx-xoff)+(T->T[i].ccy-yoff)*(T->T[i].ccy-yoff));
		if ((T->T[i].ccz-zoff)/d>r) // do not compute anything for triangles below the zenith threshold
			if (TriangleAziRange(T->T[i], T->x, T->y, T->T[i].ccz-zoff, xoff, yoff, &a1, &a2)) // compute azimuthal range
				RizeHorizon(H, a1, a2, d/(T->T[i].ccz-zoff)); // for now store the ratio, do atan2's on the horizon array at the end
		i--;
	}	
}


horizon MakeHorizon(sky_grid *sky, topology *T, double xoff, double yoff, double zoff)
{
	horizon H={0,0,NULL};
	H=InitHorizon(sky->Nz);
	if (ssdp_error_state)
		return H;
	ComputeHorizon(&H, T, M_PI/4.0/((double)sky->Nz), xoff, yoff, zoff);
	AtanHorizon(&H);
	return H;
}


int HorizonSobolSet(topogrid *T, int n, double decay)
{
		if ((T->horizon_nsample==n) &&
			(fabs(T->horizon_decay-decay)<1e-5))
				return 0;
		T->horizon_nsample = n;
		T->horizon_decay = decay;

		if (T->horizon_nsample < 0) return -2;
		free(T->horizon_sample);
		T->horizon_sample=calloc((2*T->Nx-1)*(2*T->Ny-1), sizeof(T->horizon_sample));
		if (NULL == T->horizon_sample) goto ecalloc;

		int i, x, y, N=(T->Nx>T->Ny)?T->Nx:T->Ny;
		double p[2];
		struct sobolseq* sq;
		if (NULL==(sq=sobolseq_init(2))) goto esq;

		for (i=0; i < n; ++i) {
				if (sobolseq_gen(sq, p)) goto egen;
				p[0] *= 2*M_PI;
				p[1] = N * pow(p[1], decay);
				x = (int) (p[1]*cos(p[0]));
				y = (int) (p[1]*sin(p[0]));

				if ((abs(x)>=T->Nx)||(abs(y)>=T->Ny))
						continue;

				x = (T->Nx - 1 + x)*(2*T->Ny-1) + T->Ny - 1 + y;
				T->horizon_sample[x] = 1;
		}

		x=0;
		for (i=0; i < 4*T->Nx*T->Ny; ++i)
				if (T->horizon_sample[i]) ++x;

		sobolseq_free(sq);
		return x;
egen:
		T->horizon_nsample = -1;
		sobolseq_free(sq);
esq:
		free(T->horizon_sample);
		T->horizon_sample = NULL;
ecalloc:
		return -1;
}

// take a topology and rize the horizon accordingly
void ComputeGridHorizon(horizon *H, topogrid *T, double minzen, double xoff, double yoff, double zoff)
// the minzen parameter can save alot of work
// it specifies the minimum zenith angle to consider a triangle for the horizon
// I would set it to 0.5 times the zenith step in the sky. Especially for large topographies it reduces the amount of work considerably as far away triangles are less likely of consequence
{
		int i, kl0, k, l, m, n;
		double d, dx, dy, Dx, Dy, a1, a2, r;
		r=tan(minzen);// compute threshold height over distance ratio
		dx=(T->x2-T->x1)/T->Nx;
		dy=(T->y2-T->y1)/T->Ny;
		k=(int)round((xoff-T->x1)/dx);
		l=(int)round((yoff-T->y1)/dy);

		/*
		  if (k<0)
		  {

		  if (k<-1)
		  {
		  // ERRORFLAG LOCNOTINTOPO  "Error: location outside topography"
		  AddErr(LOCNOTINTOPO);
		  return;
		  }
		  k=0;

		  }
		  if (k>=T->Nx)
		  {
		  if (k>T->Nx)
		  {
		  AddErr(LOCNOTINTOPO);
		  return;
		  }
		  k=T->Nx-1;
		  }
		  if (l<0)
		  {
		  if (l<-1)
		  {
		  AddErr(LOCNOTINTOPO);
		  return;
		  }
		  l=0;
		  }
		  if (l>=T->Ny)
		  {
		  if (l>T->Ny)
		  {
		  AddErr(LOCNOTINTOPO);
		  return;
		  }
		  l=T->Ny-1;
		  }*/

		kl0 = (T->Nx-1-k)*(2*T->Ny-1) + T->Ny-1-l;
		i=T->Nx*T->Ny-1;
		while (i>=0) {
				if ((T->horizon_nsample > 0) &&
					(0==T->horizon_sample[kl0+T->horizon_idx[i]])) {
						--i;
						continue;
				}

				// go backwards through the list till the triangles are lower than the projection point
				if (T->z[T->sort[i]]<=zoff) break;

				m=XINDEX(T->sort[i], T->Ny)-k;
				n=YINDEX(T->sort[i], T->Ny)-l;

				if ((abs(m)>=T->Nx)||(abs(n)>=T->Ny))
				{
						--i;
						continue;
				}
				Dx=dx*(double)m;
				Dy=dy*(double)n;
				d=sqrt(Dx*Dx+Dy*Dy);

				if ((T->z[T->sort[i]]-zoff)/d>r) // do not compute anything for triangles below the zenith threshold
						if (Arange(m, n, &a1, &a2, T->A1, T->A2)) // compute azimuthal range
						{
								RizeHorizon(H, (double)a1, (double)a2, d/(T->z[T->sort[i]]-zoff)); // for now store the ratio, do atan2's on the horizon array at the end
						}
				--i;
		}
}
