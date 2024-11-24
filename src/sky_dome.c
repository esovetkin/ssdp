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
#include <hdf5.h>
#include "vector.h"
#include "sky_dome.h"
#include "h5io.h"
#include "error.h"
#include "print.h"


#define UNUSED(x) (void)(x)

// helper routines for indexing hex meshes
int i_sqrt(int x)	// integer sqrt
{
	if(x<=0)
		return 0;
	else
		return (int)sqrt((float)(x));
}
int NNZ(int n)	// number of elements in a mesh given a number of levels (zenith discretizations)
{
	return 3*n*(n+1)+1;
}
int NZN(int n)	// inverse of above, i.e. give it a index and it returns the level (zenith) index
{
	if (n==0)
		return 0;
	return (1+(i_sqrt(12*(n-1)+9)-3)/6);
}

// compute angles for a given index, i, and mesh with Nz zenith disdcretizations
sky_pos GridPos(int Nz, int i)
{
	sky_pos r;
	int Na;
	int nz;
	if (i<=0)
	{
		r.a=0;
		r.z=0;
		return r;
	}
	if (i>NNZ(Nz))
	{
		r.a=0;
		r.z=M_PI/2.0;
		return r;
	}
	nz=NZN(i);
	Na=nz*6;
	r.z=((double)nz)*M_PI/2.0/((double)Nz+0.5);
	//if (nz==Nz)
	//	r.z-=M_PI/8.0/((double)Nz);
	r.a=fmod(2*((double)(i-NNZ(nz-1)))*M_PI/(double)Na, 2*M_PI);
	return r;
}

// compute index for a point p in a mesh with Nz zenith discretizations
// is a bit approximate though, you may end up in the patch next to the one you want.
int GridIndex(int Nz, sky_pos p)
{
	int nz,Na,i;
	double zstep, astep;
	if (p.z<=0)
		return 0;
	if (p.z>M_PI/2)
		return NNZ(Nz)-1;
		
	zstep=M_PI/2.0/((double)Nz+0.5);
	nz=(int)round(p.z/zstep);
	if (nz==0)
		return 0;
	if (nz>Nz)
		nz=Nz;
	Na=nz*6; 
	i=NNZ(nz-1);	
	astep=2*M_PI/((double)Na);
	if (p.a>=0)
		i+=(int)round(p.a/astep);
	else
		i+=(int)round((p.a+2*M_PI)/astep);
	
	if (i>NNZ(Nz)-1)
		i=NNZ(Nz)-1;
	if (i<0)
		i=0;
	return i;
}


// routines to connect patches
static int NextIsoL(int Nz, int index) // same level next
{
	UNUSED(Nz);
	int nz;
	if (index==0)
		return -1;

	nz=NZN(index);
	if (nz==NZN(index+1))
		return index+1;
	return index-nz*6+1;// last of this level, must go to first of this level
}

static int PrevIsoL(int Nz, int index) // same level previous
{
	UNUSED(Nz);
	int nz;
	if (index==0)
		return -1;
		
	nz=NZN(index);
	if (nz==NZN(index-1))// previous at same level
		return index-1;
	return index+nz*6-1;// last of this level, must go to first of this level
}

void NextL(int Nz, int index, int *N) // next level
{
	int nz, i, aoff, n=0;
	nz=NZN(index); // this level
	i=NNZ(nz-1);   // start index of this level
	aoff=index-i;  // offset at this level
	
	if (index==0) // hexpatch 0 only has next level neighbors
	{
		N[n++]=1;
		N[n++]=2;
		N[n++]=3;
		N[n++]=4;
		N[n++]=5;
		N[n++]=6;
		N[n]=-1;
		return;
	}
	if (nz==Nz) // last level has no next
	{
		N[n]=-1;
		return;
	}
	// from the center of the grid there are 3 axis where the hexpatches are perfectly aligned.
	// this means that there are 6 azimuth angles that always repeat for every level
	// this identifies whether we are on one of the 6 major azimuth angles
	if (aoff%nz==0)
	{
		N[n++]=index+(NNZ(nz)-i)+aoff/nz;
		N[n]=N[n-1]+1;
		n++;
		N[n]=N[n-2]-1;
		n++;
		if (NZN(N[n-1])==nz) // this one is two levels down...
			N[n-1]=NNZ(nz+1)-1;
	}
	else
	{
		N[n++]=index+(NNZ(nz)-i)+aoff/nz;
		N[n]=N[n-1]+1;
		n++;
	}
	N[n]=-1;
}

void PrevL(int Nz, int index, int *N) // previous level
{
	UNUSED(Nz);
	int nz, i, aoff, n=0;
	nz=NZN(index); // this level
	i=NNZ(nz-1);   // start index of this level
	aoff=index-i;  // offset at this level

	if (index==0) // hexpatch 0 only has no previous level
	{
		N[n++]=-1;
		return;
	}
	// from the center of the grid there are 3 axis where the hexpatches are perfectly aligned.
	// this means that there are 6 azimuth angles that always repeat for every level
	// this identifies whether we are on one of the 6 major azimuth angles
	if (aoff%nz==0)
	{
		if (nz==1)
			N[n++]=0;
		else
			N[n++]=index-(nz-1)*6-aoff/nz;
		N[n++]=-1;
	}
	else
	{
		N[n++]=index-(nz-1)*6-aoff/nz-1;
		N[n]=N[n-1]+1;
		n++;
		if (NZN(N[n-1])==nz)
			N[n-1]=NNZ(nz-1);
	}
	N[n]=-1;
}


// for debugging meshing some routines to export connections
void Neighbors(int Nz, int index, int *N) // collect all neighbors
{
	int *n;
	n=N;
	(*n)=PrevIsoL(Nz, index);
	if(*n>=0)
		n++;
	
	(*n)=NextIsoL(Nz, index);
	if (*n>=0)
		n++;	
		
	PrevL(Nz, index, n);
	while (*n>=0)
		n++;	
		
	NextL(Nz, index, n);
	while (*n>=0)
		n++;	
}

// used for degugging, see if patches know their neighbors
void PlotConn(int Nz, int N[], int index)
{
	sky_pos p0,p;
	vec v0, v;
	int i;
	p0=GridPos(Nz, index);
	v0=unit(p0);
	i=0;
	while (N[i]>=0)
	{
		p=GridPos(Nz, N[i]);
		v=diff(unit(p), v0);
		// gnuplot compatible vectors (splot '...' u 1:2:3:4:5:6 w vectors)
		printf("%e %e %e %e %e %e\n", v0.x, v0.y, v0.z, v.x, v.y, v.z);
		i++;
	}
}	
void Connectivity(int Nz)
{
	int i;
	int N[7];
	for (i=0;i<NNZ(Nz);i++)
	{
		Neighbors(Nz, i, N);
		PlotConn(Nz,N,i);
	}
}


// findpatch routine
int FindPatch(sky_grid *sky, sky_pos p)
{
	int r;
	int i,j;
	int dmin, d, rmin;
	r=GridIndex(sky->Nz, p);
	rmin=r;
	// result of GridIndex(...) may be a bit approximate
	// look through the neigbors for the best patch
	dmin=(p.a-sky->P[r].p.a)*(p.a-sky->P[r].p.a)+(p.z-sky->P[r].p.z)*(p.z-sky->P[r].p.z);
	i=0;
	while ((j=sky->P[r].NL[i])>=0)
	{
		d=(p.a-sky->P[j].p.a)*(p.a-sky->P[j].p.a)+(p.z-sky->P[j].p.z)*(p.z-sky->P[j].p.z);
		if (d<dmin)
		{
			dmin=d;
			rmin=j;
		}
		i++;
	}
	i=0;
	while ((j=sky->P[r].PL[i])>=0)
	{
		d=(p.a-sky->P[j].p.a)*(p.a-sky->P[j].p.a)+(p.z-sky->P[j].p.z)*(p.z-sky->P[j].p.z);
		if (d<dmin)
		{
			dmin=d;
			rmin=j;
		}
		i++;
	}
	if ((j=sky->P[r].PI)>=0)
	{
		d=(p.a-sky->P[j].p.a)*(p.a-sky->P[j].p.a)+(p.z-sky->P[j].p.z)*(p.z-sky->P[j].p.z);
		if (d<dmin)
		{
			dmin=d;
			rmin=j;
		}
	}
	if ((j=sky->P[r].NI)>=0)
	{
		d=(p.a-sky->P[j].p.a)*(p.a-sky->P[j].p.a)+(p.z-sky->P[j].p.z)*(p.z-sky->P[j].p.z);
		if (d<dmin)
		{
			dmin=d;
			rmin=j;
		}
	}
	return rmin;	
}


double SolidAngle(int N, int Nz, int i)
{
	int nz;
	double z1,z2, dz;
	/* correction factor for the solid angle of a hexpatch
	 * The zenith step is:
	 *  
	 * dθ=π/(2Nz)
	 * 
	 * The solid angle of the complete band of patches at one particular 
	 * zenith angle θ equals:
	 * 
	 * Ω=2π(cos(θ-dθ/2)-cos(θ+dθ/2))
	 * 
	 * This solid angle is the solid angle of 6 nz hexpatches, i.e.
	 * the hexpatches have a solid angle of 
	 * 
	 * Ω=π(cos(θ-dθ/2)-cos(θ+dθ/2))/(3 nz),
	 * 
	 * where
	 * 
	 * θ=nz π/(2Nz)
	 * 
	 * For the patch at nz=0 we get
	 * 
	 * Ω=2π(1-cos(dθ/2))
	 * 
	 * In this routine we return the normalized solid angle w.r.t. 
	 * the average solid angle 2π/N with N the total number of patches
	 * 
	 */ 
	// obviously not efficient to recompute thinbs every time
	// however, we do not need to initialize skies so often
	if (i>N)
		return 0;
	nz=NZN(i);
	dz=M_PI/2.0/((double)Nz+0.5);
	if (nz==0)
		return ((double)N)*(1-cos(dz/2));
	/*if (nz==Nz)
	{
		z1=M_PI/2-dz/2;
		return ((double)N)*cos(z1)/6.0/((double)Nz);
	}*/
	
	z1=((double)nz-0.5)*dz;
	z2=z1+dz;
	return ((double)N)*(cos(z1)-cos(z2))/6.0/((double)nz);
}

// as the number of elements rises quadratically we better set a
// practical limit on the Nz values
// if you ever change the 1 Mio below you should consider modifying 
// ZENEPS in project.c.
#define MAXNZ NZN(1000000)
sky_grid InitSky(int Nz)
{
	sky_grid sky;
	sky_pos sun={0,0};// default sun, straight above
	int i;
	if (Nz>MAXNZ)
	{
		Print(WARNING,"Warning: number of zenith discretizations %d too large, using %d instead\n", Nz, MAXNZ);
		Nz=MAXNZ;
	}
	sky.N=NNZ(Nz);	
	sky.Nz=Nz;
	sky.sp=sun; // The default sun: is strainght above
	sky.sI=0;	// and pitch black
	sky.suni=0;
	// ERRORFLAG MALLOCFAILSKYDOME  "Error memory allocation failed in creating a sky-dome"
	if ((sky.P=malloc((sky.N+1)*sizeof(hexpatch)))==NULL)
	{
		AddErr(MALLOCFAILSKYDOME);
		sky.Nz=0;
		sky.N=0;
		return sky;
	}
	sky.h5b_P.len = sky.N+1;
	sky.h5b_P.p = sky.P;
	if ((sky.cosz=malloc((sky.N+1)*sizeof(double)))==NULL)
	{
		free(sky.P);
		sky.P=NULL;
		AddErr(MALLOCFAILSKYDOME);
		sky.Nz=0;
		sky.N=0;
		return sky;
	}
	sky.h5b_cosz.len = sky.N+1;
	sky.h5b_cosz.p = sky.cosz;
	if ((sky.sa=malloc((sky.N+1)*sizeof(double)))==NULL)
	{
		free(sky.P);
		sky.P=NULL;
		AddErr(MALLOCFAILSKYDOME);
		sky.Nz=0;
		sky.N=0;
		return sky;
	}
	sky.h5b_sa.len = sky.N+1;
	sky.h5b_sa.p = sky.sa;
	sky.icosz=0;
	
	for (i=0;i<sky.N;i++)
	{
		sky.P[i].I=0;
		sky.P[i].p=GridPos(Nz, i);
		NextL(Nz, i, sky.P[i].NL);
		PrevL(Nz, i, sky.P[i].PL);
		sky.P[i].NI=NextIsoL(Nz, i);
		sky.P[i].PI=PrevIsoL(Nz, i);
		sky.cosz[i]=cos(sky.P[i].p.z);
		sky.sa[i]=SolidAngle(sky.N, Nz, i);
		sky.icosz+=sky.cosz[i]*sky.sa[i]; // normalization constant
	}
	return sky;
}

void free_sky_grid(sky_grid *sky)
{
	if (sky->P)
		free(sky->P);
	sky->P=NULL;
	sky->h5b_P.len=0;
	sky->h5b_P.p=NULL;
	if (sky->cosz)
		free(sky->cosz);
	sky->cosz=NULL;	
	sky->h5b_cosz.len=0;
	sky->h5b_cosz.p=NULL;
	if (sky->sa)
		free(sky->sa);
	sky->sa=NULL;	
	sky->h5b_sa.len=0;
	sky->h5b_sa.p=NULL;
	sky->Nz=0;
	sky->N=0;
}


void ZenithRange(const sky_grid *sky, double minz, double maxz,
				 int* imin, int *imax)
{
		// here it is guaranteed that all patches
		//
		// for i >= imax => patch[i].z > maxz
		// for i <= imin => patch[i].z < minz

		minz *= 2*((double)sky->Nz + 0.5) / M_PI;
		maxz *= 2*((double)sky->Nz + 0.5) / M_PI;

		*imin = NNZ((int)floor(minz)) - 1;
		*imax = NNZ((int)ceil(maxz) - 1);

		if (*imin < 0) *imin = -1;
		if (*imax < 0) *imax = 0;
		if (*imax >= sky->N) *imax = sky->N;
		if (*imin >= sky->N) *imin = sky->N;
}


hid_t h5t_hexpatch()
{
		hid_t t_nl, t_pl, t_res, t_sky_pos;

		t_nl=H5Tarray_create(H5T_NATIVE_INT,1,(hsize_t[]){7});
		if (H5I_INVALID_HID==t_nl) goto et_nl;
		t_pl = H5Tarray_create(H5T_NATIVE_INT,1,(hsize_t[]){3});
		if (H5I_INVALID_HID==t_pl) goto et_pl;
		t_sky_pos = h5t_sky_pos();
		if (H5I_INVALID_HID==t_sky_pos) goto et_sky_pos;

		t_res = H5Tcreate(H5T_COMPOUND, sizeof(hexpatch));
		if (H5I_INVALID_HID==t_res) goto et_res;
		if (0>H5Tinsert(t_res, "I", HOFFSET(hexpatch, I), H5T_NATIVE_DOUBLE)) goto einsert;
		if (0>H5Tinsert(t_res, "p", HOFFSET(hexpatch, p), t_sky_pos)) goto einsert;
		if (0>H5Tinsert(t_res, "NL", HOFFSET(hexpatch, NL), t_nl)) goto einsert;
		if (0>H5Tinsert(t_res, "PL", HOFFSET(hexpatch, PL), t_pl)) goto einsert;
		if (0>H5Tinsert(t_res, "NI", HOFFSET(hexpatch, NI), H5T_NATIVE_INT)) goto einsert;
		if (0>H5Tinsert(t_res, "PI", HOFFSET(hexpatch, PI), H5T_NATIVE_INT)) goto einsert;

		H5Tclose(t_sky_pos);
		H5Tclose(t_pl);
		H5Tclose(t_nl);
		return t_res;
einsert:
		H5Tclose(t_res);
et_res:
		H5Tclose(t_sky_pos);
et_sky_pos:
		H5Tclose(t_pl);
et_pl:
		H5Tclose(t_nl);
et_nl:
		return H5I_INVALID_HID;
}


hid_t h5t_sky_grid()
{
		hid_t t_vld, t_res, t_hexpatch, t_sky_pos, t_vl_hexpatch;
		t_vld = H5Tvlen_create(H5T_NATIVE_DOUBLE);
		if (H5I_INVALID_HID==t_vld) goto et_vld;
		t_hexpatch=h5t_hexpatch();
		if (H5I_INVALID_HID==t_hexpatch) goto et_hexpatch;
		t_vl_hexpatch=H5Tvlen_create(t_hexpatch);
		if (H5I_INVALID_HID==t_vl_hexpatch) goto et_vl_hexpatch;
		t_sky_pos = h5t_sky_pos();
		if (H5I_INVALID_HID==t_sky_pos) goto et_sky_pos;

		t_res = H5Tcreate(H5T_COMPOUND, sizeof(sky_grid));
		if (H5I_INVALID_HID==t_res) goto et_res;
		if (0>H5Tinsert(t_res, "P", HOFFSET(sky_grid, h5b_P), t_vl_hexpatch)) goto einsert;
		if (0>H5Tinsert(t_res, "sp", HOFFSET(sky_grid, sp), t_sky_pos)) goto einsert;
		if (0>H5Tinsert(t_res, "sI", HOFFSET(sky_grid, sI), H5T_NATIVE_DOUBLE)) goto einsert;
		if (0>H5Tinsert(t_res, "suni", HOFFSET(sky_grid, suni), H5T_NATIVE_INT)) goto einsert;
		if (0>H5Tinsert(t_res, "cosz", HOFFSET(sky_grid, h5b_cosz), t_vld)) goto einsert;
		if (0>H5Tinsert(t_res, "sa", HOFFSET(sky_grid, h5b_sa), t_vld)) goto einsert;
		if (0>H5Tinsert(t_res, "icosz", HOFFSET(sky_grid, icosz), H5T_NATIVE_DOUBLE)) goto einsert;
		if (0>H5Tinsert(t_res, "N", HOFFSET(sky_grid, N), H5T_NATIVE_INT)) goto einsert;
		if (0>H5Tinsert(t_res, "Nz", HOFFSET(sky_grid, Nz), H5T_NATIVE_INT)) goto einsert;

		H5Tclose(t_sky_pos);
		H5Tclose(t_vl_hexpatch);
		H5Tclose(t_hexpatch);
		H5Tclose(t_vld);
		return t_res;
einsert:
		H5Tclose(t_res);
et_res:
		H5Tclose(t_sky_pos);
et_sky_pos:
		H5Tclose(t_vl_hexpatch);
et_vl_hexpatch:
		H5Tclose(t_hexpatch);
et_hexpatch:
		H5Tclose(t_vld);
et_vld:
		return H5I_INVALID_HID;
}


int h5write_sky_grid(sky_grid* data, const char* ofn, const char *dataset)
{
		hid_t dsp, file, dst, t_sky_grid;

		if (H5I_INVALID_HID == (file=h5io_fopen(ofn, 0))) goto efile;
		if (0 != h5_datasetisin(file, dataset)) goto edataset;

		dsp = H5Screate_simple(1, (hsize_t[]){1}, NULL);
		if (H5I_INVALID_HID == dsp) goto edsp;
		t_sky_grid = h5t_sky_grid();
		if (H5I_INVALID_HID == t_sky_grid) goto et_sky_grid;

		dst = H5Dcreate(file, dataset, t_sky_grid, dsp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (H5I_INVALID_HID == dst) goto edst;
		H5Oset_comment(dst, "SSDP data: write_sky");

		if (0>H5Dwrite(dst, t_sky_grid, H5S_ALL, H5S_ALL,	H5P_DEFAULT, data)) goto ewrite;

		H5Dclose(dst);
		H5Tclose(t_sky_grid);
		H5Sclose(dsp);
		H5Fclose(file);
		return 0;
ewrite:
		H5Dclose(dst);
edst:
		H5Tclose(t_sky_grid);
et_sky_grid:
		H5Sclose(dsp);
edsp:
edataset:
		H5Fclose(file);
efile:
		return -1;
}


int h5read_sky_grid(sky_grid* data, const char* ifn, const char *dataset)
{
		hid_t file, dst, t_sky_grid, t_ftype;

		file=H5Fopen(ifn, H5F_ACC_RDONLY, H5P_DEFAULT);
		if (H5I_INVALID_HID == file) goto efile;
		if (1 > h5_datasetisin(file, dataset)) goto edataset;

		dst = H5Dopen(file, dataset, H5P_DEFAULT);
		if (H5I_INVALID_HID == dst) goto edst;

		if (H5I_INVALID_HID==(t_ftype=H5Dget_type(dst))) goto et_ftype;
		if (H5I_INVALID_HID==(t_sky_grid=h5t_sky_grid())) goto et_sky_grid;

		htri_t status = H5Tequal(t_sky_grid, t_ftype);
		if (status < 0) goto efaileq;
		if (0 == status) goto enonmatch;

		// data must be freed before
		// NO! free_sky_grid(data);
		if (0>H5Dread(dst, t_sky_grid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) goto eread;

		data->P = data->h5b_P.p;
		data->cosz = data->h5b_cosz.p;
		data->sa = data->h5b_sa.p;

		H5Tclose(t_sky_grid);
		H5Tclose(t_ftype);
		H5Dclose(dst);
		H5Fclose(file);
		return 0;
eread:
enonmatch:
efaileq:
		H5Tclose(t_sky_grid);
et_sky_grid:
		H5Tclose(t_ftype);
et_ftype:
		H5Fclose(dst);
edst:
edataset:
		H5Fclose(file);
efile:
		return -1;
}


#ifdef RUNTEST

#include <signal.h>
#include <assert.h>

/** comment #include <assert.h> and use this one for gdb

void assert(int v)
{
		if (v)
				return;

		raise(SIGTRAP);
}
*/

static void test_range(sky_grid* sky, double min, double max, int imin, int imax)
{
		int i;
		assert(imin < imax);
		assert(imin >= -1);
		assert(imax <= sky->N);

		for (i=0; i <= imin; ++i)
				assert((sky->P[i].p.z < min));
		assert(sky->P[i].p.z >= min);

		for (i=imax; i < sky->N; ++i)
				assert((sky->P[i].p.z > max));
		assert(sky->P[imax-1].p.z <= max);
}


static void test_zenithrange(int Nz, int nzens, int nx)
{
		double dx = 1.0 / ((double) nzens + 1);
		int i, imin, imax;
		sky_grid sky = InitSky(Nz);

		for (i=1; i < nzens; ++i) {
				double x = ((double)i) * dx;
				ZenithRange(&sky, x, x + nx*dx, &imin, &imax);
				test_range(&sky, x, x+nx*dx, imin, imax);
		}

		free_sky_grid(&sky);
}


int main(void)
{
		signal(SIGTRAP, SIG_IGN);

		printf("testing sky_dome ...\n");
		test_zenithrange(21, 10, 0);
		test_zenithrange(21, 10, 1);
		test_zenithrange(21, 10, 2);
		test_zenithrange(21, 100, 2);
		test_zenithrange(21, 100, 20);
		test_zenithrange(21, 10, 10);

		printf("PASSED\n");
}

#endif
