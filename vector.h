#ifndef _VECTOR_H
#define _VECTOR_H
// basic 3D vector calculation
typedef struct sky_pos {
	double z,a;
} sky_pos;

typedef struct vec {
	double x,y,z;
} vec;
vec cross(vec a, vec b);
vec sum(vec a, vec b);
vec scalevec(vec a, double b);
double dot(vec a, vec b);
double norm(vec v);
vec unit(sky_pos a);
sky_pos vecdir(vec a);
double amean(double a1, double a2);
double adiff(double a1, double a2);
sky_pos rrf(sky_pos p, sky_pos axis, double beta);
double rad2degr(double rad);
double degr2rad(double degr);
#endif
