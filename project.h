#ifndef _RROJECT_H
#define _RROJECT_H
double DiffusePlaneOfArray(sky_grid sky, double tilt, double a, int mask);
double PlaneOfArray(sky_grid sky, double tilt, double a, int mask);
double DiffusteHorizontal(sky_grid sky, int mask);
double GlobalHorizontal(sky_grid sky, int mask);
double POA_Albedo(sky_grid sky, double albedo, double tilt, double a, int mask);
#endif
