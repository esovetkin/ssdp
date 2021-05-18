#include <math.h>
// fast atan2 approx
double fatan2(double y, double x)
{
	float fx, fy, a, s, r;
	if ((y==0)&&(x>=0))
		return 0;
	if ((y==0)&&(x<0))
		return M_PI;
	if ((x==0)&&(y>=0))
		return M_PI/2;
	if ((x==0)&&(y<0))
		return 3*M_PI/2;		
	fx=fabs((float)x);
	fy=fabs((float)y);
	if (fx<fy)
		a=fx/fy;
	else
		a=fy/fx;
	s=a*a;
	r=((-0.0464964749 * s + 0.15931422) * s - 0.327622764) * s * a + a;
	if (fy>fx)
		r=M_PI/2-r;
	if (x<0)
		r=M_PI-r;
	if (y<0)
		r=-r;
	return (double)r;
}
double fatan(double x)
{
	float a, s, r;
	a=(float)fabs(x);
	s=a*a;
	r=((-0.0464964749 * s + 0.15931422) * s - 0.327622764) * s * a + a;
	if (x<0)
		r=-r;
	return (double)r;
}
