#include <math.h>
#include <float.h>
// fast atan2 approx
// subatan: approximate atan(x) on the interval x on [0-1]
float subatan(float a)
{
	// approximation from some forum
	float s;
	s=a*a;
	return ((-0.0464964749 * s + 0.15931422) * s - 0.327622764) * s * a + a;
}
float subatan__(float a)
{
	// home brew second order approximation
	return (-2.7874098272920872e-01*a+1.0641391461266569e+00)*a;
}
float subatan_(float a)
{
	// first order
	return 7.85398163397448309628e-01*a;
}

double fatan2(double y, double x)
{
	float fx, fy, a, r;
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
	r=subatan(a);
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
	float a, r;
	if (x>FLT_MAX)
		return M_PI/2;
	if (x<-FLT_MAX)
		return -M_PI/2;
		
	a=(float)fabs(x);
	r=subatan(a);
	if (!isfinite(r))
	{
		if (x>0)
			return M_PI/2;
		else
			return -M_PI/2;
	}	
	if (x<0)
		r=-r;
	return (double)r;
}
