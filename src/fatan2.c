#include <math.h>
#include <float.h>
// fast atan2 approx
/* my approximation approach
 * We only need atan on an interval from 0-1 (the subatan function)
 * I use polynomials of the form a1*x^n + a2*x^n-1 + ... a1 x + a0 
 * I fix a0=0 and a1=1 as that is probably not bad anyway
 * Then I fix the number of terms and maximum order and simply try all 
 * combinations and select the best precision
*/
float subatan_ref(float a)
{
	float s;
	// my reference: approx with 3 terms, max order 6 Emax: 0.0001401683 
	s=a*a;	
	return ((-0.095375152306*a+0.212415310139)*s -0.331782140628)*s*a+a;
}

float subatan_4_9(float a)
{
	float s2, s3;
	// 4 terms, max order 9 Emax: 0.0000035196 (Erel: 0.025110 trel: 1.17)
	s2=a*a;
	s3=s2*a;
	return (((0.006304618756*s3 -0.120631272109)*a + 0.234063762657)*s2 -0.334335841774)*s3+a;
}

float subatan_5_9(float a)
{
	float s;
	// 5 terms, max order 9 Emax: 0.0000003907 (Erel: 0.002787 trel: 1.28)
	s=a*a;
	return ((((-0.033340276431*a+0.140924407553)*a -0.191655659596)*s + 0.202922388832)*s -0.333452613666)*s*a+a;
}

#define subatan subatan_5_9
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
	float fx, a, r;
	if (x>FLT_MAX)
		return M_PI/2;
	if (x<-FLT_MAX)
		return -M_PI/2;		
		
	fx=(float)fabs(x);
	if (fx>1.0)
		a=1.0/x;
	else
		a=fx;
	r=subatan(a);
	if (fx>1.0)
		r=M_PI/2-r;
	if (x<0)
		r=M_PI-r;
	return (double)r;
}
