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
// 2 terms order 4
static inline double subatan_2_4(double a) // 3 digits
{
	return (0.152578515549353*a -0.366840729424806)*a*a*a+a;
}

static inline double subatan_5_9(double a) // 6 digits
{
	return ((((-0.033340276431366*a+0.140924407553154)*a -0.191655659596189)*a*a + 0.202922388831609)*a*a -0.333452613665954)*a*a*a+a;
}

static inline double subatan_9_10(double a) // 8 digits
{
	return ((((((((0.021320677603614*a -0.130029802390657)*a + 0.321042974564416)*a -0.368370428974842)*a + 0.096998715843775)*a + 0.173898218191358)*a + 0.004227037834855)*a -0.333702137561143)*a + 0.000012909580497)*a*a+a;
}

#define subatan subatan_5_9
double fatan2(double y, double x)
{
	double fx, fy, a, r;
	if ((y==0)&&(x>=0))
		return 0;
	if ((y==0)&&(x<0))
		return M_PI;
	if ((x==0)&&(y>=0))
		return M_PI/2;
	if ((x==0)&&(y<0))
		return 3*M_PI/2;		
	fx=fabs(x);
	fy=fabs(y);
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
	return r;
}
double fatan(double x)
{
	double fx, a, r;
	if (x>FLT_MAX)
		return M_PI/2;
	if (x<-FLT_MAX)
		return -M_PI/2;		
		
	fx=fabs(x);
	if (fx>1.0)
		a=1.0/x;
	else
		a=fx;
	r=subatan(a);
	if (fx>1.0)
		r=M_PI/2-r;
	if (x<0)
		r=M_PI-r;
	return r;
}
