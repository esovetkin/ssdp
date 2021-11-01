#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* simple model for angle of incidence effects on PV modules */

double snellius(double theta_i, double n1, double n2)
{
	return asin(n1*sin(theta_i)/n2);
}

double fresnell(double n1,double n2,double theta_i)
{
	double theta_t, Rs, Rp;
	theta_t=snellius(theta_i, n1, n2);
	Rs=fabs((n1*cos(theta_i)-n2*cos(theta_t))/(n1*cos(theta_i)+n2*cos(theta_t)));
	Rs=Rs*Rs;
	Rp=fabs((n1*cos(theta_t)-n2*cos(theta_i))/(n1*cos(theta_t)+n2*cos(theta_i)));
	Rp=Rp*Rp;
	return (Rs+Rp)/2;
}
// compute transmission 
double Transmission(double n1,double n2,double theta_i)
{
	return 1-fresnell(n1,n2,theta_i);
}

// compute transmission with an antireflection coating
double Transmission_ar(double n1,double n2,double n3, double theta_i)
{
	double I0;
	double r1, r, theta_2;
	I0=1-fresnell(n1,n2,theta_i); // intensity in AR layer after first principal transmission n1->n2
	theta_2=snellius(theta_i, n1, n2); // angle with in AR 	
	// reflection at n2/n3 interface
	r1=fresnell(n2,n3,theta_2);
	r=r1*fresnell(n2,n1,theta_2);/* fraction intensity remaining after one back-forth bounce */
	return I0*(1-r1)/(1-r); // total transmission
}
