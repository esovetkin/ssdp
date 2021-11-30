#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "error.h"
#include "vector.h"

/* implementation of the PSA algorithm as described in
 * Blanco-Muriel, Manuel, et al. "Computing the solar vector." Solar energy 70.5 (2001): 431-441.
 * Another implementation is at: http://www.psa.es/sdg/sunpos.htm
 * I reimplemented as I am not sure about the license conditions for that version
 */
// Julian day noon 1 January 2000 Universal Time
#define JD0 2451545.0
#define EARTH_R 6371.01 //km
#define ASTR_UNIT 149597890 //km

typedef struct JulianDate {
	double JD, dJD, hour;
} JulianDate;

// compute the julain date from epoch time
JulianDate MakeJulianDate(time_t t)
{
	struct tm *ut;
	long int jumon;
	long int juday;
	JulianDate JD;
	
	ut=gmtime(&t);	
	// ERRORFLAG GMTIMENULL  "Error: gmtime returned NULL"
	if (!ut)
	{
		AddErr(GMTIMENULL);
		JD.JD=JD0;
		JD.dJD=0.0;
		JD.hour=0;
		return JD;
	}
	// Calculate time of the day in UT decimal hours
	JD.hour = (double)ut->tm_hour + ((double)ut->tm_min+((double)ut->tm_sec)/60.0)/60.0;
	// Calculate current Julian Day
	// note tm_mon is months since january, i.e. it goes from 0 to 11.
	jumon =(ut->tm_mon-13)/12;
	// note tm_year is years since 1900
	juday=(1461*(ut->tm_year+1900+4800+jumon))/4+(367*(ut->tm_mon-1-12*jumon))/12-(3*((ut->tm_year+1900+4900+jumon)/100))/4+ut->tm_mday-32075;
	JD.JD=((double)(juday))-0.5+JD.hour/24.0;
	// Day since JD0 (1 January 2000)
	JD.dJD=JD.JD-JD0;
	return JD;
}

double Omega(JulianDate J)
{
	return 2.1429-0.0010394594*J.dJD;
}

double MeanLong(JulianDate J)
{
	return 4.8950630+0.017202791698*J.dJD;
}

double MeanAnomaly(JulianDate J)
{
	return 6.2400600+0.0172019699*J.dJD;
}

double EclipticLong(JulianDate J)
{
	double ma=MeanAnomaly(J);
	return MeanLong(J)+0.03341607*sin(ma)+0.00034894*sin(2*ma)-0.0001134-0.0000203*sin(Omega(J));
}

double EclipticOblique(JulianDate J)
{
	return 0.4090928-6.2140e-9*J.dJD+0.0000396*cos(Omega(J));
}

double RightAscension(JulianDate J)
{
	double el=EclipticLong(J);
	double ra;
	ra=atan2(cos(EclipticOblique(J))*sin(el),cos(el));
	if (ra<0)
		ra+=2*M_PI;
	return ra;
}

double Declination(JulianDate J)
{
	return asin(sin(EclipticOblique(J))*sin(EclipticLong(J)));
}

double lmst(JulianDate J, double longitude)
{ 
	//Greenwich Mean Sidereal Time
	double gmst=6.6974243242+0.0657098283*J.dJD+J.hour;
	//Local Mean Sidereal Time
	return deg2rad(gmst*15)+longitude;
}

double hourangle(JulianDate J, double longitude)
{ 
	return lmst(J, longitude)-RightAscension(J);
}

sky_pos sunpos(time_t t, double lat, double lon)
{
	sky_pos s;
	JulianDate J;
	double ha, de;
	J=MakeJulianDate(t);
	ha=hourangle(J,lon);
	de=Declination(J);
	s.z=acos(cos(lat)*cos(ha)*cos(de)+sin(de)*sin(lat));
	s.a=atan2(-sin(ha),tan(de)*cos(lat)-sin(lat)*cos(ha));
	s.a=fmod(s.a,2*M_PI);
	s.z+=sin(s.z)*EARTH_R/ASTR_UNIT; // add parallax
	return s;
}



// bisection within a day to second accuracy is guaranteed to finish 
// within 17 iterations
// computes the nearest solar noon
time_t solar_noon(time_t t, double lat, double lon, sky_pos *P)
{
	sky_pos pmin, pmax, po, p;
	time_t tmin, tmax, to;
	to=t;
	tmin=to-60*60*12;
	tmax=to+60*60*12;
	pmin=sunpos(tmin, lat, lon);
	pmax=sunpos(tmax, lat, lon);
	po=sunpos(to, lat, lon);
	if ((pmin.z<po.z)&&(pmin.z<pmax.z))
	{
		tmax=to;
		pmax=po;
		to=tmin;
		po=pmin;
		tmin-=60*60*12;
		pmin=sunpos(tmin, lat, lon);
	}
	else if ((pmax.z<po.z)&&(pmax.z<pmin.z))
	{
		tmin=to;
		pmin=po;
		to=tmax;
		po=pmax;
		tmax+=60*60*12;
		pmax=sunpos(tmax, lat, lon);
	}
	while(tmax-tmin>1)
	{
		t=(tmin+to)/2;
		p=sunpos(t, lat, lon);
		if (p.z<po.z)
		{
			pmax=po;
			tmax=to;
			po=p;
			to=t;
		}
		else
		{
			pmin=p;
			tmin=t;
		}
		t=(tmax+to)/2;
		p=sunpos(t, lat, lon);
		if (p.z<po.z)
		{
			pmin=po;
			tmin=to;
			po=p;
			to=t;
		}
		else
		{
			pmax=p;
			tmax=t;
		}
	}
	if (P)
		(*P)=po;
	return to;
}

// computes the next sun set takes the solar noon time
time_t sunset(time_t t, double lat, double lon, sky_pos *P)
{
	sky_pos  p;
	time_t tmin, tmax;
	tmin=t;
	tmax=tmin+60*60*12;
	
	while(tmax-tmin>1)
	{
		t=(tmin+tmax)/2;
		p=sunpos(t, lat, lon);
		if (fabs(p.z)<M_PI/2)
			tmin=t;
		else
			tmax=t;
	}
	if (P)
		(*P)=p;
	return t;
}

// computes the previous sun rise, takes solar noon time
time_t sunrise(time_t tnoon, double lat, double lon, sky_pos *P)
{
	sky_pos p;
	time_t tmin, tmax, t;
	tmax=tnoon;
	tmin=tmax-60*60*12;
	
	t=(tmin+tmax)/2;
	while(tmax-tmin>1)
	{
		p=sunpos(t, lat, lon);
		if (fabs(p.z)<M_PI/2)
			tmax=t;
		else
			tmin=t;
		t=(tmin+tmax)/2;
	}
	if (P)
		(*P)=p;
	return t;
}
// computes nearest solar noon, sunrise, and sunset
void SolarTimes(time_t t, double lat, double lon, time_t * trise, time_t *tnoon, time_t *tset, sky_pos *prise, sky_pos *pnoon, sky_pos *pset)
{
	time_t tr, tn, ts;
	tn=solar_noon(t, lat, lon, pnoon);
	tr=sunrise(tn, lat, lon, prise);
	ts=sunset(tn, lat, lon, pset);
	if (tnoon)
		(*tnoon)=tn;
	if (trise)
		(*trise)=tr;
	if (tset)
		(*tset)=ts;
}
