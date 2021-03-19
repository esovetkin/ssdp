/* Sloppy Sky Dome Projector
 * 
 */
#ifndef _SKY_DOME_H
#define _SKY_DOME_H

typedef struct hexpatch {
	double I;
	sky_pos p;
	char mask; // mark elements behind the horizon
	int NL[7]; // next level neighbors
	int PL[3]; // previous level neighbors
	int NI;    // iso level next
	int PI;    // iso level previous
} hexpatch;

typedef struct sky_grid {
	hexpatch *P;
	// get sun out of the hexpatch dome so we can separate contributions
	// if we want
	sky_pos sp;	// solar position
	double sI;	// solar intensity
	char smask;
	int N;		
	int Nz;
} sky_grid;

void Connectivity(int Nz);
sky_grid InitSky(int Nz);
void free_sky_grid(sky_grid *sky);
int FindPatch(sky_grid sky, sky_pos p);
#define SUNSR 6.807e-5
#endif
