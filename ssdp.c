/*
    Simple Sky-Dome Projector
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
/*****************************************************************              
 *  INSTITUT FUER ENERGIE- UND KLIMAFORSCHUNG                    *              
 *  IEK-5 PHOTOVOLTAIK                                           *              
 *                                                               *              
 *        ########                _   _                          *              
 *     ##########                |_| |_|                         *              
 *    ##########     ##         _ _   _ _     ___ ____ _   _     *              
 *   ##########     ####       | | | | | |   |_ _/ ___| | | |    *              
 *   #########     #####    _  | | | | | |    | | |   | |_| |    *              
 *   #    ###     ######   | |_| | |_| | |___ | | |___|  _  |    *              
 *    #          ######     \___/ \___/|_____|___\____|_| |_|    *              
 *     ##      #######      F o r s c h u n g s z e n t r u m    *              
 *       ##########                                              *              
 *                                                               *              
 *   http://www.fz-juelich.de/iek/iek-5/DE/Home/home_node.html   *              
 *****************************************************************
 *                                                               *
 *    Dr. Bart E. Pieters 2021                                   *
 *                                                               *             
 *****************************************************************/                                                                             

/*****************************************************************           
 * FUNCTION:                                                     *                
 * main program						                             *
 *                                                               *            
 *****************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "libssdp.h"
#include "io.h"
#include "util.h"

#define rad2degr(rad) ((rad)*180/M_PI)
#define degr2rad(degr) ((degr)*M_PI/180)
#define MAXSTRLEN 1028
topology LoadTopo(char *fn)
{
	char c, *line;
	int k, Na=50;
	FILE *f;
	int ddef=0;
	if ((f=fopen(fn,"rb"))==NULL)
		Fatal("Cannot open %s for reading\n", fn);
	topology T;	
	line=malloc(MAXSTRLEN*sizeof(char));
    fgets(line, MAXSTRLEN-1, f);
	T.N=0;
	T.x=malloc(Na*sizeof(double));
	T.y=malloc(Na*sizeof(double));
	T.z=malloc(Na*sizeof(double));
	while(feof(f)==0)
	{
		
    	k=sscanf(line, " %c", &c);
		if((k==1)&&(c!='#'))
		{
			k=sscanf(line, " %le %le %le", T.x+T.N, T.y+T.N, T.z+T.N);
			if(k==3)
			{
				T.N++;
				if (Na-1==T.N)
				{
					Na+=50;
					T.x=realloc(T.x, Na*sizeof(double));
					T.y=realloc(T.y, Na*sizeof(double));
					T.z=realloc(T.z, Na*sizeof(double));	
				}
			}
		}
		else if (T.N==0)
		{
			k=sscanf(line, "# d=%le", &T.d);
			if (k==1)
				ddef=1;
				
		}
    	fgets(line, MAXSTRLEN-1, f);
	}
	free(line);
	fclose(f);
	if (!ddef)
	{
		Warning("Warning: no point diameter d specified, setting to 1\n");
		T.d=1;
	}
	return T;
}
// tic-toc timer
clock_t tic=-1;
#define TIC() (tic=clock())
#define TOC() ((double)(clock()-tic)/CLOCKS_PER_SEC)
int main()
{
	double GHI=230.0, DHI=200.0, t;
	int N;
	sky_grid sky;
	sky_pos sun={degr2rad(20), degr2rad(180)};
	topology T;
	
	Verbosity=VERBOSE;
	//Connectivity(10);
	// return 0;
	TIC();
	sky=ssdp_init_sky(50);
	ssdp_make_perez_all_weather_sky(&sky, sun, GHI, DHI, 300);
	
	Print(NORMAL, "GHI:         %e\n", ssdp_total_sky_horizontal(sky,1));
	Print(NORMAL, "POA (sky):   %e\n", ssdp_total_sky_poa(sky, degr2rad(30), degr2rad(180),1));
	Print(NORMAL, "POA (total): %e\n", ssdp_total_poa(sky,0.25,degr2rad(30), degr2rad(180),1));
	
	
	T=LoadTopo("testtopo.dat");
	ssdp_mask_horizon(&sky,T,0.0,0.0,1);
	Print(NORMAL, "GHI:         %e\n", ssdp_total_sky_horizontal(sky,1));
	Print(NORMAL, "POA (sky):   %e\n", ssdp_total_sky_poa(sky, degr2rad(30), degr2rad(180),1));
	Print(NORMAL, "POA (total): %e\n", ssdp_total_poa(sky,0.25,degr2rad(30), degr2rad(180),1));
	
	WriteDome3D("3D_sky.dat", sky, 0, 1);
	WriteDome4D("4D_sky.dat", sky, 0, 1);
	t=TOC();
	Print(VERBOSE, "used %e s\n", t);
	ssdp_free_sky(&sky);
	ssdp_free_topology(&T);
	return 0;
}
