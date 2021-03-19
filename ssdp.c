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
#include "vector.h"
#include "sky_dome.h"
#include "sky_model.h"
#include "project.h"
#include "ground.h"
#include "io.h"
#include "util.h"
// tic-toc timer
clock_t tic=-1;
#define TIC() (tic=clock())
#define TOC() ((double)(clock()-tic)/CLOCKS_PER_SEC)
int main()
{
	double GHI=230.0, DHI=200.0, t;
	int i;
	sky_grid sky;
	sky_pos sun={degr2rad(20), degr2rad(180)};
	topology T;
	Verbosity=VERBOSE;
	//Connectivity(10);
	// return 0;
	TIC();
	sky=InitSky(50);
	PerezSky(&sky, sun, GHI, DHI, 300);
	Print(NORMAL, "GHI: %e\n", GlobalHorizontal(sky,1));
	Print(NORMAL, "POA: %e\n", PlaneOfArray(sky, degr2rad(30), degr2rad(180),1));
	// DumpSky(sky, 0, 1);
	T=LoadTopo("testtopo.dat");
	MakeHorizon(&sky, T, 1);
	Print(NORMAL, "GHI: %e\n", GlobalHorizontal(sky,1));
	Print(NORMAL, "POA: %e\n", PlaneOfArray(sky, degr2rad(30), degr2rad(180),1));
	Print(NORMAL, "POA ground reflected: %e\n", POA_Albedo(sky, 0.25, degr2rad(30), degr2rad(180),1));
	
	WriteDome3D("3D_sky.dat", sky, 0, 1);
	WriteDome4D("4D_sky.dat", sky, 0, 1);
	t=TOC();
	Print(VERBOSE, "used %e s\n", t);
	free_sky_grid(&sky);
	free_topo(&T);
	return 0;
}
