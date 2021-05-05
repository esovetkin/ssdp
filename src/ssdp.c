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
#include <string.h>
#include <time.h>
#include "libssdp.h"
#include "util.h"
#include "parser.h"
#include "readlineshell.h"
#include "variables.h"

#define MAXLINELEN 4096
void StripEndline(char *line)
{
	int i;
	i=strlen(line)-1;
	while ((i>0)&&((line[i]=='\n')||(line[i]=='\r')))
	{
		line[i]='\0';
		i--;
	}
}
void ParseFile(char *fn)
{
	FILE *f;
	char *line;	
	if ((f=fopen(fn,"r"))==NULL)
	{
		
		Warning("Cannot open input file %s for reading", fn);
		return;
	}
	
	ParseLineNr=1;
	ParseFileStr=fn;
	
	line=malloc(MAXLINELEN*sizeof(char));
    fgets(line, MAXLINELEN-1, f);
	while(feof(f)==0)
	{
		StripEndline(line);
		if (ParseComm(line))
			break;
		fgets(line, MAXLINELEN-1, f);
		ParseLineNr++;
	}
	ParseFileStr=NULL;
	ParseLineNr=0;
	free(line);
	fclose(f);
}

void PrintHeader()
{
	printf("               |       \n");
	printf("  __|  __|  _` | __ \\  \n");
	printf("\\__ \\\\__ \\ (   | |   | \n");
	printf("____/____/\\__,_| .__/  \n");
	printf("                _|     \n");
	printf("linked libssdp: ");
	ssdp_print_version();
	printf("compile date  :  %s\n",__DATE__);
	printf("----------------------------\n");
}

int main(int argc, char **argv)
{
	int i=1;
	InitVars();
	if (argc>1)
	{
		if (strncmp("-q", argv[1], 3)==0)
			i++;
		else
			PrintHeader();
	}
	else 
		PrintHeader();
	if (i>=argc)
		shell();
	for (;i<argc;i++)
	{
		if (argv[i][0]=='-')
		{
			switch (argv[i][1])
			{
				case 'f':
					i++;
					if (i<argc)
						ParseFile(argv[i]);
					else
						fprintf(stderr,"Error: missing file after -f option\n");
					break;
				case 'i':
					shell();
					break;
				default:
					fprintf(stderr,"Error: unknown option %s\n", argv[i]);
			}
		}
		else
			if (ParseComm(argv[i]))
				break;
	}
	ClearVars();
	return 0;
}
/*
int main()
{
	double GHI=230.0, DHI=200.0, t;
	sky_grid sky;
	sky_pos sun={degr2rad(20), degr2rad(180)};
	topology T;
	
	TIC();
	sky=ssdp_init_sky(50);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		return 1;
	}
	ssdp_make_perez_all_weather_sky(&sky, sun, GHI, DHI, 300);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		return 1;
	}
	
	printf( "GHI:         %e\n", ssdp_total_sky_horizontal(&sky,1));
	printf( "POA (sky):   %e\n", ssdp_total_sky_poa(&sky, degr2rad(30), degr2rad(180),1));
	printf( "POA (total): %e\n", ssdp_total_poa(&sky,0.25,degr2rad(30), degr2rad(180),1));
	
	//T=ssdp_make_rand_topology(100,100,10,0.01, 12000);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_sky(&sky);
		return 1;
	}
	//WriteTopo("testtopo.dat", T);
	T=LoadTopo("testtopo.dat");
	//T=LoadTopo("simpletopo.dat");
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_sky(&sky);
		ssdp_free_topology(&T);
		return 1;
	}
	ssdp_mask_horizon_z_to_ground(&sky,&T,0.0,0.0,0.1, NULL);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_sky(&sky);
		ssdp_free_topology(&T);
		return 1;
	}
	printf( "GHI:         %e\n", ssdp_total_sky_horizontal(&sky,1));
	printf( "POA (sky):   %e\n", ssdp_total_sky_poa(&sky, degr2rad(30), degr2rad(180),1));
	printf( "POA (total): %e\n", ssdp_total_poa(&sky,0.25,degr2rad(30), degr2rad(180),1));
	printf( "POA (total): %e\n", ssdp_total_poa_effective(&sky,0.25,degr2rad(30), degr2rad(180),1.5,1));
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_sky(&sky);
		ssdp_free_topology(&T);
		return 1;
	}
	
	WriteDome3D("3D_sky.dat", &sky, 0, 1);
	WriteDome4D("4D_sky.dat", &sky, 0, 1);
	t=TOC();
	printf( "used %e s\n", t);
	ssdp_unmask_horizon(&sky);	
	
	TIC();
	RasterPOA("POAraster.dat", &sky, &T, 0.25, 0.3, degr2rad(180), degr2rad(0), -40, -40, 40, 40, 100, 100);
	t=TOC();
	printf( "used %e s\n", t);
	TIC();
	RasterTopology("TOPOraster.dat", &T, -40, -40, 40, 40, 100, 100);
	t=TOC();
	printf( "used %e s\n", t);
	ssdp_free_sky(&sky);
	ssdp_free_topology(&T);
	return 0;
}*/
