/*
    Simple Sky-Dome Projector Library
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
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include "util.h"
#include "parser.h"
#ifdef OPENMP
  #include <omp.h>
#endif


void Fatal( const char *format_str, ...)
{
	va_list ap;
	va_start (ap, format_str);
	vfprintf(stderr,format_str, ap); 
	exit(1);
}


void Warning( const char *format_str, ...)
{
	va_list ap;
	va_start (ap, format_str);
	
	if ((ParseLineNr>0)&&(ParseFileStr))
		fprintf(stderr, "WARNING: (%s:%d)\n", ParseFileStr,ParseLineNr);
	vfprintf(stderr,format_str, ap); 
}

static int _ProgressBar(int pcn, int pco, int len, int tics)
/* pcn: new percentage complete */
/* pco: old percentage complete */
/* len: length of the progress bar */
/* tics: number of tics */
{
	int pc_n=(pcn*len)/100;
	int pc=pco, tic;
	int i;
	pco=(pco*len)/100;
	if (pco==len)
		return pc;
	
	tic=len/tics;
	
	for (i=pco+1;i<=pc_n;i++)
	{
		if(i%tic == 0)
		{
			if (i==0)
				printf("|");
			else
			{
				printf("\b\b\b\b");
				printf("|");
			}
		}
		else
		{
			printf("\b\b\b\b");
			printf("=");
		}
		pc=pcn;
		printf("%3i%%",pc);
	}
	if (pcn>pc)
	{
		printf("\b\b\b\b");
		printf("%3i%%",pcn);
		pc=pcn;	
	}
		
	if (pc_n==len)
		printf("\n");	
	fflush(stdout);
	
	return pc;
}


void ProgressBar(int pcn, int* pco, int len, int tics)
{
#ifdef OPENMP
		if (omp_get_thread_num()==0)
				*pco = _ProgressBar(pcn, *pco, len, tics);
#else
		*pco = _ProgressBar(pcn, *pco, len, tics);
#endif
}
