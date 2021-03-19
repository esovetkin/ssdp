#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include "util.h"
VERB Verbosity=NORMAL;



void Print(VERB v, const char *format_str, ...)
{
	va_list ap;
	va_start (ap, format_str);
	if (v<=Verbosity)
		vprintf(format_str, ap); 
}
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
	vfprintf(stderr,format_str, ap); 
}

int ProgressBar(int pcn, int pco, int len, int tics)
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
