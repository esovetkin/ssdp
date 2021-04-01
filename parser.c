#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
/* local includes */
#include "parsedef.h"
#include "libssdp.h"
#include "io.h"
#include "util.h"
#include "variables.h"

#define rad2degr(rad) ((rad)*180/M_PI)
#define degr2rad(degr) ((degr)*M_PI/180)

char * GetWord(char *in, char *word)
/* take string in and allocated string word
 * copy the first word of strin in to word
 * return a pointer to the rest of the string
 */
{
	char *end=in;
	char *out=word;
	int squoted=0, dquoted=0;
	int go=1;
	while (*end && go)
	{
		if (*end=='\'')
			squoted=!squoted;
		if (*end=='\"')
			dquoted=!dquoted;
			
		if (!(squoted||dquoted)&&(isspace(*end)))
		{
			go=0;
		}
		else 
		{
			if ((*end!='\'')||(*end=='\"'))
			{
				*out=*end;
				out++;
			}
		}
		end++;
	}
	*out='\0';
	if (squoted)
		Warning("unmatched single quotes\n");
	if (dquoted)
		Warning("unmatched double quotes\n");
	go=1; /* skip blank chars */
	while (*end && go)
	{
		if (!isspace(*end))
		{
			go=0;
		}
		else
			end++;
	}
	return end;
}

ParserFun LookupComm(char *key)
{
	int i=0;
	int lk;
	lk=strlen(key);
	while (KeyTable[i].key)
	{
		if (strlen(KeyTable[i].key)==lk)
			if (strncmp(KeyTable[i].key, key, lk)==0)
				return (KeyTable[i].fun);
		i++;
	}
	return NULL;
}


int ParseComm(char *in)
{
	char *key;
	char *arg;
	int go=1;
	ParserFun fun;
	/* skip space chars */
	if (*in=='#')
		return 0;
	while (*in && go)
	{
		if (!isspace(*in))
		{
			go=0;
		}
		else
			in++;
	}
	if (!(*in))
		return 0;
	key=malloc((strlen(in)+1)*sizeof(char));
	arg=GetWord(in, key);
	if (strncmp(key, "exit", 5)==0)
	{
		free(key);
		return 1;
	}
		
	fun=LookupComm(key);
	if (fun)
		fun(arg);
	else
		Warning("Command %s not defined\n", key);	
	free(key);
	return 0;
}

int GetOption(char *in, char *opt, char *word)
{
	char *start;
	char *opti;
	int len;
	len=(strlen(opt)+2);
	opti=malloc(len*sizeof(char));
	snprintf(opti,len,"%s=",opt);
	start=strstr(in, opti);
	
	if (!start)
	{
		*word='\0';
		free(opti);	
		return 0;
	}	
	GetWord(start+len-1, word);
	free(opti);	
	return 1;
}

int GetArg(char *in, char *opt, char *word)
{
	if (!GetOption(in, opt, word))
	{		
		Warning("Mandatory argument %s missing\n", opt);	
		return 0;
	}
	return 1;
}



/* parse routines
 * for each parse routine add a parseflag:
 * '// PARSEFLAG <keyword> <function-name>'
 * Note: you can create aliases by adding several lines like this with 
 * different keywords but the same function.
 * a parse routine must be of type void and take one character string as input
 * The gen_parseflags will take care of all the definitions in parsedef.h
 */
 
// PARSEFLAG load_topography TopoLoad "T=<topology-variable> file=<filename>"
void TopoLoad(char *in)
{
	topology T;
	char *name; 
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "T", name);
	if (GetArg(in, "file", file))
		T=LoadTopo(file);
	free(file);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_topology(&T);
		free(name);
		return;
	}
	else
	{
		printf("Defining topology \"%s\"\n", name);
		if (AddTopology(name, T))
		{
			ssdp_free_topology(&T);
			free(name);
		}
	}	
}

// PARSEFLAG save_topography TopoSave "T=<topology-variable> file=<filename>"
void TopoSave(char *in)
{
	topology *T;
	char *name;
	char *file; 
	
	name=malloc((strlen(in)+1)*sizeof(char));
	file=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "T", name);
	GetArg(in, "file", file);
	if (!LookupTopology(name, &T))
	{
		Warning("Topology %s not available\n",name); 
		free(name);
		free(file);
		return;
	}
	else
		WriteTopo(file, T);
	free(name);
	free(file);
}
// PARSEFLAG make_vector MakeVec "v=<vector-variable> z=<zenith> a=<azimuth>"
void MakeVec(char *in)
{
	sky_pos v;
	char *name;
	
	name=malloc((strlen(in)+1)*sizeof(char));
	if (GetArg(in, "z", name))
		v.z=degr2rad(atof(name));
	else
	{
		free(name);
		return;
	}
	if (GetArg(in, "a", name))
		v.a=degr2rad(atof(name));
	else
	{
		free(name);
		return;
	}
	
	if (GetArg(in, "v", name))
	{
		if (AddVec(name, v))
			free(name);
		return;
	}
	free(name);
}
// PARSEFLAG show_vector ShowVec "v=<vector-variable>"
void ShowVec(char *in)
{
	sky_pos *v;
	char *name;
	
	name=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "v", name);
	if (!LookupVec(name, &v))
	{
		Warning("Vector %s not available\n",name); 
		free(name);
		return;
	}
	printf("vector:  %s\nazimuth: %e\nzenith:  %f\n", name, rad2degr(v->a), rad2degr(v->z));
	free(name);
}
// PARSEFLAG solar_position SunPosition "v=<vector-variable> date=<jjjj-mm-dd-hh:mm:ss> lon=<longitude> lat=<latitude>"
void SunPosition(char *in)
{
	sky_pos v;
	double lat, lon;
	char *name;
	char *date;
	char *p;
	struct tm T={0,0,0,1,1,2021,5,0,0};
	
	name=malloc((strlen(in)+1)*sizeof(char));
	date=malloc((strlen(in)+1)*sizeof(char));
	if (GetArg(in, "date", date))
	{
		p=strptime(date, "%Y-%m-%d-%H:%M:%S", &T); 
		if (!p)
		{
			Warning("Could not process date string %s\n", date); 
			free(name);
			free(date);
			return;
		}
		if (*p)
		{
			Warning("Error processing date string, lost you here %s\n", p); 
			free(name);
			free(date);
			return;
		}
		// success processing date string, continute
		free(date); // no longer needed
		if (GetArg(in, "lat", name))
			lat=atof(name);
		else
		{
			Warning("No latitude provided to SunPos\n"); 
			free(name);
			return;
		}
		if (GetArg(in, "lon", name))
			lon=atof(name);
		else
		{
			Warning("No longitude provided to SunPos\n"); 
			free(name);
			return;
		}
		
		if (GetArg(in, "v", name))
		{
			v=ssdp_sunpos(mktime(&T), lat, lon);	
			printf("Defining solar position in \"%s\"\n", name);
			if (AddVec(name, v))
				free(name);	
			return;
		}
	}
	else
	{
		Warning("No date string provided to SunPos\n"); 
		free(name);
		free(date);
		return;
	}
}

// PARSEFLAG init_sky InitSky "S=<sky-dome-variable> N=<Nz>"
void InitSky(char *in)
{
	sky_grid S;
	int Nz=-1;
	char *name;
	char *nz;
	
	name=malloc((strlen(in)+1)*sizeof(char));
	nz=malloc((strlen(in)+1)*sizeof(char));
	GetArg(in, "S", name);
	if (GetArg(in, "N", nz))
		Nz=atoi(nz);
	free(nz);
		
	if (Nz>1)
	{
		S=ssdp_init_sky(Nz);
		if (ssdp_error_state)
		{
			ssdp_print_error_messages();
			ssdp_free_sky(&S);
			free(name);
			return;
		}
		else
		{
			printf("Defining sky \"%s\"\n", name);
			if (AddSkyDome(name, S))
			{
				ssdp_free_sky(&S);
				free(name);
			}
			return;
		}
	}
	free(name);
	
}

// PARSEFLAG perez_fill_sky PerezSkyCoord "S=<sky-dome-variable> GHI=<GHI> DHI=<DHI> date=<jjjj-mm-dd-hh:mm:ss> lon=<longitude> lat=<latitude>"
void PerezSkyCoord(char *in)
{
	sky_grid *S;
	double GHI, DHI;
	double lat, lon;
	char *word;
	char *p;
	struct tm T={0,0,0,1,1,2021,5,0,0};
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (GetArg(in, "S", word))
	{
		if (!LookupSkyDome(word, &S))
		{
			Warning("Sky Dome %s not available\n",word); 
			free(word);
			return;
		}
	}
	else
	{
		free(word);
		return;
	}
	if (GetArg(in, "GHI", word))
	{
		GHI=atof(word);
		if (GHI<0)
		{
			Warning("Global Horizontal Irradiance cannot be negative\n",word); 
			free(word);
			return;
		}
	}
	else
	{
		free(word);
		return;
	}
	if (GetArg(in, "DHI", word))
	{
		DHI=atof(word);
		if (DHI<0)
		{
			Warning("Diffuse Horizontal Irradiance cannot be negative\n",word); 
			free(word);
			return;
		}
	}
	else
	{
		free(word);
		return;
	}	
	if (GetArg(in, "date", word))
	{
		p=strptime(word, "%Y-%m-%d-%H:%M:%S", &T); 
		if (!p)
		{
			Warning("Could not process date string %s\n", word); 
			free(word);
			return;
		}
		if (*p)
		{
			Warning("Error processing date string, lost you here %s\n", p); 
			free(word);
			return;
		}
		// success processing date string, continute
		if (GetArg(in, "lat", word))
			lat=atof(word);
		else
		{
			Warning("No latitude provided to perez_sky_coordinate\n"); 
			free(word);
			return;
		}
		if (GetArg(in, "lon", word))
			lon=atof(word);
		else
		{
			Warning("No longitude provided to perez_sky_coordinate\n"); 
			free(word);
			return;
		}
		free(word);
	}
	else
	{
		free(word);
		return;
	}
	ssdp_make_perez_all_weather_sky_coordinate(S, mktime(&T), degr2rad(lon), degr2rad(lat), GHI, DHI);
}
// PARSEFLAG perez_fill_sky_explicit PerezSky "S=<sky-dome-variable> v=<sun-position-vector> GHI=<GHI> DHI=<DHI> yday=<day-of-year>"
void PerezSky(char *in)
{
	sky_grid *S;
	sky_pos *v;
	double GHI, DHI;
	int yday;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (GetArg(in, "S", word))
	{
		if (!LookupSkyDome(word, &S))
		{
			Warning("Sky Dome %s not available\n",word); 
			free(word);
			return;
		}
	}
	if (GetArg(in, "v", word))
	{
		if (!LookupVec(word, &v))
		{
			Warning("Vector %s not available\n",word); 
			free(word);
			return;
		}
	}
	if (GetArg(in, "GHI", word))
	{
		GHI=atof(word);
		if (GHI<0)
		{
			Warning("Global Horizontal Irradiance cannot be negative\n",word); 
			free(word);
			return;
		}
	}
	if (GetArg(in, "DHI", word))
	{
		DHI=atof(word);
		if (DHI<0)
		{
			Warning("Diffuse Horizontal Irradiance cannot be negative\n",word); 
			free(word);
			return;
		}
	}
	if (GetArg(in, "yday", word))
	{
		yday=atoi(word);
		if (yday<0)
		{
			Warning("Day of the year cannot be negative\n",word); 
			free(word);
			return;
		}
		if (yday>366)
		{
			Warning("Day of the year cannot be larger than 366\n",word); 
			free(word);
			return;
		}
	}
	ssdp_make_perez_all_weather_sky(S, *v, GHI, DHI, yday);
}


void SplitWords(char *in, char *ident)
{
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	while(*in)
	{
		in=GetWord(in, word);
		printf("%s%s\n", ident, word);
	}
	free(word);
}
// PARSEFLAG help Help "[-l/command]"
void Help(char *in)
{
	int i=0;
	int lk;
	lk=strlen(in);
	if ((strncmp(in,"-l", 3)==0)||(!*in))
	{
		while(Usage[i])
		{
			printf("Command %s:\n",KeyTable[i].key);
			SplitWords(Usage[i], "\t");
			printf("\n");
			i++;
		}
		return;
	}
			
	while (KeyTable[i].key)
	{
		if (strlen(KeyTable[i].key)==lk)
			if (strncmp(KeyTable[i].key, in, lk)==0)
			{
				printf("Command %s:\n",KeyTable[i].key);
				SplitWords(Usage[i], "\t");
				printf("\n");
				return;
			}
		i++;
	}
	printf("Command %s not known\n",in);
}


char * keyword_generator(const char *text, int state)
{
    static int list_index, len;
    char *name;

    if (!state) {
        list_index = 0;
        len = strlen(text);
    }

    while ((name = KeyTable[list_index++].key)) {
        if (strncmp(name, text, len) == 0) {
            return strdup(name);
        }
    }

    return NULL;
}


