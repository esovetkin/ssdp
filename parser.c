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

char * GetWord(const char *in, char *word);
/* core parsing routines */
/* LookupComm finds takes a keyword and retuns a pointer to the corresponding parser routine */
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

/* ParseComm takes a line of input and finds the keyword at the beginning
 * then calls LookupComm to find the right parser, then calls this parser with the remainder of the line */
int ParseComm(const char *in)
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

/* string manimulation routines (mostly about getting the right chunk out of a string) */
char * GetWord(const char *in, char *word)
/* take string in and allocated string word
 * copy the first word of strin in to word
 * return a pointer to the rest of the string
 */
{
	char *end=(char *)in;
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


int GetOption(const char *in, const char *opt, char *word)
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

int GetArg(const char *in, const char *opt, char *word)
{
	if (!GetOption(in, opt, word))
	{		
		Warning("Mandatory argument %s missing\n", opt);	
		return 0;
	}
	return 1;
}


void SplitWords(const char *in, const char *ident)
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

/* readline command completion */
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

/* parse routines
 * for each parse routine add a parseflag:
 * '// PARSEFLAG <keyword> <function-name>'
 * Note: you can create aliases by adding several lines like this with 
 * different keywords but the same function.
 * a parse routine must be of type void and take one character string as input
 * The gen_parseflags will take care of all the definitions in parsedef.h
 */
int FetchConfig(const char *in, const char *pat, char *str, simulation_config **a)
{
	if (GetArg(in, pat, str))
	{
		if (!LookupSimConf(str, a))
		{
			Warning("Simulation config %s is not available\n",str);
			return 1;
		}
		return 0;
	}
	return 1;
}
int FetchArray(const char *in, const char *pat, char *str, array **a)
{
	if (GetArg(in, pat, str))
	{
		if (!LookupArray(str, a))
		{
			Warning("array %s is not available\n",pat);
			return 1;
		}
		return 0;
	}
	return 1;
}
int FetchFloat(const char *in, const char *pat, char *str, double *a)
{
	if (GetArg(in, pat, str))
	{
		(*a)=atof(str);
		return 0;
	}
	return 1;
}
int FetchInt(const char *in, const char *pat, char *str, int *a)
{
	if (GetArg(in, pat, str))
	{
		(*a)=atoi(str);
		return 0;
	}
	return 1;
}


// PARSEFLAG init_sim_config InitConfig "C=<config-variable>"
void InitConfig(char *in)
{
	simulation_config C;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (GetArg(in, "C", word))
	{
		C=InitConf();
		printf("Defining simulation configuration \"%s\"\n", word);
		if (AddSimConf(word, C)) // note InitConf does not allocate anything, no need to free
			free(word);
		return;
	}
	free(word);
}
// PARSEFLAG config_aoi ConfigAOI "C=<config-variable> model=<none/front-cover/anti-reflect/user> [nf=<front-cover-refractive-index> [nar=<antireflection-refractive-index>]] [file=<user-defined-aoi>]"
void ConfigAOI(char *in)
{
	simulation_config *C;
	AOI_Model M;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
		
	if (!GetArg(in, "model", word))
	{
		free(word);
		return;
	}
	if (strncmp(word,"none",12)==0)
		M=AOI_NONE;
	else if (strncmp(word,"front-cover",12)==0)
		M=AOI_GLASS;
	else if (strncmp(word,"anti-reflect",12)==0)
		M=AOI_GLASS_AR;
	else if (strncmp(word,"user",12)==0)
		M=AOI_USER;
	else
	{
		Warning("Unknown AOI model %s\n",word); 
		free(word);
		return;
	}
	
	switch (M)
	{
		case AOI_GLASS:
		{
			if (FetchFloat(in, "nf", word, &(C->M.ng)))
			{
				free(word);
				return;
			}
			printf("Setting AOI model to \"front-cover\" with a refractive index of %e\n", C->M.ng);
			C->M.M=M;
			return;
		}
		case AOI_GLASS_AR:
		{
			if (FetchFloat(in, "nf", word, &(C->M.ng)))
			{
				free(word);
				return;
			}
			if (FetchFloat(in, "nar", word, &(C->M.nar)))
			{
				free(word);
				return;
			}
			printf("Setting AOI model to \"anti-reflect\" with a cover refractive index of %e\n", C->M.ng);
			printf("and a anti-reflecxtion refrective index of %e\n", C->M.nar);
			C->M.M=M;
			return;
		}
		case AOI_USER:
		{
			if (!GetArg(in, "file", word))
			{
				free(word);
				return;
			}
			// read and allocate theta-effT data
			Warning("Not implemented yet\n"); 
			C->M.M=M;
			return;
		}
		default:
			break;		
	}
	free(word);
}

void FreeConfigMask(simulation_config *C)
{
	int i;
	if (C->mask)
	{
		for (i=0;C->Nl;i++)
			ssdp_free_sky_mask(C->mask+i);
	}
	free(C->mask);
	C->mask=NULL;
}

void FreeConfigLocation(simulation_config *C)
{
	if (C->loc_init)
	{
		if (C->x)
			free(C->x);
		if (C->y)
			free(C->y);
		if (C->z)
			free(C->z);
		if (C->o)
			free(C->o);
	}
	FreeConfigMask(C);
	C->Nl=0;
	C->loc_init=0;
}

void InitConfigMask(simulation_config *C)
{
	int i;
	if ((C->sky_init)&&(C->topo_init)&&(C->loc_init))
	{
		FreeConfigMask(C); // make sure we are clear to allocate new memory
		C->mask=malloc(C->Nl*sizeof(sky_mask));
		for (i=0;i<C->Nl;i++)
			C->mask[i]=ssdp_mask_horizon(&(C->S),&(C->T),C->x[i],C->y[i],C->z[i]);
	}
}

// PARSEFLAG config_sky ConfigSKY "C=<config-variable> N=<zenith-steps>"
void ConfigSKY(char *in)
{
	simulation_config *C;
	char *word;
	int N;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
		
	if (FetchInt(in, "N", word, &N))
	{
		free(word);
		return;
	}
	if (N>1)
	{
		if (C->sky_init)
			ssdp_free_sky(&C->S);
		else
			C->sky_init=1;
		C->S=ssdp_init_sky(N);
		InitConfigMask(C);
		printf("Configuring sky with %d zenith discretizations\n", N);
	}
	else
		Warning("Number of zenith discretizations must me larger than 1\n");
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_sky(&C->S);
		C->sky_init=0;
		ssdp_reset_errors();
	}
	free(word);
	return;
}

// PARSEFLAG config_topology ConfigTOPO "C=<config-variable> x=<x-array-variable> y=<y-array-variable> z=<z-array-variable>"
void ConfigTOPO (char *in)
{
	simulation_config *C;
	array *x, *y, *z;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
		
	if (FetchArray(in, "x", word, &x))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "y", word, &y))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "z", word, &z))
	{
		free(word);
		return;
	}
	free(word);
	if ((x->N!=y->N)||(x->N!=z->N))
	{
		Warning("x- y- and z-arrays must be of equal length\n"); 
		return;
	}
	
	if (C->topo_init)
	{
		ssdp_free_topology(&C->T);
	}
	else
		C->topo_init=1;
	printf("Configuring topology with %d points\n", x->N);
	C->T=ssdp_make_topology(x->D, y->D, z->D, x->N);
	if (ssdp_error_state)
	{
		ssdp_print_error_messages();
		ssdp_free_topology(&C->T);
		C->topo_init=0;
		ssdp_reset_errors();
	}
	InitConfigMask(C);
	return;
}
// PARSEFLAG config_locations ConfigLoc "C=<config-variable> x=<x-array-variable> y=<y-array-variable> z=<z-array-variable> azimuth=<azimuth-array-variable> zenith=zenith-array-variable>"
void ConfigLoc (char *in)
{
	simulation_config *C;
	array *x, *y, *z, *az, *ze;
	int i;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
		
	if (FetchArray(in, "x", word, &x))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "y", word, &y))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "z", word, &z))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "azimuth", word, &az))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "zenith", word, &ze))
	{
		free(word);
		return;
	}
	free(word);
	if (x->N!=y->N)
	{
		Warning("x- y- arrays must be of equal length\n"); 
		return;
	}
	if ((z->N!=x->N)&&(z->N!=1))
	{
		Warning("z- array must be either of equal length as x- and y- or of length 1\n"); 
		return;
	}
	if (az->N!=ze->N)
	{
		Warning("azimuth and zenith arrays must be of equal length\n"); 
		return;
	}
	if ((az->N!=x->N)&&(az->N!=1))
	{
		Warning("azimuth and zenith arrays must either have length 1 or be of equal length ans z,y,z\n"); 
		return;
	}
	
	FreeConfigLocation(C);
	printf("Configuring locations with %d points\n", x->N);
	
	C->Nl=x->N;
	C->x=malloc(C->Nl*sizeof(double));
	C->y=malloc(C->Nl*sizeof(double));
	C->z=malloc(C->Nl*sizeof(double));
	C->o=malloc(C->Nl*sizeof(sky_pos));
	
	if ((z->N==1)&&(az->N==1))
	{
		for (i=0;i<C->Nl;i++)
		{
			C->x[i]=x->D[i];
			C->y[i]=y->D[i];
			C->z[i]=z->D[0];
			C->o[i].a=az->D[0];
			C->o[i].z=ze->D[0];
		}
	}
	else if (z->N==1) // many orientations, one z
	{
		for (i=0;i<C->Nl;i++)
		{
			C->x[i]=x->D[i];
			C->y[i]=y->D[i];
			C->z[i]=z->D[0];
			C->o[i].a=az->D[i];
			C->o[i].z=ze->D[i];
		}
	}
	else if (az->N==1) // many z's one orientation
	{
		for (i=0;i<C->Nl;i++)
		{
			C->x[i]=x->D[i];
			C->y[i]=y->D[i];
			C->z[i]=z->D[i];
			C->o[i].a=az->D[0];
			C->o[i].z=ze->D[0];
		}
	}
	else // just many
	{
		for (i=0;i<C->Nl;i++)
		{
			C->x[i]=x->D[i];
			C->y[i]=y->D[i];
			C->z[i]=z->D[i];
			C->o[i].a=az->D[i];
			C->o[i].z=ze->D[i];
		}
	}
	C->loc_init=1;
	
	InitConfigMask(C);
	
	return;
}


/* add some parsers:
 * 1: compute z by specifying height and sample the topology for a grid
 * 2: compute orientation from direction
 * 3: modify orientation by surface normal
 * 4: make ground albedo editable
 * 4: if all works, can we save and load a config?
 */ 


// PARSEFLAG sim_static SimStatic "C=<config-variable> t=<array-variable> GHI=<array-variable> DHI=<array-variable> POA=<out-array>"
void SimStatic(char *in)
{
	int i, j; // loop through space and time
	char *word;
	simulation_config *C;
	array *t, *GH, *DH, out;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}		
	if (!C->sky_init)
	{		
		Warning("Simulation config has no sky initialized\n");
		free(word);
		return;
	}	
	if (!C->topo_init) // TODO: in this case just omit the horizon and compute only one location?
	{	
		Warning("Simulation config has no topology initialized\n");
		free(word);
		return;
	}
	if (!C->loc_init) // TODO: in this case just omit the horizon and compute only one location?
	{	
		Warning("Simulation config has no locations initialized\n");
		free(word);
		return;
	}		
	if (FetchArray(in, "t", word, &t))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "GHI", word, &GH))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "DHI", word, &DH))
	{
		free(word);
		return;
	}	
	if ((t->N!=GH->N)||(t->N!=DH->N))
	{
		Warning("Length of t-, GHI-, and DHI-arrays do not match\n");
		free(word);
		return;
	}
	// fetch name of output var
	if (!GetArg(in, "POA", word))
	{
		free(word);
		return;
	}
	out.D=malloc(t->N*C->Nl*sizeof(double));
	if (out.D==NULL)
	{
		free(word);
		return;
	}	
	out.N=t->N*C->Nl;	
	for (j=0;j<t->N;j++)
	{
		// compute sky at evert time instance
		ssdp_make_perez_all_weather_sky_coordinate(&(C->S), (time_t) t->D[j], C->lon, C->lat, GH->D[j], DH->D[j]);
		for (i=0;i<C->Nl;i++)
		{
			// compute POA at evert location  instance
			out.D[j*C->Nl+i]=ssdp_total_poa(&(C->S), C->albedo, C->o[i], &(C->M),C->mask+i);
		}
	}
	if(AddArray(word, out))
	{
		free(word); // failed to make array
		free(out.D);
	}	
}
// _PARSEFLAG sim_route SimRoute "C=<config-variable> x=<array-variable> y=<array-variable> t=<array-variable> GHI=<array-variable> DHI=<array-variable> POA=<out-array>"
/*void SimStatic(char *in)
{
	sky_mask *mask; // a mask for every grid point
	int i, j; // loop through space and time
	
	
}*/


int GetNumOption(char *in, char *opt, int i, char *word)
{
	char *start;
	char *opti;
	int len, k;
	len=1;
	if (i<0)
		len++;
	k=abs(i);
	while((k=k/10)>0)
		len++;
		
	len+=(strlen(opt)+2);
	opti=malloc(len*sizeof(char));
	snprintf(opti,len,"%s%d=",opt,i);
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
// PARSEFLAG read_array ReadArraysFromFile "a0=<array0> a1=<array1> .. aN=<arrayN> file=<file>"
void ReadArraysFromFile(char *in)
{
	char **names;
	char *word;
	char *file;
	double **data;
	int i, k, Na=4, N;
	file=malloc((strlen(in)+1)*sizeof(char));
	if (!GetArg(in, "file", file))
	{
		free(file);
		return;
	}
	word=malloc((strlen(in)+1)*sizeof(char));
	i=0;
	names=malloc(Na*sizeof(char *));
	while(GetNumOption(in, "a", i, word))
	{
		names[i]=word;
		i++;
		word=malloc((strlen(in)+1)*sizeof(char));
		if (i==Na-1)
		{
			Na+=4;
			names=realloc(names, Na*sizeof(char *));
		}
	}
	free(word);
	if (i==0)
	{
		Warning("Cannot define arrays from file, no array arguments recognized\n"); 
		free(names);
		free(file);
		return;
	}
	data=ReadArrays(file, i, &N);
	free(file);
	if (N>0)
	{
		array a;
		a.N=N;
		// create i arrays
		for (k=0;k<i;k++)
		{
			a.D=data[k];
			if(AddArray(names[k], a))
			{
				free(names[k]); // failed to make array
				free(a.D);
			}
		}
	}
	else
	{
		for (k=0;k<i;k++)
		{
			free(data[k]);
			free(names[k]);
		}
	}
	free(names);
	free(data);
}
// PARSEFLAG write_array WriteArraysToFile "a0=<array0> a1=<array1> .. aN=<arrayN> file=<file>"
void WriteArraysToFile(char *in)
{
	char *word;
	char *file;
	double **data;
	array *a;
	int i, Na=4, N=-1;
	file=malloc((strlen(in)+1)*sizeof(char));
	if (!GetArg(in, "file", file))
	{
		free(file);
		return;
	}
	word=malloc((strlen(in)+1)*sizeof(char));
	i=0;
	data=malloc(Na*sizeof(double *));
	while(GetNumOption(in, "a", i, word))
	{
		
		if (!LookupArray(word, &a))
		{
			Warning("Array %s not available\n",word); 
			free(word);
			free(data);
			free(file);
			return;
		}
		data[i]=a->D;
		if (N<0)
			N=a->N;
		else if (a->N!=N)
		{
			Warning("Error: arrays must be of equal length\n"); 
			free(word);
			free(data);
			free(file);
			return;
		}
		i++;
		if (i==Na-1)
		{
			Na+=4;
			data=realloc(data, Na*sizeof(double *));
		}
	}
	free(word);
	if (i==0)
	{
		Warning("Cannot define arrays from file, no array arguments recognized\n"); 
		free(data);
		free(file);
		return;
	}
	WriteArrays(file,data,i,N);
	free(file);
	free(data);
}
// PARSEFLAG make_array MakeArray "x=<array-variable> x1=<startval> x2=<endval> Nx=<x-steps>"
void MakeArray(char *in)
{
	char *word;
	double x1, x2, dx;
	array a;
	int i, Nx;
	word=malloc((strlen(in)+1)*sizeof(char));
		
	if (FetchFloat(in, "x1", word, &x1))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "x2", word, &x2))
	{
		free(word);
		return;
	}
	if (FetchInt(in, "Nx", word, &Nx))
	{
		free(word);
		return;
	}
	Nx=abs(Nx);
	a.D=malloc((Nx+1)*sizeof(double));
	if (!a.D)
	{
		Warning("memory allocation failed\n");
		free(word);
		return;
	}
	a.N=(Nx+1);
	if (Nx==0)
		dx=(x2-x1);
	else
		dx=(x2-x1)/((double)Nx);
	for (i=0;i<=Nx;i++)
		a.D[i]=x1+((double)i)*dx;
	
	if (!GetArg(in, "x", word))
	{
		free(word);
		return;
	}	
	printf("Creating array %s\n",word);
	if(AddArray(word, a))
	{
		Warning("Failed to create array %s\n",word);
		free(a.D);	
		free(word);
	}
	return;	
}

// PARSEFLAG make_grid MakeGrid "x=<array-variable> x=<array-variable> x1=<start-x> x2=<end-x> y1=<start-y> y2=<end-y> Nx=<x-stapes> Ny=<y-steps>"
void MakeGrid(char *in)
{
	char *word;
	double x1, x2, y1, y2, dx, dy;
	array a, b;
	int i, j, k, Nx, Ny;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchFloat(in, "x1", word, &x1))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "x2", word, &x2))
	{
		free(word);
		return;
	}
	if (FetchInt(in, "Nx", word, &Nx))
	{
		free(word);
		return;
	}
	Nx=abs(Nx);
	if (FetchFloat(in, "y1", word, &y1))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "y2", word, &y2))
	{
		free(word);
		return;
	}
	if (FetchInt(in, "Ny", word, &Ny))
	{
		free(word);
		return;
	}
	Ny=abs(Ny);
	
	a.D=malloc((Nx+1)*(Ny+1)*sizeof(double));
	if (!a.D)
	{
		Warning("memory allocation failed\n");
		free(word);
		return;
	}
	a.N=(Nx+1)*(Ny+1);
	b.D=malloc(a.N*sizeof(double));
	if (!b.D)
	{
		Warning("memory allocation failed\n");
		free(word);
		return;
	}
	b.N=a.N;
	if (Nx==0)
		dx=(x2-x1);
	else
		dx=(x2-x1)/((double)Nx);
	if (Ny==0)
		dy=(y2-y1);
	else
		dy=(y2-y1)/((double)Ny);
	k=0;
	for (i=0;i<=Nx;i++)
		for (j=0;j<=Ny;j++)
		{
			a.D[k]=x1+((double)i)*dx;
			b.D[k]=y1+((double)j)*dy;
			k++;
		}
	
	if (!GetArg(in, "x", word))
	{
		free(word);
		return;
	}	
	if(AddArray(word, a))
	{
		free(a.D);	
		free(word);
	}
	word=malloc((strlen(in)+1)*sizeof(char)); // allocate new, word is swallowed into the variable list by AddArray
	if (!GetArg(in, "y", word))
	{
		free(word);
		return;
	}	
	if(AddArray(word, b))
	{
		free(b.D);
		free(word);
	}
	return;	
}


// _PARSEFLAG sim_static SimStatic "x=<array-variable> y=<array-variable> z=<height-value> t=<array-variable> GHI=<array-variable> DHI=<array-variable> O=<orientation-vector> C=<config-variable> POA=<out-array>"
/* let us add everything to the config variable, it takes no space 
 * how about the spatial verus temporal coordinates?
 * runtime check whether x,y,GHI,DHI,t all have the same length, or
 * GHI,DHI,t have the same length and x,y not?
 * Maybe no check, just let for x rund and t run and they end the same time, no interpolation, just nearest
 */

/*void SimStatic(char *in)
{
	sky_mask *mask; // a mask for every grid point
	int i, j; // loop through space and time
	
	
}*/
// _PARSEFLAG sim_route  SimRoute  "x=<array-variable> y=<array-variable> z=<height-value> t=<array-variable> GHI=<array-variable> DHI=<array-variable> O=<orientation-vector> C=<config-variable> POA=<out-array>"



/* how to get a useful simulator command set:
 * need to be able to specify:
 * GHI,DHI vs time
 * x,y vs time
 * x,y arrays
 * sim_static <(x,y)> <(GHI,DHI)(t)>
 * 		compute horizon for every (x,y)
 * 		compute sky S(t)
 * 		compute poa for every (x,y,t)
 * 
 * sim_route <(x,y)(t)> <(GHI,DHI)(t)>
 *  	compute horizon for every (x,y)
 * 		compute sky S(t)
 * 		compute poa for every (x(t),y(t),t)
 */		


/*
/ PARSEFLAG project_diffuse_sky_poa ProjectDiffSkyPOA "C=<config-variable> S=<sky> pn=<panel-surface-normal>"
/ PARSEFLAG project_direct_sky_poa ProjectDirSkyPOA "C=<config-variable> S=<sky> pn=<panel-surface-normal"
/ PARSEFLAG project_ground_albedo_poa ProjectGroundAlbedoPOA "C=<config-variable> S=<sky> pn=<panel-surface-normal"
/ PARSEFLAG project_total_poa ProjectTotalPOA "C=<config-variable> S=<sky> pn=<panel-surface-normal"
/ PARSEFLAG project_diffuse_sky_h ProjectDiffSkyH "C=<config-variable> S=<sky>"
/ PARSEFLAG project_direct_sky_h ProjectDirSkyH "C=<config-variable> S=<sky>"\
/ PARSEFLAG project_total_h ProjectTotalH"C=<config-variable> S=<sky>"
void ProjectDiffSkyPOA(char *in)
{
	simulation_config *C;
	sky_grid *S;
	sky_pos *v;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (GetArg(in, "C", word))
	{
		if (!LookupSimConf(word, &C))
		{
			Warning("Simulation config %s not available\n",word); 
			free(word);
			return;
		}
	}
	else
	{
		free(word);
		return;
	}
	if (GetArg(in, "S", word))
	{
		if (!LookupSkyDome(word, &S))
		{
			Warning("Sky-Dome %s not available\n",word); 
			free(word);
			return;
		}
	}
	else
	{
		free(word);
		return;
	}
	if (GetArg(in, "pn", word))
	{
		if (!LookupVec(word, &v))
		{
			Warning("Vector %s not available\n",word); 
			free(word);
			return;
		}
	}
	else
	{
		free(word);
		return;
	}
	ssdp_diffuse_sky_poa(S, *v, C->M, 1);
	// surface normal, how to get it in here?
}

void ssdp_poa_to_surface_normal(sky_pos pn0, sky_pos sn, sky_pos *pn); // orient module w.r.t ground orientation

AOI_Model_Data ssdp_init_aoi_model(AOI_Model model,double nf, double nar,double *theta, double *effT, int N);
double ssdp_diffuse_sky_poa(sky_grid * sky, sky_pos pn, AOI_Model_Data *M, int mask);					// diffuse contribution
double ssdp_direct_sky_poa(sky_grid * sky, sky_pos pn, AOI_Model_Data *M, int mask);					// direct contribution
double ssdp_total_sky_poa(sky_grid * sky, sky_pos pn, AOI_Model_Data *M, int mask);					// all sky contributions together
double ssdp_groundalbedo_poa(sky_grid * sky, double albedo, sky_pos pn, AOI_Model_Data *M, int mask);	// ground albedo contribution (with crude assumptions)
double ssdp_total_poa(sky_grid * sky, double albedo, sky_pos pn, AOI_Model_Data *M, int mask);		// sky+ground
double ssdp_diffuse_sky_horizontal(sky_grid * sky, AOI_Model_Data *M, int mask);
double ssdp_direct_sky_horizontal(sky_grid * sky, AOI_Model_Data *M, int mask);
double ssdp_total_sky_horizontal(sky_grid * sky, AOI_Model_Data *M, int mask);
*/
// PARSEFLAG list_vars VarLister ""
void VarLister(char *in)
{
	ListVars();
}
// PARSEFLAG help Help "[-l/command]"
void Help(char *in)
{
	int i=0;
	int lk;
	lk=strlen(in);
	while ((in[lk-1]==' ')&&(lk>0))
		lk--;
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



