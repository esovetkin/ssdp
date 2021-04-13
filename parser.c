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

clock_t tic;
#define TOC() ((double)(clock()-tic)/CLOCKS_PER_SEC)
#define TIC() (tic=clock())

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
// PARSEFLAG config_coord ConfigCoord "C=<config-variable> lat=<latitude> lon=<longitude>"
void ConfigCoord (char *in)
{
	simulation_config *C;
	char *word;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "lon", word, &(C->lon)))
	{
		free(word);
		return;
	}
	C->lon=degr2rad(C->lon);
	if (FetchFloat(in, "lat", word, &(C->lat)))
	{
		free(word);
		return;
	}
	C->lat=degr2rad(C->lat);
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
	if (C->ST)
	{
		for (i=0;C->Nl;i++)
			ssdp_free_sky_transfer(C->ST+i);
	}
	free(C->ST);
	C->ST=NULL;
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
#define ProgressLen 40
#define ProgressTics 4
void InitConfigMask(simulation_config *C)
{
	int i;
	if ((C->sky_init)&&(C->topo_init)&&(C->loc_init))
	{
		double dt;
		int pco=0;
		sky_transfer ST;
		FreeConfigMask(C); // make sure we are clear to allocate new memory
		C->ST=malloc(C->Nl*sizeof(sky_transfer));
		TIC();
		for (i=0;i<C->Nl;i++)
		{
			ST=ssdp_mask_horizon(&(C->S),&(C->T),C->x[i],C->y[i],C->z[i]);
			C->ST[i]=ssdp_total_transfer(&(C->S), C->albedo, C->o[i], &(C->M), &ST);
			ssdp_free_sky_transfer(&ST);
			pco=ProgressBar((100*(i+1))/C->Nl, pco, ProgressLen, ProgressTics);
		}
		dt=TOC();
		printf("\n");
		printf("%d horizons computed in %g s (%g s/horizons)\n", C->Nl, dt, dt/((double)C->Nl));
	}
}
void InitConfigMaskNoH(simulation_config *C) // same as above but without horizon calculation
{
	int i;
	if ((C->sky_init)&&(C->loc_init))
	{
		double dt;
		int pco=0;
		FreeConfigMask(C); // make sure we are clear to allocate new memory
		C->ST=malloc(C->Nl*sizeof(sky_transfer));
		TIC();
		for (i=0;i<C->Nl;i++)
		{
			//C->ST[i]=ssdp_total_transfer(&(C->S), C->albedo, C->o[i], &(C->M), NULL);
			C->ST[i]=ssdp_sky_transfer(&(C->S), C->o[i], &(C->M), NULL);
			pco=ProgressBar((100*(i+1))/C->Nl, pco, ProgressLen, ProgressTics);
		}
		dt=TOC();
		printf("\n");
		printf("%d transfer functions computed in %g s (%g s/horizons)\n", C->Nl, dt, dt/((double)C->Nl));
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
// add albedo to this routine I suppose (i.e. enable locally varying albedo values)
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
	
	if ((z->N==1)&&(az->N==1)) // one z, one orientation
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
typedef enum arrayops{ARR_PLUS,ARR_MINUS,ARR_MULT,ARR_DIV} arrayops;
// PARSEFLAG array_array_comp array_array_comp "a=<a-array-variable> op=<operator:+,-,*,/> b=<b-array-variable> c=<c-output-array>"
void array_array_comp(char *in)
{
	int i;
	char *word;
	array *a, *b, c;
	arrayops OP;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchArray(in, "a", word, &a))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "b", word, &b))
	{
		free(word);
		return;
	}
	if (!GetArg(in, "op", word))
	{
		free(word);
		return;
	}
	else
	{
		if (strlen(word)!=1)
		{
			Warning("Unknown Operator %s\n", word);
			free(word);
			return;
		}
		switch(*word)
		{
			case '+':
				OP=ARR_PLUS;
				break;
			case '-':
				OP=ARR_MINUS;
				break;
			case '*':
				OP=ARR_MULT;
				break;
			case '/':
				OP=ARR_DIV;
				break;
			default:
			{
				Warning("Unknown Operator %s\n", word);
				free(word);
				return;
			}
		}
	}	
	if (a->N!=b->N)
	{
		free(word);
		Warning("Cannot Add Arrays, arrays not of same length");
		return;
	}
	c.D=malloc(a->N*sizeof(double));
	c.N=a->N;
	switch (OP)
	{
		case ARR_PLUS:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]+b->D[i];
			break;
		case ARR_MINUS:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]-b->D[i];
			break;
		case ARR_MULT:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]*b->D[i];
			break;
		case ARR_DIV:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]/b->D[i];
			break;
		default:
			Warning("the large Hadron collider finally did destroy the world");
			free(word);
			return;
	}	
	if (!GetArg(in, "c", word))
	{
		free(word);
		return;
	}
	printf("creating array %s\n", word);
	if(AddArray(word, c))
	{
		free(word); // failed to make array
		free(c.D);
	}	
}

// PARSEFLAG array_scalar_comp array_scalar_comp "a=<a-array-variable> op=<operator:+,-,*,/> b=<b-array-variable> c=<c-output-array>"
void array_scalar_comp(char *in)
{
	int i;
	char *word;
	array *a, c;
	double b;
	arrayops OP;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchArray(in, "a", word, &a))
	{
		free(word);
		return;
	}	
	if (FetchFloat(in, "b", word, &b))
	{
		Warning("Could not get froat value from %s\n", word);
		free(word);
		return;
	}
	if (!GetArg(in, "op", word))
	{
		free(word);
		return;
	}
	else
	{
		if (strlen(word)!=1)
		{
			Warning("Unknown Operator %s\n", word);
			free(word);
			return;
		}
		switch(*word)
		{
			case '+':
				OP=ARR_PLUS;
				break;
			case '-':
				OP=ARR_MINUS;
				break;
			case '*':
				OP=ARR_MULT;
				break;
			case '/':
				OP=ARR_DIV;
				break;
			default:
			{
				Warning("Unknown Operator %s\n", word);
				free(word);
				return;
			}
		}
	}	
	c.D=malloc(a->N*sizeof(double));
	c.N=a->N;
	switch (OP)
	{
		case ARR_PLUS:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]+b;
			break;
		case ARR_MINUS:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]-b;
			break;
		case ARR_MULT:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]*b;
			break;
		case ARR_DIV:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]/b;
			break;
		default:
			Warning("the large Hadron collider finally did destroy the world");
			free(word);
			return;
	}
	
	if (!GetArg(in, "c", word))
	{
		free(word);
		return;
	}
	printf("creating array %s\n", word);
	if(AddArray(word, c))
	{
		free(word); // failed to make array
		free(c.D);
	}	
}
// PARSEFLAG sample_topo SampleTopography "C=<config-variable>  x=<x-array-variable> y=<y-array-variable> z=<z-output-array> azimuth=<azimuth-output-array> zenith=<zenith-output-array>"
void SampleTopography(char *in)
{
	int i;
	char *word;
	simulation_config *C;
	array *x, *y, z, azi, zen;
	sky_pos sn;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}	
	if (!C->topo_init) // TODO: in this case just omit the horizon and compute only one location?
	{	
		Warning("Simulation config has no topology initialized\n");
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
	if (x->N!=y->N)
	{
		Warning("Length of x- and y-arrays do not match\n");
		free(word);
		return;
	}
	z.D=malloc(x->N*sizeof(double));
	z.N=x->N;
	azi.D=malloc(x->N*sizeof(double));
	azi.N=x->N;
	zen.D=malloc(x->N*sizeof(double));
	zen.N=x->N;
	for (i=0;i<x->N;i++)
	{
		z.D[i]=ssdp_sample_topology(x->D[i], y->D[i], &(C->T),&sn);
		azi.D[i]=sn.a;
		zen.D[i]=sn.z;
	}
	
	if (GetOption(in, "z", word))
	{
		printf("Creating array %s\n",word);
		if(AddArray(word, z))
		{
			Warning("Failed to create array %s\n",word);
			free(z.D);	
		}
		else
			word=malloc((strlen(in)+1)*sizeof(char));
	}
	else
		free(z.D);
		
	if (GetOption(in, "azimuth", word))
	{
		printf("Creating array %s\n",word);
		if(AddArray(word, azi))
		{
			Warning("Failed to create array %s\n",word);
			free(azi.D);	
		}
		else
			word=malloc((strlen(in)+1)*sizeof(char));
	}
	else
		free(azi.D);
		
	if (GetOption(in, "zenith", word))
	{
		printf("Creating array %s\n",word);
		if(AddArray(word, zen))
		{
			Warning("Failed to create array %s\n",word);
			free(zen.D);	
		}
		else
			word=NULL;
	}
	else
		free(zen.D);
	if (word)
		free(word);
}
// PARSEFLAG offset_topo OffsetTopography "C=<config-variable>  o=<offset-value> x=<x-array-variable> y=<y-array-variable> xoff=<offset-x-output-array> yoff=<offset-x-output-array> zoff=<offset-z-output-array>"
void OffsetTopography(char *in)
{
	int i;
	char *word;
	simulation_config *C;
	array *x, *y, zoff, yoff, xoff;
	double o;
	sky_pos sn;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}	
	if (!C->topo_init)
	{	
		Warning("Simulation config has no topology initialized\n");
		free(word);
		return;
	}	
	if (FetchFloat(in, "o", word, &o))
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
	if (x->N!=y->N)
	{
		Warning("Length of x- and y-arrays do not match\n");
		free(word);
		return;
	}
	zoff.D=malloc(x->N*sizeof(double));
	zoff.N=x->N;
	yoff.D=malloc(x->N*sizeof(double));
	yoff.N=x->N;
	xoff.D=malloc(x->N*sizeof(double));
	xoff.N=x->N;
	for (i=0;i<x->N;i++)
	{
		zoff.D[i]=ssdp_sample_topology(x->D[i], y->D[i], &(C->T),&sn)+o*cos(sn.z);
		xoff.D[i]=x->D[i]+o*sin(sn.z)*cos(sn.a);
		yoff.D[i]=y->D[i]+o*sin(sn.z)*sin(sn.a);
	}
	
	
	if (!GetArg(in, "zoff", word))
	{
		free(word);
		return;
	}
	printf("Creating array %s\n",word);
	if(AddArray(word, zoff))
	{
		Warning("Failed to create array %s\n",word);
		free(zoff.D);	
	}
	else
		word=malloc((strlen(in)+1)*sizeof(char));
		
	if (!GetArg(in, "xoff", word))
	{
		free(word);
		return;
	}	
	printf("Creating array %s\n",word);
	if(AddArray(word, xoff))
	{
		Warning("Failed to create array %s\n",word);
		free(xoff.D);	
	}
	else
		word=malloc((strlen(in)+1)*sizeof(char));
		
	if (!GetArg(in, "yoff", word))
	{
		free(word);
		return;
	}
	printf("Creating array %s\n",word);
	if(AddArray(word, yoff))
	{
		Warning("Failed to create array %s\n",word);
		free(yoff.D);	
	}
	else
		word=NULL;
	
	if (word)
		free(word);
}
// PARSEFLAG rotate_POAto_surface RotatePOA "poa_a=<azimuth-value>  poa_z=<zenith-value> surf_a=<azimuth-array-variable> surf_z=<zenith-array-variable> out_a=<azimuth-output-array> out_z=<zenith-output-array>"
void RotatePOA(char *in)
{
	int i;
	char *word;
	array *rot_a, *rot_z, azi, zen;
	sky_pos poa, poa0;
	sky_pos sn;
	word=malloc((strlen(in)+1)*sizeof(char));
	
	if (FetchFloat(in, "poa_a", word, &poa0.a))
	{
		free(word);
		return;
	}
	if (FetchFloat(in, "poa_z", word, &poa0.z))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "surf_a", word, &rot_a))
	{
		free(word);
		return;
	}
	if (FetchArray(in, "surf_z", word, &rot_z))
	{
		free(word);
		return;
	}
	if (rot_a->N!=rot_z->N)
	{
		Warning("Length of azimuth and zenoth of the surface normal do not match\n");
		free(word);
		return;
	}
	azi.D=malloc(rot_a->N*sizeof(double));
	azi.N=rot_a->N;
	zen.D=malloc(rot_a->N*sizeof(double));
	zen.N=rot_a->N;
	for (i=0;i<rot_a->N;i++)
	{
		sn.a=rot_a->D[i];
		sn.z=rot_z->D[i];
		ssdp_poa_to_surface_normal(poa0, sn, &poa); // orient module w.r.t ground orientation
		azi.D[i]=poa.a;
		zen.D[i]=poa.z;
	}
	
	
	if (!GetArg(in, "out_a", word))
	{
		free(word);
		return;
	}
	printf("Creating array %s\n",word);
	if(AddArray(word, azi))
	{
		Warning("Failed to create array %s\n",word);
		free(azi.D);	
	}
	else
		word=malloc((strlen(in)+1)*sizeof(char));
		
	if (!GetArg(in, "out_z", word))
	{
		free(word);
		return;
	}	
	printf("Creating array %s\n",word);
	if(AddArray(word, zen))
	{
		Warning("Failed to create array %s\n",word);
		free(zen.D);	
	}
	else
		word=NULL;
	if (word)
		free(word);
}


// PARSEFLAG sim_static SimStatic "C=<config-variable> t=<array-variable> GHI=<array-variable> DHI=<array-variable> POA=<out-array>"
void SimStatic(char *in)
{
	int i, j; // loop through space and time
	char *word;
	simulation_config *C;
	array *t, *GH, *DH, out;
	clock_t tsky0, tpoa0;
	clock_t tsky=0, tpoa=0;
	double ttsky, ttpoa;
	int pco=0;
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
	if (!C->topo_init) 
	{	
		Warning("No topological data available, omitting horizon\n");
		InitConfigMaskNoH(C);
	}
	if (!C->loc_init) 
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
		tsky0=clock();
		ssdp_make_perez_all_weather_sky_coordinate(&(C->S), (time_t) t->D[j], C->lon, C->lat, GH->D[j], DH->D[j]);
		tsky+=clock()-tsky0;
		tpoa0=clock();
		for (i=0;i<C->Nl;i++)
		{
			// compute POA at evert location  instance
			// this operation can be somewhat expensive
			// we should consider to compute a transfer efficiency for every sky element and every location
			// this could speed things up quite a bit. The transfer efficiency should entail
			// cos(z) w.r.t the POA, AOI effects and ground albedo
			out.D[j*C->Nl+i]=ssdp_total_poa(&(C->S),C->ST+i);
		}
		pco=ProgressBar((100*((j+1)*C->Nl))/(t->N*C->Nl), pco, ProgressLen, ProgressTics);
		tpoa+=clock()-tpoa0;
	}
	ttsky=(double)tsky/CLOCKS_PER_SEC;
	ttpoa=(double)tpoa/CLOCKS_PER_SEC;
	printf("Computed %d skies in %g s (%g s/sky)\n", t->N, ttsky, ttsky/((double)t->N));
	printf("Computed %d POA Irradiances in %g s (%g s/POA)\n", t->N*C->Nl, ttpoa, ttpoa/((double)(t->N*C->Nl)));
	if(AddArray(word, out))
	{
		free(word); // failed to make array
		free(out.D);
	}	
}
// _PARSEFLAG sim_route SimRoute "C=<config-variable> x=<array-variable> y=<array-variable> t=<array-variable> GHI=<array-variable> DHI=<array-variable> POA=<out-array>"
/*void SimStatic(char *in)
{
	sky_ST *ST; // a mask for every grid point
	int i, j; // loop through space and time
	
	
}*/

// PARSEFLAG solpos SolarPos "t=<array-variable> lon=<longitude> lat=<latitude> azimuth=<output-array> zenith=<output-array>"
void SolarPos(char *in)
{
	int i;
	char *word;
	sky_pos s;
	array *t, azi, zen;
	double lon, lat;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchArray(in, "t", word, &t))
	{
		free(word);
		return;
	}
	
	if (FetchFloat(in, "lon", word, &lon))
	{
		free(word);
		return;
	}
	lon=degr2rad(lon);
	if (FetchFloat(in, "lat", word, &lat))
	{
		free(word);
		return;
	}
	lat=degr2rad(lat);
	azi.D=malloc(t->N*sizeof(double));
	if (azi.D==NULL)
	{
		free(word);
		return;
	}	
	azi.N=t->N;
	zen.D=malloc(t->N*sizeof(double));
	if (zen.D==NULL)
	{
		free(word);
		return;
	}	
	zen.N=t->N;
	for (i=0;i<t->N;i++)
	{
		s=ssdp_sunpos((time_t)t->D[i], lat, lon);
		azi.D[i]=s.a;
		zen.D[i]=s.z;
	}
	if (!GetArg(in, "azimuth", word))
	{
		free(word);
		return;
	}
	if(AddArray(word, azi))
	{
		free(word); // failed to make array
		free(azi.D);
	}	
	word=malloc((strlen(in)+1)*sizeof(char));
	if (!GetArg(in, "zenith", word))
	{
		free(word);
		return;
	}
	if(AddArray(word, zen))
	{
		free(word); // failed to make array
		free(zen.D);
	}	
}
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
	int i, k=0, Na=4, N;
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
	if (N>0)
	{
		array a;
		a.N=N;
		// create i arrays
		for (k=0;k<i;k++)
		{
			a.D=data[k];
			printf("Creating array %s\n", names[k]);
			if(AddArray(names[k], a))
			{
				free(names[k]); // failed to make array
				free(a.D);
			}
		}
	}
	else
	{
		printf("Could not parse file %s\n", file);
		for (k=0;k<i;k++)
		{
			free(data[k]);
			free(names[k]);
		}
	}
	free(file);
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



