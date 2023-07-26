#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libssdp.h"
#include "variables.h"
#include "util.h"
typedef enum {SIMCONF, ARRAY} VarType;
typedef struct {
	array a;
	simulation_config C;
	char *name;
	VarType type;
} Var;

#define BLOCK 10
Var *variables=NULL;
int Nvar=0;


// the default location ios the Forschungszentrum Jülich
// gotta be somewhere and I happen to be here alot - bp
#define FZLAT 50.902996388
#define FZLON 6.407165038
simulation_config InitConf()
{
	simulation_config C;
	C.M.M=AOI_NONE;
	C.M.ng=1.5;
	C.M.nar=1.5;
	C.M.theta=NULL;
	C.M.effT=NULL;
	C.M.N=0;
	C.sky_init=0;
	C.topo_init=0;
	C.grid_init=0;
	C.albedo=0.0; 
	C.lon=FZLON;
	C.lat=FZLAT;
	C.E=0;
	C.loc_init=0;
	C.x=NULL;
	C.y=NULL;
	C.z=NULL;
	C.o=NULL;
	C.L=NULL;
	C.Nl=0;
	return C;
}
void FreeConf(simulation_config *C)
{
	int i;
	C->M.N=0;
	C->M.M=AOI_NONE;
	if (C->M.theta)
	{
		free(C->M.theta);
		C->M.theta=NULL;
	}
	if (C->M.effT)
	{
		free(C->M.effT);
		C->M.effT=NULL;
	}
	if (C->sky_init)
	{
		ssdp_free_sky(&C->S);
		C->sky_init=0;
	}
	if (C->topo_init)
	{
		ssdp_free_topology(&C->T);
		C->topo_init=0;
	}
	if (C->grid_init)
	{
		ssdp_free_topogrid(&C->Tx);
		C->grid_init=0;
	}
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
		if (C->L)
		{
			for (i=0;i<C->Nl;i++)
				ssdp_free_location(C->L+i);
			free(C->L);
			C->L=NULL;
		}
		C->Nl=0;
		C->loc_init=0;
	}
	
}

void FreeArray(array *a)
{
	free(a->D);
	a->D=NULL;
	a->N=0;
}
	
void InitVars()
{
	if (variables)
		return; /* will only initialize once */
	variables=calloc(BLOCK, sizeof(Var));
	Nvar=0;
}
void ClearVars()
{
	int i;
	if (!variables)
		return;
	for (i=0;i<Nvar;i++)
	{
		switch(variables[i].type)
		{
			case SIMCONF:
				FreeConf(&(variables[i].C));
				break;
			case ARRAY:
				FreeArray(&(variables[i].a));
				break;
			default:
				break;
		}
		free(variables[i].name);
	}
	free(variables);
	variables=NULL;
	Nvar=0;
}
void ListVars()
{
	int i=0;
	for (i=0;i<Nvar;i++)
	{
		switch(variables[i].type)
		{
			case SIMCONF:
				printf("Simulation Config :  ");
				break;
			case ARRAY:
				printf("Array             :  ");
				break;
			default:
				break;
		}
		printf("\"%s\"", variables[i].name);
		printf("\n");
	}
}


int lookupvar(char *name)
{
	int i=0;
	int lk;
	lk=strlen(name);

    if (name[0] <= '@')
        return -1;

	for (i=0;i<Nvar;i++)
	{
		if (strlen(variables[i].name)==lk)
			if (strncmp(variables[i].name, name, lk)==0)
				return i;
	}
	return -1;	
}

/* to define search functions we use a macro
 * This creates functions of the form:
 * int <fn_name>(char *name, <t_name> **<e_name>);
 * Theres functions serve to serve a pointer to data in the variable array
 * The macro also needs the kay-name <k_name> which is the key from the 
 * VarType typedef
 */
#define LookUpVar(fn_name,t_name,e_name, k_name) \
int fn_name(char *name, t_name **e_name) \
{ \
	int i; \
	if ((i=lookupvar(name))<0) \
		return 0; \
	else \
	{ \
		if (variables[i].type==k_name) \
		{ \
			(*e_name)=&variables[i].e_name; \
			return 1; \
		} \
	} \
	return 0;\
}
/* likewise a macro to define add functions
 * This creates functions of the form:
 * int <fn_name>(char *name, <t_name> <e_name>);
 * There functions serve to put data in the variable array
 * The macro also needs the kay-name <k_name> which is the key from the 
 * VarType typedef, and a free command <freecomm> of the form
 * <free-function>(&(variables[i].<e_name>));
 * If no fee command is needed leave this argument blank
 */
#define AddVar(fn_name,t_name,e_name, k_name, freecomm) \
int fn_name(char *name, t_name e_name) \
{ \
	int i;\
	if ((i=lookupvar(name))<0)\
	{\
		variables[Nvar].name=name;\
		variables[Nvar].e_name=e_name;\
		variables[Nvar].type=k_name;\
		Nvar++;\
		if (Nvar%BLOCK==0)\
			variables=realloc(variables, (Nvar+BLOCK)*sizeof(Var));\
		return 0;\
	}\
	if (variables[i].type==k_name)\
	{\
		freecomm\
		free(variables[i].name);\
		variables[i].name=name;\
		variables[i].e_name=e_name;\
		return 0;\
	}\
	else\
		Warning("Variable %s exists as a different type\n", name);\
	return 1;\
}
// use the macros to create the functions we need
LookUpVar(LookupSimConf,simulation_config,C, SIMCONF)
LookUpVar(LookupArray,array,a, ARRAY)

AddVar(AddSimConf,simulation_config,C, SIMCONF, FreeConf(&(variables[i].C));)
AddVar(AddArray,array,a, ARRAY,FreeArray(&(variables[i].a));)

int RMVar(char *name)
{
	int i;
	if ((i=lookupvar(name))>-1)
	{
		switch(variables[i].type)
		{
			case SIMCONF:
				FreeConf(&(variables[i].C));
				break;
			case ARRAY:
				FreeArray(&(variables[i].a));
				break;
			default:
				break;
		}
		free(variables[i].name);
		i++;
		for (;i<Nvar;i++)
			variables[i-1]=variables[i];
		Nvar--;
		if (Nvar%BLOCK==BLOCK-1)
			variables=realloc(variables, Nvar+1);
		return 0;
	}
	Warning("Cannot remove variable %s, variable does not exist\n", name);
	return 1;
}


