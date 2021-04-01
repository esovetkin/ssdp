#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libssdp.h"
#include "variables.h"
#include "util.h"
typedef enum {TOPOVAR, SKYDOMEVAR, VECVAR, SIMCONF} VarType;
typedef struct {
	topology T;
	sky_grid S;
	sky_pos v;
	simulation_config C;
	char *name;
	VarType type;
} Var;

#define BLOCK 10
Var *variables=NULL;
int Nvar=0;

simulation_config InitConf()
{
	simulation_config C;
	C.M.M=AOI_NONE;
	C.M.ng=1.5;
	C.M.nar=1.5;
	C.M.theta=NULL;
	C.M.effT=NULL;
	C.M.N=0;
	C.P=POA_ZREL;
	return C;
}
void FreeConf(simulation_config *C)
{
	C->M.N=0;
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
			case TOPOVAR:
				ssdp_free_topology(&(variables[i].T));
				break;
			case SKYDOMEVAR:
				ssdp_free_sky(&(variables[i].S));
				break;
			case SIMCONF:
				FreeConf(&(variables[i].C));
				break;
			case VECVAR:
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
			case TOPOVAR:
				printf("Topology          :  ");
				break;
			case SKYDOMEVAR:
				printf("Kky DOme          :  ");
				break;
			case VECVAR:
				printf("Vector            :  ");
				break;
			case SIMCONF:
				printf("Simulation Config :  ");
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
	for (i=0;i<Nvar;i++)
	{
		if (strlen(variables[i].name)==lk)
			if (strncmp(variables[i].name, name, lk)==0)
				return i;
	}
	return -1;	
}

int LookupTopology(char *name, topology **T)
{
	int i;
	if ((i=lookupvar(name))<0)
		return 0;
	else
	{
		if (variables[i].type==TOPOVAR)
		{
			(*T)=&variables[i].T;
			return 1;
		}
	}
	return 0;
} 

int LookupSkyDome(char *name, sky_grid **S)
{
	int i;
	if ((i=lookupvar(name))<0)
		return 0;
	else
	{
		if (variables[i].type==SKYDOMEVAR)
		{
			(*S)=&variables[i].S;
			return 1;
		}
	}
	return 0;
} 

int LookupVec(char *name, sky_pos **v)
{
	int i;
	if ((i=lookupvar(name))<0)
		return 0;
	else
	{
		if (variables[i].type==VECVAR)
		{
			(*v)=&variables[i].v; // return a pointer directly into the veriable
			return 1;
		}
	}
	return 0;
} 
int LookupSimConf(char *name, simulation_config **C)
{
	int i;
	if ((i=lookupvar(name))<0)
		return 0;
	else
	{
		if (variables[i].type==SIMCONF)
		{
			(*C)=&variables[i].C;
			return 1;
		}
	}
	return 0;
} 

int AddTopology(char *name, topology T)
{
	int i;
	if ((i=lookupvar(name))<0)
	{
		variables[Nvar].name=name;
		variables[Nvar].T=T;
		variables[Nvar].type=TOPOVAR;
		Nvar++;
		if (Nvar%BLOCK==0)
			variables=realloc(variables, Nvar+BLOCK);
		return 0;
	}
	if (variables[i].type==TOPOVAR)	
	{
		ssdp_free_topology(&(variables[i].T));
		variables[i].T=T;
		return 0;
	}
	else
		Warning("Variable %s exists as a different type\n", name);
	return 1;
}

int AddSkyDome(char *name, sky_grid S)
{
	int i;
	if ((i=lookupvar(name))<0)
	{
		variables[Nvar].name=name;
		variables[Nvar].S=S;
		variables[Nvar].type=SKYDOMEVAR;
		Nvar++;
		if (Nvar%BLOCK==0)
			variables=realloc(variables, Nvar+BLOCK);
		return 0;
	}
	if (variables[i].type==SKYDOMEVAR)	
	{
		ssdp_free_sky(&(variables[i].S));
		variables[i].S=S;
		return 0;
	}
	else
		Warning("Variable %s exists as a different type\n", name);
	return 1;
}

int AddVec(char *name, sky_pos v)
{
	int i;
	if ((i=lookupvar(name))<0)
	{
		variables[Nvar].name=name;
		variables[Nvar].v=v;
		variables[Nvar].type=VECVAR;
		Nvar++;
		if (Nvar%BLOCK==0)
			variables=realloc(variables, Nvar+BLOCK);
		return 0;
	}
	if (variables[i].type==VECVAR)	
	{
		variables[i].v=v;
		return 0;
	}
	else
		Warning("Variable %s exists as a different type\n", name);
	return 1;
}

int AddSimConf(char *name, simulation_config C)
{
	int i;
	if ((i=lookupvar(name))<0)
	{
		variables[Nvar].name=name;
		variables[Nvar].C=C;
		variables[Nvar].type=SIMCONF;
		Nvar++;
		if (Nvar%BLOCK==0)
			variables=realloc(variables, Nvar+BLOCK);
		return 0;
	}
	if (variables[i].type==SIMCONF)	
	{
		FreeConf(&(variables[i].C));
		variables[i].C=C;
		return 0;
	}
	else
		Warning("Variable %s exists as a different type\n", name);
	return 1;
}

int RMVar(char *name)
{
	int i;
	if ((i=lookupvar(name))>-1)
	{
		switch(variables[i].type)
		{
			case TOPOVAR:
				ssdp_free_topology(&(variables[i].T));
				break;
			case SKYDOMEVAR:
				ssdp_free_sky(&(variables[i].S));
				break;
			case SIMCONF:
				FreeConf(&(variables[i].C));
				break;
			case VECVAR:
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


