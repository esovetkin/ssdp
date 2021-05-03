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
#include "parser.h"
#include "parserutils.h"

clock_t tic;
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

/* 
 * To automate the documentation process a bit every parser gets a 
 * description in the source which is automatically added to the ssdp 
 * man page. The format is as follows:
 * In comments use the BEGIN_DESCRIPTION/END_DESCRIPTION flahs to mark the beginning
 * and end of a man page entry.
 * The parseflag was already used to generate the parsedef.h data and is also used here
 * I like to put them together.
 * the DESCRIPTION flag is followed by a description of the command (one line!)
 * The ARGUMENT flag described one argument
 * The OUTPUT flag describes an output of the function
 */
/*
BEGIN_DESCRIPTION
SECTION General
PARSEFLAG list_vars VarLister ""
DESCRIPTION List all variables
END_DESCRIPTION
*/
void VarLister(char *in)
{
	ListVars();
}

/*
BEGIN_DESCRIPTION
SECTION General
PARSEFLAG help Help "[-l/command]"
DESCRIPTION Print basic help. If no arguments or the -l argument is given this will list all availablke help data. You can also request help on a specific command.
END_DESCRIPTION
*/
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
			if (KeyTable[i].key)
			{
				printf("Command %s:\n",KeyTable[i].key);
				SplitWords(Usage[i], "\t");
				printf("\n");
			}
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



