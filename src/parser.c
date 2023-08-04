#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
/* local includes */
#include "parsedef.h"
#include "libssdp.h"
#include "lio.h"
#include "util.h"
#include "variables.h"
#include "parser.h"
#include "parserutil.h"

#ifdef OPENMP
	double tic=0;
#else // OPENMP
	clock_t tic;
#endif // OPENMP

int ParseLineNr=0;
char *ParseFileStr=NULL;

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
			Warning("Input config %s=%s, %s is not available\n",pat,str, str);
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
			Warning("Input array %s=%s, %s is not available\n",pat, str, str);
			return 1;
		}
		return 0;
	}
	Warning("Input array %s is not specified\n",pat);
	return 1;
}
int FetchOptFloat(const char *in, const char *pat, char *str, double *a)
{
        array *x;
        if (! GetOption(in, pat, str))
                goto noarg;

        if (LookupArray(str, &x))
                goto array;

        (*a)=atof(str);
        return 0;
array:
        if (1 != x->N)
                Warning("non-unit array provided for %s\n", pat);
        (*a) = x->D[0];
        return 0;
noarg:
        return 1;
}
int FetchFloat(const char *in, const char *pat, char *str, double *a)
{
        int x = FetchOptFloat(in, pat, str, a);
        if (x)
                Warning("Error: input float %s is not specified\n", pat);
        return x;
}
int FetchOptInt(const char *in, const char *pat, char *str, int *a)
{
        array *x;
        if (! GetOption(in, pat, str))
                goto noarg;

        if (LookupArray(str, &x))
                goto array;

        (*a)=atoi(str);
        return 0;
array:
        if (1 != x->N)
                Warning("Warning: non-unit array provided for %s\n", pat);
        (*a) = (int) x->D[0];
        return 0;
noarg:
        return 1;
}
int FetchInt(const char *in, const char *pat, char *str, int *a)
{
        int x = FetchOptInt(in, pat, str, a);
        if (x)
                Warning("Error: input integer %s is not specified\n", pat);
        return x;
}

/* 
 * To automate the documentation process a bit every parser gets a 
 * description in the source which is automatically added to the ssdp 
 * man page. The format is as follows:
 * In comments use the BEGIN_DESCRIPTION/END_DESCRIPTION flags to mark the beginning
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
DESCRIPTION Print basic help. If no arguments or the -l argument is given this will list all available help data. You can also request help on a specific command.
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


int GetNumOption(const char *in, const char *opt, int i, char *word)
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


struct cvec* cvec_init(int a) {
        struct cvec *self;
        self = malloc(sizeof(*self));
        if (NULL == self) goto cvec_init_eself;

        self->astep = a > 0 ? a : 1;
        self->n = 0;
        self->s = malloc(self->astep * sizeof(*self->s));
        if (NULL == self->s) goto cvec_init_es;

        return self;
cvec_init_es:
        free(self);
cvec_init_eself:
        return NULL;
}


void cvec_free(struct cvec *self)
{
        if (NULL == self)
                return;

        if (self->s) {
                int i;
                for (i=0; i<self->n; ++i)
                        free(self->s[i]);
                free(self->s);
        }
        self->s = NULL;
        free(self);
        self = NULL;
}


int cvec_push(struct cvec *self, char *word)
{
        if (self->n == self->a) {
                self->s = realloc(self->s, (self->a + self->astep)*sizeof(*self->s));
                if (NULL == self->s) goto cvec_push_erealloc;
                self->a += self->astep;
        }

        self->s[self->n++] = word;
        return 0;
cvec_push_erealloc:
        return -1;
}

int read_filelist(char *fn, struct cvec *dst)
{
        FILE *fd;
        char *x = NULL, *line = NULL;
        size_t len = 0, read;

        if (NULL == (fd=fopen(fn, "rb"))) goto efopen;
        while (-1 != (read = getline(&line, &len, fd))) {
                x = malloc((read+1)*sizeof(*x));
                if (NULL == x) goto epush;
                strncpy(x, line, read);
                x[strcspn(x, "\r\n")] = 0;
                if (cvec_push(dst, x)) goto epush;
        }

        fclose(fd);
        return 0;
epush:
        free(x);
        fclose(fd);
efopen:
        return -1;
}
