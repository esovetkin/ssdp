#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
/* local includes */
#include "libssdp.h"
#include "lio.h"
#include "util.h"
#include "variables.h"
#include "parser.h"
#include "parserutil.h"
#include "HDF5/H5FileIO.h"
#include "HDF5/H5Datatypes.h"
#include "HDF5/H5Enums.h"

typedef enum arrayops{ARR_PLUS,ARR_MINUS,ARR_MULT,ARR_DIV} arrayops;
/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG array_eval array_comp "a=<in-array> op=<operator:+,-,*,/> b=<in-array> c=<out-array>"
DESCRIPTION Basic operations on array variables: c = a <op> b. Works both for element wise array-array operations aswell as for scalar-array operations.
ARGUMENT a input array
ARGUMENT op Operator character (+,-,*,/)
ARGUMENT b input array
OUTPUT c output array
END_DESCRIPTION
*/
void array_comp(char *in)
{
	int i;
	char *word;
	array *a, *b, c;
	int swap=0;
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
	if ((a->N!=b->N)&&(a->N!=1)&&(b->N!=1))
	{
		free(word);
		Warning("Array lengths so not match in array_eval");
		return;
	}
	if (a->N==1)
	{
		array *d;
		// swap a and b
		d=b;
		b=a;
		a=d;
		swap=1;
	}
	c.D=malloc(a->N*sizeof(double));
	c.N=a->N;
	switch (OP)
	{
		case ARR_PLUS:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]+b->D[i%b->N];
			break;
		case ARR_MINUS:
			for (i=0;i<a->N;i++)
			{
				if (swap)
					c.D[i]=b->D[i%b->N]-a->D[i];
				else
					c.D[i]=a->D[i]-b->D[i%b->N];
			}
			break;
		case ARR_MULT:
			for (i=0;i<a->N;i++)
				c.D[i]=a->D[i]*b->D[i%b->N];
			break;
		case ARR_DIV:
			for (i=0;i<a->N;i++)
			{
				if (swap)
					c.D[i]=b->D[i%b->N]/a->D[i];
				else
					c.D[i]=a->D[i]/b->D[i%b->N];
			}
			break;
		default:
			Warning("the large Hadron collider finally did destroy the world (or is it a bug?)");
			free(word);
			return;
	}	
	if (!GetArg(in, "c", word))
	{
		free(word);
		return;
	}
	printf("array_eval: creating array %s\n", word);
	if(AddArray(word, c))
	{
		free(word); // failed to make array
		free(c.D);
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
/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG read_array ReadArraysFromFile "file=<file-str> a0=<out-array> a1=<out-array> .. aN=<out-array>"
DESCRIPTION Reads columns from a file and stores them in arrays. The i-th column is stored in the i-th output array. Note that you cannot skip columns!
ARGUMENT file input filename
OUTPUT ai the i-th output array
END_DESCRIPTION
*/
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
/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG write_array WriteArraysToFile "a0=<in-array> a1=<in-array> .. aN=<in-array> file=<file-str>"
DESCRIPTION Writes arrays in columns of a file. The i-th array is written to the i-th column in the file. Note that you cannot skip columns!
ARGUMENT ai the i-th input array
OUTPUT file output filename
END_DESCRIPTION
*/
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
			printf("0: %d %d:%d\n", N, i, a->N);
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
	printf("writing arrays to file %s\n", file);
	WriteArrays(file,data,i,N);
	free(file);
	free(data);
}
/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG write_array_to_H5 WriteArraysToH5 "a0=<in-array> a1=<in-array> .. aN=<in-array> file=<file-str>"
DESCRIPTION Writes arrays in columns of a HDF5 File in a dataset called 'data'. The i-th array is written to the i-th column in the file. Note that you cannot skip columns!
ARGUMENT ai the i-th input array
OUTPUT file output filename
END_DESCRIPTION
*/
void WriteArraysToH5(char *in)
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
			printf("0: %d %d:%d\n", N, i, a->N);
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
	printf("writing arrays to file %s\n", file);
	// TODO 
	// add writing to h5 file here
	//WriteArrays(file,data,i,N);
	struct H5FileIOHandler *handler = H5FileIOHandler_init(file, W);
	if (NULL == handler){
		Warning("Error creating HDF5 file! Try using .h5 file extention.\n"); 
		free(file);
		free(data);
		return;
	}
	hid_t small_float = H5T_define_16bit_float();
	ErrorCode err;
	err = H5FileIOHandler_write_array(handler, "data", &(data[0][0]), N, i, N, small_float);
	if (SUCCESS != err){
		Warning("Error writing to HDF5 file! Try using .h5 file extention.\n");
	}
	free(handler);
	free(file);
	free(data);
}
/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG make_array MakeArray "x=<out-array> x1=<float> x2=<float> Nx=<int>"
DESCRIPTION Creates an array with Nx+1 elements ranging from x1 to and including x2.
ARGUMENT x1 Start float value
ARGUMENT x2 End float value
ARGUMENT Nx Number of steps (i.e. number of elements is Nx+1)
OUTPUT x output array
END_DESCRIPTION
*/
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


/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG make_grid MakeGrid "x=<out-array> y=<out-array> x1=<float> x2=<float> y1=<float> y2=<float> Nx=<int> Ny=<int>"
DESCRIPTION Creates two arrays with a regular grid with (Nx+1)*(Ny+1) elements.
ARGUMENT x1 Start float x value
ARGUMENT x2 End float x value
ARGUMENT y1 Start float y value
ARGUMENT y2 End float y value
ARGUMENT Nx Number of x-steps
ARGUMENT Ny Number of y-steps
OUTPUT x output array
OUTPUT y output array
END_DESCRIPTION
*/
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
	printf("Creating array %s\n",word);
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
	printf("Creating array %s\n",word);
	if(AddArray(word, b))
	{
		free(b.D);
		free(word);
	}
	return;	
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG get_grid GetGrid "C=<in-config> x=<out-array> y=<out-array>"
DESCRIPTION Extract the grid from the configured topogrid
ARGUMENT C  Simulation config with a configures topogrid
OUTPUT x output array
OUTPUT y output array
END_DESCRIPTION
*/
void GetGrid(char *in)
{
	char *word;
	array x, y;
	simulation_config *C;
	int i, j;
	double dx, dy, xx;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchConfig(in, "C", word, &C))
	{
		free(word);
		return;
	}
	if (C->grid_init==0)
	{
		Warning("Simulation config does not cintain a topogrid\n");
		free(word);
		return;
	}
	x.N=(C->Tx.Nx*C->Tx.Ny);
	x.D=malloc(x.N*sizeof(double));
	
	if (!x.D)
	{
		Warning("memory allocation failed\n");
		free(word);
		return;
	}
	y.N=x.N;
	y.D=malloc(x.N*sizeof(double));
	if (!y.D)
	{
		Warning("memory allocation failed\n");
		free(word);
		return;
	}
	dx=(C->Tx.x2-C->Tx.x1)/C->Tx.Nx;
	dy=(C->Tx.y2-C->Tx.y1)/C->Tx.Ny;
	for (i=0;i<C->Tx.Nx;i++)
	{
		xx=C->Tx.x1+i*dx;
		for (j=0;j<C->Tx.Ny;j++)
		{
			x.D[i*C->Tx.Ny+j]=xx;
			y.D[i*C->Tx.Ny+j]=C->Tx.y1+j*dy;
		}
	}	
	if (!GetArg(in, "x", word))
	{
		free(word);
		return;
	}	
	printf("Creating array %s\n",word);
	if(AddArray(word, x))
	{
		free(x.D);	
		free(word);
	}
	word=malloc((strlen(in)+1)*sizeof(char)); // allocate new, word is swallowed into the variable list by AddArray
	if (!GetArg(in, "y", word))
	{
		free(word);
		return;
	}	
	printf("Creating array %s\n",word);
	if(AddArray(word, y))
	{
		free(y.D);
		free(word);
	}
	return;	
}
/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG make_scalar MakeScalar "x=<out-array> val=<float>"
DESCRIPTION Creates an array with length 1 (simply a shorter way to create a 1 valued array than using the make_array command)
ARGUMENT val float value
OUTPUT x output array
END_DESCRIPTION
*/
void MakeScalar(char *in)
{
	char *word;
	double v;
	array a;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchFloat(in, "val", word, &v))
	{
		free(word);
		return;
	}
	if (!GetArg(in, "x", word))
	{
		free(word);
		return;
	}	
	
	a.D=malloc(sizeof(double));
	if (!a.D)
	{
		Warning("memory allocation failed\n");
		free(word);
		return;
	}
	a.N=1;
	a.D[0]=v;
	
	printf("Creating scalar %s\n",word);
	if(AddArray(word, a))
	{
		free(a.D);	
		free(word);
	}
	return;	
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG rad2deg ArrayRad2Deg "x=<in/out-array>"
DESCRIPTION Converts radians to degrees. The conversion is in-place (i.e. input and output are the same array).
ARGUMENT x input array
OUTPUT x output array
END_DESCRIPTION
*/
void ArrayRad2Deg(char *in)
{
	char *word;
	array *x;
	int i;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchArray(in, "x", word, &x))
	{
		free(word);
		return;
	}
	
	printf("converting %s to degrees\n",word);
	for(i=0;i<x->N;i++)
		x->D[i]=rad2deg(x->D[i]);
	return;	
}
/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG deg2rad ArrayDeg2Rad "x=<in/out-array>"
DESCRIPTION Converts degrees to radians. The conversion is in-place (i.e. input and output are the same array).
ARGUMENT x input array
OUTPUT x output array
END_DESCRIPTION
*/
void ArrayDeg2Rad(char *in)
{
	char *word;
	array *x;
	int i;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchArray(in, "x", word, &x))
	{
		free(word);
		return;
	}
	
	printf("converting %s to radians\n",word);
	free(word);
	for(i=0;i<x->N;i++)
		x->D[i]=deg2rad(x->D[i]);
	return;	
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG sin Sin "phi=<in-array> o=<out-array>"
DESCRIPTION Computes the sine function.
ARGUMENT phi input array with angle
OUTPUT o output array win sin(phi)
END_DESCRIPTION
*/
void Sin(char *in)
{
	char *word;
	array *p, o;
	int i;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchArray(in, "phi", word, &p))
	{
		free(word);
		return;
	}
	o.D=malloc(p->N*sizeof(double));
	o.N=p->N;
	
	for(i=0;i<p->N;i++)
		o.D[i]=sin(p->D[i]);
		
	if (!GetArg(in, "o", word))
	{
		free(word);
		return;
	}	
	printf("Creating array %s\n",word);
	if(AddArray(word, o))
	{
		free(o.D);	
		free(word);
	}
	return;	
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG cos Cos "phi=<in-array> o=<out-array>"
DESCRIPTION Computes the cosine function.
ARGUMENT phi input array with angle
OUTPUT o output array win cos(phi)
END_DESCRIPTION
*/
void Cos(char *in)
{
	char *word;
	array *p, o;
	int i;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchArray(in, "phi", word, &p))
	{
		free(word);
		return;
	}
	o.D=malloc(p->N*sizeof(double));
	o.N=p->N;
	
	for(i=0;i<p->N;i++)
		o.D[i]=cos(p->D[i]);
		
	if (!GetArg(in, "o", word))
	{
		free(word);
		return;
	}	
	printf("Creating array %s\n",word);
	if(AddArray(word, o))
	{
		free(o.D);	
		free(word);
	}
	return;	
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG perturb Perturb "x=<in/out-array> [releps=<float value>] [abseps=<float value>]"
DESCRIPTION Make infinitesimal random changes.
ARGUMENT x input array
ARGUMENT releps relative magnitude of perturbations (default 0)
ARGUMENT abseps absolute magnitude of perturbations (default 0)
OUTPUT x output array
END_DESCRIPTION
*/
void Perturb(char *in)
{
	char *word;
	double releps = 0;
	double abseps = 0;
	int i;
	array *x;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (GetOption(in, "releps", word))
	{
		releps = atof(word);
	}
	releps=fabs(releps);
	if (GetOption(in, "abseps", word))
	{
		abseps = atof(word);
	}
	abseps=fabs(abseps);	
	if (FetchArray(in, "x", word, &x))
	{
		free(word);
		return;
	}	
	printf("perturbing %s by %e (relative) and %e (absolute)\n",
		word,releps,abseps);
	
	srand(time(NULL));
	for(i=0;i<x->N;i++)
	{
		x->D[i]*=(1.0+releps*(((double)rand())/RAND_MAX-0.5));
		x->D[i]+=abseps*(((double)rand())/RAND_MAX-0.5);
	}
	return;	
}
