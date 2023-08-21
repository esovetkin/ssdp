#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
/* local includes */
#include "libssdp.h"
#include "lio.h"
#include "pngout.h"
#include "util.h"
#include "variables.h"
#include "parser.h"
#include "parserutil.h"
#include "h5io.h"

#define ARRAY_BS 4

#ifdef RUNMEMTEST
#include "random_fail_malloc.h"
#define malloc(x) random_fail_malloc(x)
#endif

typedef enum arrayops{ARR_PLUS,ARR_MINUS,ARR_MULT,ARR_DIV} arrayops;
/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG array_eval array_comp "a=<in-array> op=<operator:+,-,*,/> b=<in-array> c=<out-array>"
DESCRIPTION Basic operations on array variables: c = a <op> b. Works both for element wise array-array operations as well as for scalar-array operations.
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

static void* reallocate(void *data, size_t size, int n, int *bs)
{
        if (n < *bs - 1)
                return data;

        *bs += ARRAY_BS;
        void *tmp;
        if (NULL == (tmp=realloc(data, *bs*size))) goto erealloc;

        return tmp;
erealloc:
        return NULL;
}


static int getnames(char *in, char ***names, int *len)
{
        char *word;
        int n = 0, bs = ARRAY_BS;

        if (NULL == (word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
        if (NULL == (*names=malloc(bs*sizeof(**names)))) goto ename;
        *len = -1;

        while (GetNumOption(in, "a", n, word)) {
                (*names)[n] = word;
                n++;
                if (NULL == (word=malloc((strlen(in)+1)*sizeof(*word)))) goto eloop;
                *names = (char **) reallocate(*names, sizeof(**names), n, &bs);
                if (NULL == *names) goto eloop;
        }

        if (0 == n) {
                Warning("Error: array names arguments is missing!\n");
                goto enoargs;
        }

        *len=n;
        free(word);
        return 0;
enoargs:
eloop:
        for (int i=0; i < n; ++i)
                free((*names)[i]);
        free(*names);
ename:
        free(word);
eword:
        return -1;
}


static int getarrays(char *in, double ***data, int *len, int *narr)
{
        char *word;
        int ncols = 0, bs=ARRAY_BS;
        array *a;
        *len = -1;

        if (NULL == (word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
        if (NULL == (*data=malloc(bs*sizeof(**data)))) goto edata;

        while (GetNumOption(in, "a", ncols, word)) {
                if (!LookupArray(word, &a)) goto earray;
                (*data)[ncols] = a->D;
                if (*len < 0) *len = a->N;

                if (*len != a->N) {
                        Warning("Error: len(%s) differs from others\n", word);
                        goto elen;
                }
                ncols++;
                *data = (double**) reallocate(*data, sizeof(**data), ncols, &bs);
                if (NULL == *data) goto eralloc;
        }

        if (0 == ncols) {
                Warning("Error: array names arguments is missing!\n");
                goto enoargs;
        }

        *narr=ncols;
        free(word);
        return 0;
enoargs:
eralloc:
elen:
earray:
        free(*data);
edata:
        free(word);
eword:
        return -1;
}


static int addarray(double *data, int n, const char* name)
{
        array *a = malloc(sizeof(*a));
        if (NULL == a) goto ea;

        a->D = data;
        a->N = n;

        if (AddArray((char*)name, *a)) goto eadd;
        printf("Created array %s\n", name);

        free(a);
        return 0;
eadd:
        free(a);
ea:
        return -1;
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
        int i, arrlen, narr;
        double **data;
        char *file, **names;

        if (NULL == (file=malloc((strlen(in)+1)*sizeof(*file)))) goto eword;
        if (!GetArg(in, "file", file)) {
                Warning("Error: file argument is missing!\n");
                goto efile;
        }
        printf("Reading from %s\n", file);

        if (getnames(in, &names, &narr)) goto enames;
        TIC();
        if (NULL==(data=ReadArrays(file, narr, &arrlen))) goto eread;
        printf("Read from %s in %g s\n", file, TOC());

        for (i=0; i < narr; ++i)
                if (addarray(data[i], arrlen, names[i])) goto eadd;

        free(names);
        free(data);
        free(file);
        return;
eadd:
        for (int j=i; j < narr; ++j) {
                free(data[j]);
                free(names[j]);
        }
        free(data);
eread:
        free(names);
enames:
efile:
        free(file);
eword:
        Warning("Error: read_array failed!\n");
        return;
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
        int narr, arrlen;
        char *file;
        double **data;

        if (NULL==(file=malloc((strlen(in)+1)*sizeof(*file)))) goto eword;
        if (!GetArg(in, "file", file))
        {
                Warning("Error: file argument is missing!\n");
                goto efile;
        }
        printf("Writing arrays to file %s\n", file);

        if (getarrays(in, &data, &arrlen, &narr)) goto edata;

        TIC();
        WriteArrays(file, data, narr, arrlen);
        printf("Wrote to %s in %g s\n", file, TOC());

        free(data);
        free(file);
        return;
        free(data);
edata:
efile:
        free(file);
eword:
        Warning("Error: write_array failed!\n");
        return;
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG write_h5 WriteH5 "a0=<in-array> a1=<in-array> .. aN=<in-array> file=<file-str> [dataset=<str>] [type=<str>] [gzip=<int-value>] [chunksize=<int-value>] [cachemb=<int-value>] [cacheslots=<int-value>]"
DESCRIPTION Write arrays in columns of an HDF5-file in a dataset. The i-th array is written to the i-th column in the file. Note that you cannot skip columns! The user must provide enough variables to store each column. Each dataset handles a single basic data type.  That means different datasets must be created if one wishes to store different types.
ARGUMENT ai the i-th input array
ARGUMENT file name of file should end with .h5
ARGUMENT type optional the data type that is stored in h5. Supported datatypes are: "float{16,32,64}", "u?int{16,32,64}" (default: "float64")
ARGUMENT dataset optional name of dataset (default "data")
ARGUMENT gzip level of compression to use (default "0")
ARGUMENT chunksize optional number of arrays kept in one chunk (default "1"). If 0 the array is written contigiously.
ARGUMENT cachemb,cacheslots optional size of h5 cache (default: cachemb=64, cacheslots=12421).
OUTPUT file output filename
END_DESCRIPTION
*/
void WriteH5(char *in)
{
        int narr, arrlen;
        char *word;
        double **data;

        if (NULL == (word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
        if (!GetArg(in, "file", word)) {
                Warning("Error: file argument is missing!\n");
                goto efile;
        }
        printf("Writing to %s\n", word);
        struct h5io* io = h5io_init(word);
        if (NULL == io) goto eio;

        if (getarrays(in, &data, &arrlen, &narr)) goto edata;
        if (!GetOption(in, "dataset", word)) snprintf(word, strlen(in), "data");
        h5io_setdataset(io, word);

        if (h5io_isin(io)) {
                Warning("Error: dataset %s exists!\n", io->dataset);
                goto eexist;
        }

        if (!GetOption(in, "type", word)) snprintf(word, strlen(in), "float64");
        h5io_setdtype(io, word);
        if (FetchOptInt(in, "gzip", word, &io->compression)) io->compression = 0;
        if (FetchOptInt(in, "chunksize", word, &io->chunkarr)) io->chunkarr = 1;
        if (FetchOptInt(in, "cachemb", word, &io->cachemb)) io->cachemb = 64;
        if (FetchOptInt(in, "cacheslots", word, &io->cacheslots)) io->cacheslots = 12421;

        TIC();
        if (h5io_write(io, data, arrlen, narr)) goto ewrite;
        printf("Wrote %s in %g s\n", io->dataset, TOC());

        h5io_free(io);
        free(data);
        free(word);
        return;
ewrite:
eexist:
        free(data);
edata:
        h5io_free(io);
eio:
efile:
        free(word);
eword:
        Warning("Error: write_h5 failed!\n");
        return;
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG read_h5 ReadH5 "file=<file-str> a0=<out-array> a1=<out-array> .. aN=<out-array> [dataset=<str>]"
DESCRIPTION Read a dataset from an HDF5-file and stores the columns in arrays. The i-th column is stored in the i-th output array. Every array is converted to a 64-bit float for internal processing. Note that you cannot skip columns!
ARGUMENT file input filename
ARGUMENT dataset optional name of the dataset to read from (default: "data")
OUTPUT ai the i-th output array
END_DESCRIPTION
*/
void ReadH5(char *in)
{
        int i, arrlen, narr;
        double **data;
        char *word, **names;

        if (NULL == (word=malloc((strlen(in)+1)*sizeof(*word)))) goto eword;
        if (!GetArg(in, "file", word)) {
                Warning("Error: file argument is missing!\n");
                goto efile;
        }
        printf("Reading from %s\n", word);

        struct h5io* io = h5io_init(word);
        if (NULL == io) goto eio;

        if (!GetOption(in, "dataset", word)) snprintf(word, strlen(in), "data");
        h5io_setdataset(io, word);
        if (getnames(in, &names, &narr)) goto enames;

        TIC();
        if (h5io_read(io, &data, &arrlen, narr)) goto eread;
        printf("Read %s in %g s\n", io->dataset, TOC());

        for (i=0; i < narr; ++i)
                if (addarray(data[i], arrlen, names[i])) goto eadd;

        free(names);
        free(data);
        h5io_free(io);
        free(word);
        return;
eadd:
        // cleanup data for failed array
        for (int j=i; j < narr; ++j) {
                free(data[j]);
                free(names[j]);
        }
        free(data);
eread:
        free(names);
enames:
        h5io_free(io);
eio:
efile:
        free(word);
eword:
        Warning("Error: read_h5 failed!\n");
        return;
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG write_png WritePng "z=<in-array> nx=<int-value> ny=<int-value> ofn=<file> [normalise=<int-value>]"
DESCRIPTION Writes 2d array to a png image. The output file contains 3 channels, though only grayscale images are possible to write.
ARGUMENT z the 2d array
ARGUMENT nx the width of the image
ARGUMENT ny the height of the image
ARGUMENT ofn name of the output filename
ARGUMENT normalise optional argument. If 0 no normalisation is performed (negative values turn to 0, values larger than 255 turn to 255).
OUTPUT file output filename
END_DESCRIPTION
*/
void WritePng(char *in)
{
        int norm, nx, ny;
        char *word, *ofn;
        array *z;
        enum normalisation normalise;

        ofn = malloc((strlen(in)+1)*sizeof(*ofn));
        if (NULL == ofn) goto eofn;
        if (!GetArg(in, "ofn", ofn)) goto eofnmissing;

        word = malloc((strlen(in)+1)*sizeof(*word));
        if (NULL == word) goto eword;

        if (FetchArray(in, "z", word, &z)) goto efetch;
        if (FetchInt(in, "nx", word, &nx)) goto efetch;
        if (FetchInt(in, "ny", word, &ny)) goto efetch;
        if (FetchOptInt(in, "normalise", word, &norm))
                norm = 1;

        if (0==norm)
                normalise = NORM_NONE;
        else
                normalise = NORM_MAXMIN;

        if (z->N != nx * ny) {
                Warning("length(z) != nx*ny!\n");
                goto enxny;
        }

        TIC();
        if (write_png((const char*)ofn, z->D, nx, ny, normalise))
                Warning("Failed writing the png file!\n");
        printf("Wrote %s in %g s\n",ofn,TOC());
		
enxny:
efetch:
        free(word);
eword:
eofnmissing:
        free(ofn);
eofn:
        return;
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG make_array MakeArray "x=<out-array> x1=<float-value> x2=<float-value> Nx=<int-value>"
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
PARSEFLAG make_grid MakeGrid "x=<out-array> y=<out-array> x1=<float-value> x2=<float-value> y1=<float-value> y2=<float-value> Nx=<int-value> Ny=<int-value>"
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
PARSEFLAG get_grid GetGrid "C=<in-config> x=<out-array> y=<out-array> [nx=<out-1dim-array> ny=<out-1dim-array>]"
DESCRIPTION Extract the grid from the configured topogrid
ARGUMENT C  Simulation config with a configures topogrid
OUTPUT x output array
OUTPUT y output array
OUTPUT nx optional output
OUTPUT ny optional output
END_DESCRIPTION
*/
void GetGrid(char *in)
{
	char *word;
	array x, y, nx, ny;
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

    word=malloc((strlen(in)+1)*sizeof(char));
    if (GetOption(in, "nx", word)) {
            nx.N=1;
            nx.D=malloc(sizeof(*nx.D));
            nx.D[0] = C->Tx.Nx;

            if(AddArray(word, nx)) {
                    free(nx.D);
                    free(word);
            }
    }

    word=malloc((strlen(in)+1)*sizeof(char));
    if (GetOption(in, "ny", word)) {
            ny.N=1;
            ny.D=malloc(sizeof(*ny.D));
            ny.D[0] = C->Tx.Ny;

            if(AddArray(word, ny)) {
                    free(ny.D);
                    free(word);
            }
    }

	return;
}
/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG make_scalar MakeScalar "x=<out-array> val=<float-value>"
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
PARSEFLAG perturb Perturb "x=<in/out-array> [releps=<float-value>] [abseps=<float-value>]"
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
	double releps = 0.0;
	double abseps = 0.0;
	int i;
	array *x;
	word=malloc((strlen(in)+1)*sizeof(char));
	if (FetchOptFloat(in, "releps", word, &releps))
            releps = 0.0;
	releps=fabs(releps);

	if (FetchOptFloat(in, "abseps", word, &abseps))
            abseps = 0.0;
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
