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
#include "h5interface.h"

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
	Transpose a 2d array into a 1d array.
	This function allocates memory that must be freed.

	args:
		arr: 2d array of size ncols x nrows
		nrows: number of rows in arr
		ncols: number of columns in arr
	return:
		pointer to transposed and unraveled/flattened array or NULL if this function fails
*/
double* transpose_unravel(double** arr, int nrows, int ncols) {
    // Allocate memory for the transposed array
    double* transposed = (double*)malloc(nrows * ncols * sizeof(double));
	if(NULL == transposed){
		return NULL;
	}
    // Transpose the array
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            double yeet = arr[j][i];
			transposed[i * ncols + j] = yeet;
        }
    }
    return transposed;
}

/*
private helper function to find the names of array variables from a given input

args:
	in: user input
	n_arr_vars: address of int variable to store the number of name variables found
return:
	array of strings of length n_arr_vars if successful else NULL
*/
char ** get_array_names_from_input(char *in, int *n_arr_vars){
	char **names = NULL;
	char *word = NULL;
	int Na=4;
	*n_arr_vars = 0;
	
	
	names=malloc(Na*sizeof(char *));
	if (NULL == names){
		Warning("Error: Out of malloc memory");
		goto error;
	}
	word = malloc((strlen(in)+1)*sizeof(char));
	if (NULL == word){
		Warning("Error: Out of malloc memory");
		goto error;
	}
	while(GetNumOption(in, "a", *n_arr_vars, word))
	{	
		names[*n_arr_vars]=word;
		word = malloc((strlen(in)+1)*sizeof(char));
		if (NULL == word){
			Warning("Error: Out of malloc memory");
			goto error;
		}
		*n_arr_vars=(*n_arr_vars)+1;
		if ((*n_arr_vars)==Na-1)
		{
			Na+=4;
			char **tmp;
            tmp=realloc(names, Na*sizeof(char *));
            if(NULL == tmp){
                Warning("Error: Out of malloc memory");
                goto error;
            }else {
                names = tmp;
            }
		}
	};
	free(word);
	if ((*n_arr_vars)==0)
	{
		Warning("Cannot define arrays from file, no array arguments recognized\n"); 
		goto error;
	}
	return names;
	
error:
	for(int i = 0; i < *n_arr_vars; i++){
		free(names[i]);
	}
	free(names);
	return NULL;
}

/*
	private helper function to create array variables from an array of columns
	args:
		ncols: number of columns
		nrows: number of rows
		names: array of strings
		data: flattened 2d array of shape nrows x ncols
	return:
		SUCCESS if successful
*/
static ErrorCode create_arrays(int ncols, int nrows, char **names, double **data){
	ErrorCode out = SUCCESS;
	array *arrays = NULL;
	arrays = malloc(sizeof(*arrays)*ncols);
	if(NULL == arrays){
		Warning("Failed to create arrays. Out of malloc memory");
		out = OUTOFMEMORY;
		goto end;
	}
	// create read_ncols arrays
	for(int i = 0; i < ncols; i++){
		arrays[i].N = nrows;
		arrays[i].D = data[i];
	}
	

	for(int i = 0; i < ncols; i++){
		printf("Creating array %s\n", names[i]);
		if(AddArray(names[i], arrays[i]))
		{	
			Warning("Failed to create array variable");
			free(names[i]);
			free(arrays[i].D);
			out = FAILURE;
		}
	}
end:
	free(arrays);
	return out;
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG read_h5 ReadArraysFromH5 "file=<file-str> a0=<out-array> a1=<out-array> .. aN=<out-array> [dataset=<str>]"
DESCRIPTION Reads a dataset called `dataset` from an HDF5-file and stores the columns in arrays. The i-th column is stored in the i-th output array. Every array is converted to a 64 bit float for internal processing. Note that you cannot skip columns! This command makes use of the HDF5 file-pool see `flush_h5` for more informations.
ARGUMENT file input filename
ARGUMENT dataset optional name of the dataset to read from (default: "data")
OUTPUT ai the i-th output array
END_DESCRIPTION
*/
void ReadArraysFromH5(char *in)
{
	char **names = NULL;
	int n_arr_vars = 0;
	char *file = NULL;
	double **data = NULL;
	char *dataset_name = NULL;
	int read_nrows, read_ncols;
	struct H5FileIOHandler *handler = NULL;
	ErrorCode err = SUCCESS;
	file=malloc((strlen(in)+1)*sizeof(char));
	if(NULL == file) {
        Warning("Out of malloc memory");
        goto end;
    }
	if (!GetArg(in, "file", file)) {
		goto end;
	}

	dataset_name=malloc((strlen(in)+1)*sizeof(char));
	if (!GetOption(in, "dataset", dataset_name)) {
		snprintf(dataset_name, strlen(in)+1, "data");
		Warning("No argument `dataset` provided! Using default value: `%s`\n", dataset_name);
	}

	names = get_array_names_from_input(in, &n_arr_vars);
	if (NULL == names) {
		goto end;
	}

	handler = H5FileIOHandlerPool_get_handler(g_h5filepool, file, R);
	if (NULL == handler) {
		Warning("Error reading H5 file: Could not create handler.!\n");
		goto end;
	}
	err = H5FileIOHandler_read_array_of_columns(handler, dataset_name, &data, &read_nrows, &read_ncols);
	if(SUCCESS != err) {
		Warning("Error reading H5 file could not read dataset\n");
		goto end;
	}

	if(n_arr_vars != read_ncols) {
		// user didn't provide enough variables
		Warning("Not enough variables provided to store arrays! Provided Variables: %d Read Arrays: %d\n",
		n_arr_vars, read_ncols);
		goto end;
	}

	if (read_ncols>0) {
		// this function frees memory if it failed to create a variable
		// thus it has its own goto tag
		err = create_arrays(read_ncols, read_nrows, names, data);
		goto end_create_arrays;
	}
end:
	if(SUCCESS != err) {
		Warning("An error while parsing the file has ocurred");
		for(int i = 0; i < n_arr_vars; i++){
			if(NULL != data)
				free(data[i]);
			if(NULL != names)
				free(names[i]);
		}
	}
end_create_arrays:
	free(file);
	free(names);
	free(data);
	free(dataset_name);
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG flush_h5 FlushH5 "[f0=<file-str>, f1=<file-str>, .., fN=<file-str>]"
DESCRIPTION This library handles all HDF5-files it reads and writes using a resource pool. Whenever an HDF5-file is created and opened for writing or opened for reading it is added to the pool. Files are opened until the end of the SSDP program or until they are flushed using this command. Only once a file is flushed its contents are stored on the hard drive. If no file-string is provides flush_h5 closes all currently opened files. Otherwise only the specified files are flushed. Flushing a HDF5-file does not affect variables created by reading the file. If one wishes to read an HDF5-file previously created by SSDP in the same script it must be flushed first. If a one wishes to write multiple datasets to a file it should only be flushed after all write_h5 calls.
END_DESCRIPTION
*/
void FlushH5(char *in){
	
	char *filename;
	filename=malloc((strlen(in)+1)*sizeof(char));
	int i = 0;
	while(GetNumOption(in, "f", i, filename)){
		i++;
		printf("Writing H5 file %s.\n", filename);
		H5FileIOHandlerPool_close_file(g_h5filepool, filename);
	}
	if(i == 0){
		printf("Writing all currently opened H5 files.\n");
		H5FileIOHandlerPool_close_all_files(g_h5filepool);
	}
	free(filename);
}


/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG write_h5 WriteArraysToH5 "a0=<in-array> a1=<in-array> .. aN=<in-array> file=<file-str> [type=<type-str>] [dataset=<str>] [chunksize=<int>]"
DESCRIPTION Writes arrays in columns of a HDF5-file in a dataset called 'dataset'. The i-th array is written to the i-th column in the file. Note that you cannot skip columns! The user must provide enough variables to store each column. Each dataset handles a single basic data type. That means if one wishes to store different types different datasets are created. This command makes use of the HDF5 file-pool see `flush_h5` for more informations.
ARGUMENT ai the i-th input array
ARGUMENT file name of file should end with .h5
ARGUMENT type optional string that describes the datatype to save on the disc supported datatypes are: float16, float64, int32, int64 (default: float64)
ARGUMENT dataset optional name of dataset (default "data")
ARGUMENT chunksize in number of rows for chunked I/O. Chunksize influence the efficiency of the compression (default 1000)
OUTPUT file output filename
END_DESCRIPTION
*/
void WriteArraysToH5(char *in)
{
	char *word = NULL;
	char *file = NULL;
	char *dataset_name = NULL;
	char *type_name = NULL;
	struct supported_type type;
	double **data = NULL;
	int chunk_size;
	array *a = NULL;
	int ncols, data_buffer_size=4, nrows=-1;
	file=malloc((strlen(in)+1)*sizeof(char));
	if (!GetArg(in, "file", file)) goto error;

	dataset_name=malloc((strlen(in)+1)*sizeof(char));
	if (!GetOption(in, "dataset", dataset_name)){
		snprintf(dataset_name, strlen(in)+1, "data");
		Warning("No argument `dataset` provided! Using default value: `%s`\n", dataset_name);
	};

	type_name=malloc((strlen(in)+1)*sizeof(char));
	if (!GetOption(in, "type", type_name)){
		snprintf(type_name, strlen(in)+1, "float64");
		Warning("No argument `type` provided! Using default value: `%s`\n", type_name);
	};

	type = map_supported_types_to_h5types(type_name);
	if(H5I_INVALID_HID == type.type_id) {
		Warning("Invalid datatype type: %s\n", type_name);
		goto error;
	}
	word=malloc((strlen(in)+1)*sizeof(char));
	if(NULL == word) {
		Warning("Error: Out of malloc memory!");
		goto error;
	}
	if(FetchOptInt(in, "chunksize", word, &chunk_size)) {
		chunk_size = 1000;
		Warning("No argument `chunksize` provided! Using default value %d."
		" If HDF5 IO performance is poor consider increasing the chunk size or the size of the chunk cache\n", chunk_size);
		// TODO we may need an API call for the chunk cache	in here maybe as an optional argument
	}


	ncols=0;
	data=malloc(data_buffer_size*sizeof(*data));
	if(NULL == data) {
		Warning("Error: Out of malloc memory!");
		goto error;
	}

	while(GetNumOption(in, "a", ncols, word)) {
		if (!LookupArray(word, &a)) {
			Warning("Array %s not available\n",word);
			goto error;
		}

		data[ncols]=a->D;
		if (nrows<0) {
			nrows=a->N;
		} else if (a->N!=nrows) {
			printf("0: %d %d:%d\n", nrows, ncols, a->N);
			Warning("Error: arrays must be of equal length\n");
			goto error;
		}
		ncols++;
		if (ncols==data_buffer_size-1)
		{
			data_buffer_size+=4;
			double **tmp;
			tmp=realloc(data, data_buffer_size*sizeof(*data));
			if(NULL == tmp) {
				Warning("Error: Out of malloc memory!");
				goto error;
			} else {
				data = tmp;
			}

		}
	}
	if (ncols==0)
	{
		Warning("Cannot define arrays from file, no array arguments recognized\n");
		goto error;
	}
	printf("Preparing to write arrays to H5 file `%s` in dataset `%s`\n", file, dataset_name);
	struct H5FileIOHandler *handler = H5FileIOHandlerPool_get_handler(g_h5filepool, file, W);
	if (NULL == handler){
		Warning("Error creating HDF5 file!\n");
		goto error;
	}
	ErrorCode err;

	err = H5FileIOHandler_write_array_of_columns(handler, dataset_name, data, nrows, ncols, chunk_size, type.type_id);
	if (SUCCESS != err)
		Warning("Error writing to HDF5 file!\n");
error:
	free(file);
	free(data);
	free(type_name);
	free(dataset_name);
	free(word);
}

/*
BEGIN_DESCRIPTION
SECTION Array
PARSEFLAG write_png WritePng "z=<in-array> nx=<1-dim array> ny=<1-dim array> ofn=<file> [normalise=1]"
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
        int n;
        char *word, *ofn;
        array *z, *nx, *ny;
        enum normalisation normalise;

        ofn = malloc((strlen(in)+1)*sizeof(*ofn));
        if (NULL == ofn) goto eofn;
        if (!GetArg(in, "ofn", ofn)) goto eofnmissing;

        word = malloc((strlen(in)+1)*sizeof(*word));
        if (NULL == word) goto eword;

        if (FetchArray(in, "z", word, &z)) goto efetch;
        if (FetchArray(in, "nx", word, &nx)) goto efetch;
        if (FetchArray(in, "ny", word, &ny)) goto efetch;

        if (FetchOptInt(in, "normalise", word, &n))
                n = 1;

        if (0==n)
                normalise = NORM_NONE;
        else
                normalise = NORM_MAXMIN;

        if (nx->N != 1 || ny->N != 1) {
                Warning("nx and ny must be 1 dimensional arrays!\n");
                goto enxny;
        }

        int Nx = (int)(round(nx->D[0]));
        int Ny = (int)(round(ny->D[0]));
        if (z->N != Nx * Ny) {
                Warning("length(z) != nx*ny!\n");
                goto enxny;
        }

        if (write_png((const char*)ofn, z->D, Nx, Ny, normalise))
                Warning("Failed writing the png file!\n");

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
