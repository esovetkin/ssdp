#ifndef _H5IO_H
#define _H5IO_H

#include <hdf5.h>

#define H5IO_LNAME 1024

/**
   file: path to the h5 file

   dataset: name of the dataset to read/write (default: "data")

   dtype: dtype to write (default: "float64")

   compression: write compression level (0 no, max 9, default: 0)

   chunkarr: number of arrays in one chunk (default: 1). If 0 no
   chunking and compression is used

   cachemb: cache MB (default: 64)
*/
struct h5io {
		hid_t file;
		int compression;
		int chunkarr;
		int cachemb;
		int cacheslots;
		char dtype[H5IO_LNAME];
		char dataset[H5IO_LNAME];
};

struct h5io* h5io_init(const char *filename, int readonly);
void h5io_free(struct h5io* self);
void h5io_setdataset(struct h5io* self, const char *dataset);
void h5io_setdtype(struct h5io* self, const char* dtype);

/** Read *first* `narr` into &(**void)

	This allocates memory for the **data according to the size needed
	for `dtype` and sets the `arrlen` the length of each
	array. Specify the size of the void* with dtype.
 */
int h5io_read(struct h5io* self, void ***data, const char *dtype, int* arrlen, int narr);


/** Write (**double) to the specified dataset

	len( data) == narr
	len(*data) == arrlen
 */
int h5io_write(struct h5io* self, void **data, const char *dtype, int arrlen, int narr);


/**
	Returns:
	- 0 path does not exists
	- 1 path exists
	- -1 path is invalid (e.g. /path/1 is invalid, wheh /path is a dataset)
	- -2 malloc fails
*/
int h5io_isin(struct h5io* self);

hid_t h5io_fopen(const char *fn, int readonly);
int h5_datasetisin(hid_t file, const char* dst);

/** Set comment to the currently selected dataset.

	Return non-zero if failed.

 */
int h5io_comment(struct h5io* self, const char *cmmnt);

/** Check that size of double/int you provide have the same sizes of
 * dtype string. Return non-zero for non-equal sizes. For instance,
 *
 * h5io_checksize(sizeof(double), "float64");
 * h5io_checksize(sizeof(float), "float32");
 * h5io_checksize(sizeof(int), "int32");
 *
 */
int h5io_checksize(size_t s, const char* dtype);

#endif
