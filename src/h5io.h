#ifndef _H5IO_H
#define _H5IO_H

#include <hdf5.h>

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
        char dtype[1024];
        char dataset[1024];
};

struct h5io* h5io_init(const char *filename);
void h5io_free(struct h5io* self);
void h5io_setdataset(struct h5io* self, const char *dataset);
void h5io_setdtype(struct h5io* self, const char *dtype);

/** Read *first* `narr` into &(**double)

    This allocated memory for the **double and sets the `arrlen` the
    length of each array.
 */
int h5io_read(struct h5io* self, double ***data, int* arrlen, int narr);


/** Write (**double) to the specified dataset

    len( data) == narr
    len(*data) == arrlen
 */
int h5io_write(struct h5io* self, double **data, int arrlen, int narr);


/**
    Returns:
    - 0 path does not exists
    - 1 path exists
    - -1 path is invalid (e.g. /path/1 is invalid, wheh /path is a dataset)
    - -2 malloc fails
*/
int h5io_isin(struct h5io* self);

#endif
