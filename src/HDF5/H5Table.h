#pragma once
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <hdf5/serial/hdf5.h>
#include "H5Enums.h"

struct H5TableHandler{
    /* 
        TODO
        checkout how to do 16 bit float and use it instead of uint16_t
        float16 offers more flexibility e.g. also store coordinates in radians
        make datatype and digits optional -> Enum?
        main datatypes 16, 32 and 64 bit
    */
    time_t t;
    hsize_t read_nrows;
    hsize_t read_ncols;
    double *read_data;
    char **read_columnnames;
    hid_t datatype;
    double digit_scale;
    const char *name;
    hid_t loc;
};


/*
    Create and initialize a H5TableHandler struct.
    The struct must be freed!
    args:
        name: name of dataset should be like a unix path e.g. /foo/bar/datasetname
        loc: hid_t of the H5D file in which the dataset lives
    return:
        pointer to H5TableHandler struct
*/
struct H5TableHandler* H5TableHandler_init(const char *name, hid_t loc);

/*
    Write a continuous 1d array of doubles which is interpreted as a table into a HDF5 file using the H5TableHandler helper struct
    args:
        self: Handler object created by H5TableHandler_init
        data: pointer to array of doubles
        nrows: rows of table
        ncols: columns of table
        names: array of ncols pointers to \0 terminated strings which
        chunk_size: size of chunks for compressed io
    return:
        SUCESS if operation worked otherwise an enum with a nonzero value
*/ 
ErrorCode H5TableHandler_write_table(struct H5TableHandler *self, double* data, int nrows, int ncols, const char **names, hsize_t chunk_size);

/*
    Read a 2D Dataset with fixed sizes from a HDF5 file into the Handler struct
    This function allocated memory in H5TableHandler->read_data and H5Tablehandler->read_columnnames which needs to be freed!
    args:
        self: pointer returned by H5TableHandler_init
    return:
        SUCESS if operation worked otherwise an enum with a nonzero value
*/
ErrorCode H5TableHandler_read_table(struct H5TableHandler *self);

