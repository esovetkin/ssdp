#pragma once
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <hdf5/serial/hdf5.h>
#include "H5Enums.h"

struct H5DatasetHandler{
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
    const char *name;
    hid_t loc;
};

/*
    Create and initialize a H5DDatasetHandler struct.
    The struct must be freed!
    args:
        name: name of dataset should be like a unix path e.g. /foo/bar/datasetname
        loc: hid_t of the H5D file in which the dataset lives
    return:
        pointer to H5DatasetHandler struct
*/
struct H5DatasetHandler* H5DatasetHandler_init(const char *name, hid_t loc);

/*
    Write a continuous 1d array of doubles which is interpreted as a matrix into a HDF5 file using the H5DatasetHandler helper struct
    args:
        self: Handler object created by H5DatasetHandler_init
        data: pointer to array of doubles
        nrows: rows of matrix
        ncols: columns of matrix
        chunk_size: number of rows making up a chunk for IO purposes
    return:
        SUCESS if operation worked otherwise an enum with a nonzero value
*/ 
ErrorCode H5DatasetHandler_write_array(struct H5DatasetHandler *self, double* data, int nrows, int ncols, hid_t disk_datatype, hsize_t chunk_size);

/*
    Read a 2D Dataset with fixed sizes from a HDF5 file into the Handler struct
    This function allocated memory in H5DatasetHandler->read_data which needs to be freed!
    args:
        self: pointer returned by H5DatasetHandler_init
    return:
        SUCESS if operation worked otherwise an enum with a nonzero value
*/
ErrorCode H5DatasetHandler_read_array(struct H5DatasetHandler *self);