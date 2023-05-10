#include <hdf5/serial/hdf5.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "H5Enums.h"
#include "H5Table.h"


/*
    Create and initialize a H5TableHandler struct.
    The struct must be freed!
    args:
        name: name of dataset should be like a unix path e.g. /foo/bar/datasetname
        loc: hid_t of the H5D file in which the dataset lives
    return:
        pointer to H5TableHandler struct
*/
struct H5TableHandler* H5TableHandler_init(const char *name, hid_t loc){
    struct H5TableHandler* self;
    self = malloc(sizeof(*self));
    if (NULL == self){
        return NULL;
    }
    self->loc = loc;
    self->name = name;
    self->read_data = NULL;
    self->read_ncols = -1;
    self->read_nrows = -1;
    self->read_columnnames = NULL;
    self->digit_scale = 10;
    self->datatype = H5T_NATIVE_UINT16; // values times 10 and round to int make sure it fits into range 2^16 - 1
    return self;
}

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
ErrorCode H5TableHandler_write_table(struct H5TableHandler *self, double* data, int nrows, int ncols, const char **names, hsize_t chunk_size){
    /*
        TODOs:
            what do we do if a dataset already exists?
            check if the values are in the correct range to fit into a 16 bit int
    */
    herr_t status;
    // if we know during compiletime the structure of the table
    // we can build structs to represent records of columns with different datatypes
    // https://portal.hdfgroup.org/display/HDF5/HDF5+Table+%28H5TB%29+Interface
    size_t field_offsets[ncols];
    hid_t field_types[ncols];
    int        *fill_data  = NULL;
    int         compress   = 1; // compression flag I think 0 is no compression
    if(chunk_size <= 0){
        return FAILURE;
    }
    for(int i = 0; i<ncols;i++){
        field_offsets[i]=sizeof(uint16_t)*i;
        field_types[i]=self->datatype;
    }
    
    uint16_t *rounded_data = malloc(sizeof(uint16_t)*ncols*nrows);
    if(NULL == rounded_data){
        return OUTOFMEMORY;
    }
    for(int i = 0; i<nrows; i++){
        for(int j = 0; j<ncols; j++){
            rounded_data[i*ncols+j] = (uint16_t) round((data[i*ncols+j]*self->digit_scale));
        }
    }
    
    status = H5TBmake_table(self->name, self->loc, self->name, ncols, nrows, sizeof(uint16_t)*ncols,
    names,field_offsets, field_types, chunk_size, fill_data, compress, rounded_data);
    free(rounded_data);
    if (0 > status){
        return FAILURE;
    }
    else{
        return SUCCESS;
    }
}

/*
    Read a 2D Dataset with fixed sizes from a HDF5 file into the Handler struct
    This function allocated memory in H5TableHandler->read_data which needs to be freed!
    args:
        self: pointer returned by H5TableHandler_init
    return:
        SUCESS if operation worked otherwise an enum with a nonzero value
*/
ErrorCode H5TableHandler_read_table(struct H5TableHandler *self){
    /*
        todos
    */

    ErrorCode return_val = SUCCESS;
    herr_t status;
    size_t *field_sizes = NULL;
    size_t *field_offsets = NULL;
    size_t *type_size = NULL;
    uint16_t *buff = NULL;

    status = H5TBget_table_info(self->loc, self->name, &(self->read_ncols), &(self->read_nrows));
    if (status < 0){
        return_val = FAILURE;
        goto error;
    }
    field_offsets = malloc(sizeof(size_t)*self->read_ncols);
    field_sizes = malloc(sizeof(size_t)*self->read_ncols);
    self->read_columnnames = malloc(sizeof(char*)*(self->read_ncols));
    for(hsize_t i=0; i<self->read_ncols; i++){
        self->read_columnnames[i] = malloc(sizeof(char)*100);
    }
    status = H5TBget_field_info(self->loc, self->name, self->read_columnnames, field_sizes, field_offsets, type_size);

    if (status < 0){
        return_val = FAILURE;
        goto error;
    }
    buff = malloc((self->read_nrows)*(self->read_ncols)*sizeof(uint16_t));
    if(NULL == buff){
        return_val = OUTOFMEMORY;
        goto error;
    }
    
    status = H5TBread_table(self->loc, self->name, sizeof(uint16_t)*(self->read_ncols), field_offsets, field_sizes, buff);
    

    self->read_data = malloc((self->read_nrows)*(self->read_ncols)*sizeof(double));
    if(NULL == self->read_data){
        return_val = OUTOFMEMORY;
        goto error;
    }
    
    for(hsize_t i = 0; i<self->read_nrows; i++){
        for (hsize_t j = 0; j<self->read_ncols; j++){
            self->read_data[i*(self->read_ncols) + j] = ((double) buff[i*(self->read_ncols) + j])/self->digit_scale;
        }
    }
    
    if (status < 0){
        return_val = FAILURE;
        goto error;
    }
    else{
        goto fin;
    }

error:
    self->read_data = NULL;
    self->read_nrows = -1;
    self->read_ncols = -1;
    self->read_columnnames = NULL;
    
fin:
    free(buff);
    free(field_sizes);
    free(field_offsets);
    free(type_size);
    return return_val;

}

/*
    we dont have this free because the memory is passed to the user and therfore we dont really need a free function anymore
    void H5DatasetHandler_free(struct H5DatasetHandler **self_addr){
        struct H5DatasetHandler *self = *self_addr;
        if(NULL != self->read_data){
            free(self->read_data);
        }
        free(self);
        *self_addr = NULL;
    }
*/