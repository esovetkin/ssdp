#include <stdlib.h>
#include <math.h>
#include <hdf5/serial/hdf5.h>

#include "H5Enums.h"
#include "H5Dataset.h"
#include "H5Datatypes.h"
#include "H5Properties.h"

/*
    Create and initialize a H5DDatasetHandler struct.
    The struct must be freed!
    args:
        name: name of dataset should be like a unix path e.g. /foo/bar/datasetname
        loc: hid_t of the H5D file in which the dataset lives
    return:
        pointer to H5DatasetHandler struct
*/
struct H5DatasetHandler* H5DatasetHandler_init(const char *name, hid_t loc){
    struct H5DatasetHandler* self;
    self = malloc(sizeof(*self));
    if (NULL == self){
        return NULL;
    }
    self->loc = loc;
    self->name = name;
    self->read_data = NULL;
    self->read_ncols = -1;
    self->read_nrows = -1;
    return self;
}

/*
    Write a continuous 1d array of doubles which is interpreted as a matrix into a HDF5 file using the H5DatasetHandler helper struct
    args:
        self: Handler object created by H5DatasetHandler_init
        data: pointer to array of doubles
        nrows: rows of matrix
        ncols: columns of matrix
        disk_datatype: HDF5 Datatype to save the data on the disk
        chunk_size: number of rows making up a chunk for IO purposes
    return:
        SUCESS if operation worked otherwise an enum with a nonzero value
*/ 
ErrorCode H5DatasetHandler_write_array(struct H5DatasetHandler *self, double* data, int nrows, int ncols, hid_t disk_datatype, hsize_t chunk_size){
    /*
        TODOs:
            check if the values are in the correct range to fit into a 16 bit int
            add chunksize argument 
    */
    hsize_t dims[2];
    herr_t status;
    dims[0]=nrows;
    dims[1]=ncols;
    hid_t dataspace_id, dataset_id;
    hid_t dataspace_create_props, dataspace_access_props;
    dataspace_id = H5Screate_simple(2,dims,NULL);
    if(H5I_INVALID_HID == dataspace_id){
        return FAILURE;
    }
    const hsize_t chunk_dims[2] = {chunk_size,ncols};
    dataspace_create_props = H5P_create_dataset_proplist(2, chunk_dims);
    dataspace_access_props = H5P_create_16_MB_Chunk_Cache_access_dataset_proplist();
    dataset_id = H5Dcreate(self->loc, self->name, disk_datatype, dataspace_id, H5P_DEFAULT, dataspace_create_props, dataspace_access_props);
    if (H5I_INVALID_HID == dataset_id){
        H5Sclose(dataspace_id);
        H5Pclose(dataspace_create_props);
        return FAILURE;
    }
    
    status  = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id); // todo is it ok to close it here already?
    H5Pclose(dataspace_create_props);
    if (0 > status){
        return FAILURE;
    }
    else{
        return SUCCESS;
    }
}

/*
    Read a 2D Dataset with fixed sizes from a HDF5 file into the Handler struct
    This function allocated memory in H5DatasetHandler->read_data which needs to be freed!
    args:
        self: pointer returned by H5DatasetHandler_init
    return:
        SUCESS if operation worked otherwise an enum with a nonzero value
*/
ErrorCode H5DatasetHandler_read_array(struct H5DatasetHandler *self){
    /*
        todos
    */

    ErrorCode return_val = SUCCESS;
    hid_t dataspace_id, dataset_id;
    herr_t status;
    hsize_t dims[2];
    hsize_t maxdims[2];
    int n_dims;

    dataset_id = H5Dopen(self->loc, self->name, H5P_DEFAULT);
    if (H5I_INVALID_HID == dataset_id){
        return_val = FAILURE;
        goto error;
    }
    dataspace_id = H5Dget_space(dataset_id);
    if(H5I_INVALID_HID == dataspace_id){
        H5Dclose(dataset_id);
        return_val = FAILURE;
        goto error;
    }
   
    // Find rank and retrieve current and maximum dimension sizes.
    n_dims = H5Sget_simple_extent_dims(dataspace_id, dims, maxdims);
    if (n_dims < 0){
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        return_val = FAILURE;
        goto error;
    }

    self->read_nrows = dims[0];
    self->read_ncols = dims[1];
    
    self->read_data = malloc(sizeof(double)*dims[0]*dims[1]);
    
    if (NULL == self->read_data){
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        return_val = OUTOFMEMORY;
        goto error;
    }

    // todo create array out_data and fill from h5 file
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, self->read_data);

    if(status < 0){
        free(self->read_data);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        return_val = FAILURE;
        goto error;
    }

    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    return return_val;

error:
    self->read_data = NULL;
    self->read_nrows = -1;
    self->read_ncols = -1;
    return return_val;

}

/*
    we dont have this free because the memory is passed to the user and therfore we dont really need a free function anymore
    void H5DatasetHandler_free(struct H5DatasetHandler **self_addr){
        struct H5DatasetHandler *self = *self_addr;
        //if(NULL != self->read_data){
        //    free(self->read_data);
        //}
        free(self);
        *self_addr = NULL;
    }
*/