#pragma once
#include <hdf5/serial/hdf5.h>
#include "H5Enums.h"

/*
    This struct helps managing a HD5Dataset which holds a n \times m array
    The data is stored using 16-bit int to save diskspace => PRECISION MAY BE LOST

    create it using int
    write to file using write_array
    read from file using read_array
    when read read_X variables contain the data
*/
struct H5FileIOHandler{
    hid_t file_id;
    const char* filename;
    IOMode mode;
};


/*  
    Create a H5FileIOHandler which makes it easier to open and create HDF5 files.
    This functions allocates memory on the heap and opens resources. Free the struct using H5FileIOHandler_free!

    args:
        fn: name of the file to be opened or created
        mode: wether to create a new file, open a file for reading or open a file for editing
    return:
        pointer to initialized H5FileIOHandler struct

*/
struct H5FileIOHandler* H5FileIOHandler_init(const char* fn, IOMode mode);

/*
    Free the memory allocated and the resources opended by H5FileIOHandler_init
    args:
        self: pointer to H5FileIOHandler struct created by H5FileIOHandler_init
*/
void H5FileIOHandler_free(struct H5FileIOHandler **self_addr);

/*
    Read a 1D array of doubles which should be treated as a 2D Matrix from the dataset called dataset_name from the HDF5 file.
    This function allocates memory for it's output which must be freed by the user!
    args:
        self: pointer created by H5FileIOHandler_init
        dataset_name: name of the dataset in the HDF5 file
        out_data: pointer output array value
        out_nrows: number of rows of the Matrix
        out_ncols: number of columns of the Matrix
    return:
        SUCESS if it worked else a nonzero enum value
    TODO
        - check if self in in read mode if not return a non success errorcode
 */
ErrorCode H5FileIOHander_read_array(struct H5FileIOHandler *self, const char *dataset_name, double **out_data, int* out_nrows, int* out_ncols);

/*
    Write a 1D array of doubles which should be treated as a 2D Matrix from to a hd5 file in as a data set named dataset_name.
    
    args:
        self: pointer created by H5FileIOHandler_init
        dataset_name: name of the dataset in the HDF5 file
        data: pointer pointer to data
        nrows: number of rows of the Matrix
        ncols: number of columns of the Matrix
        chunk_size: number of rows making up a chunk for IO purposes
        disk_datatype: H5 Datatype which will be saved in the file
    return:
        SUCESS if it worked else a nonzero enum value
 */
ErrorCode H5FileIOHandler_write_array(struct H5FileIOHandler *self, const char *dataset_name, double *data, int nrows, int ncols, hsize_t chunk_size, hid_t disk_datatype);

/*
    Read a 1D array of doubles which should be treated as a table from the dataset called dataset_name from the HDF5 file.
    This function allocates memory for it's output which must be freed by the user!
    args:
        self: pointer created by H5FileIOHandler_init
        dataset_name: name of the dataset in the HDF5 file
        out_data: pointer to output array
        out_nrows: number of rows of the Matrix
        out_ncols: number of columns of the Matrix
        out_column_names: pointer to array of column names
    return:
        SUCESS if it worked else a nonzero enum value
    TODO
        - check if self in in read mode if not return a non success errorcode
 */
ErrorCode H5FileIOHander_read_table(struct H5FileIOHandler *self, const char *dataset_name, double **out_data, int* out_nrows, int* out_ncols, char*** out_column_names);

/*
    Write a 1D array of doubles which should be treated as a Table to a hd5 file in as a dataset named dataset_name.
    
    args:
        self: pointer created by H5FileIOHandler_init
        dataset_name: name of the dataset in the HDF5 file
        data: pointer pointer to data
        nrows: number of rows of the Matrix
        ncols: number of columns of the Matrix
        column_names: array of column names
    return:
        SUCESS if it worked else a nonzero enum value
 */
ErrorCode H5FileIOHandler_write_table(struct H5FileIOHandler *self, const char *dataset_name, double *data, int nrows, int ncols, const char** columns_names, hsize_t chunk_size);

/*
    TODO
    - make first extention of ssdp in new branch
    - function to add columns
    - function to add records
    - set type for each column maybe use basic dataset and attributes for it

*/