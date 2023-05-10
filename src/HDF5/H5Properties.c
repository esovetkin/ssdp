#include <hdf5/serial/hdf5.h>
#include <math.h>
#include <stdlib.h>
#include "H5Properties.h"

hid_t H5P_create_dataset_proplist(int chunk_ndims, const hsize_t *chunk_dims){
    herr_t status;
    hid_t dset_create_props;
    dset_create_props = H5Pcreate (H5P_DATASET_CREATE);

    status = H5Pset_chunk (dset_create_props, chunk_ndims, chunk_dims);
    if(status < 0) {
        goto error;
    }
    status = H5Pset_nbit(dset_create_props);
    if(status < 0) {
        goto error;
    }
    status = H5Pset_deflate(dset_create_props, 9);
    if(status < 0) {
        goto error;
    }
    
    return dset_create_props;
error:
    H5Pclose(dset_create_props);
    return H5I_INVALID_HID;
}
/*
linearly interpolate between start and stop
args:
    start: some value
    stop: some value grater than start
    factor: a value between 0 and 1
return:
    a value between start and stop if sucessfull otherwise NaN
*/
double linear_interpolate(double start, double stop, double factor){
    if(stop < start){
        return NAN;
    }
    if(factor <= 0 || factor >= 1){
        return NAN;
    }
    return start * (1.0 - factor) + (stop * factor);
}

hid_t H5P_create_16_MB_Chunk_Cache_access_dataset_proplist(){
    /*
    Example for creating a dataset acess plist with a 16 MB chunk cache
    taken from
    https://docs.hdfgroup.org/hdf5/v1_14/group___d_a_p_l.html#ga104d00442c31714ee073dee518f661f1
    */
    herr_t status;
    hid_t dset_access_props;
    dset_access_props = H5Pcreate (H5P_DATASET_ACCESS);
    status = H5Pset_chunk_cache(dset_access_props, 12421, 16*1024*1024, H5D_CHUNK_CACHE_W0_DEFAULT);
     if(status < 0) {
        goto error;
    }
    
    return dset_access_props;
error:
    H5Pclose(dset_access_props);
    return H5I_INVALID_HID;
}