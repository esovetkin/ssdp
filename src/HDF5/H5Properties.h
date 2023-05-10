#include <hdf5/serial/hdf5.h>

/*
    Create a property list for writing datasets
    args:
        chunk_ndims: number of dimens the chunk this must be the length of chunk_dims
        chunk_dims: the size of the chunk in each dim
    
    e.g. if the dataset is a nrows x ncols matix and you wish to chunk together 2 columns and 5 rows
    you choose chunk_ndims = 2 and chunk_dims = {}
*/
hid_t H5P_create_dataset_proplist(int chunk_ndims, const hsize_t *chunk_dims);

hid_t H5P_create_16_MB_Chunk_Cache_access_dataset_proplist();

/*
TODO
I do not know how to solve it right currently
Jenya is not happy with this but the HDF5 API function H5Pset_chunk_cache is too cluncky
for later optimizations of performance so if the need arrives we have to start looking around
here for ideas how to fix it

Create a dataset access property list which has
the chunk cache set to be chunk_cache_size bytes large.

args:
    chunk_size: size in bytes of the chunks of the dataset
    chunk_cache_size: size in bytes of the chunk_cache. The bigger the better the performance and the higher the memory demand
    performance_val:
                    Value between 0 and 1 this will influence the size of the chunk slats in the chunk cache
                    If it is 0 the performance will be the smallest but the memory demand will be smaller
                    If is 1 the performance will be highest but the memory demand will be highe
                    The chunk slots are chose in a range according to the documentation
                    found here https://docs.hdfgroup.org/hdf5/v1_14/group___d_a_p_l.html#ga104d00442c31714ee073dee518f661f1
return:
    hid_t of the dataset access property list if sucessfull otherwise H5I_INVALID_HID
*/
//hid_t H5P_create_access_dataset_proplist(size_t chunk_size, size_t chunk_cache_size, double performance_val);