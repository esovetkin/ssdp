#include "HDF5/hd5extention/src/H5FileIO.h"
#include "HDF5/hd5extention/src/H5Datatypes.h"
#include "HDF5/hd5extention/src/H5Enums.h"
#include <stdlib.h>

struct supported_type {
	char *type_str;
	hid_t type_id;
	int is_custom;
};

extern struct H5FileIOHandlerPool* g_h5filepool;
extern struct supported_type *g_supported_h5types;

/*
    This function initializes the H5FileIOHandlerPool and all Datatypes via the supported_type struct.
    It must be called at the beginning of the ssdp main.

    The function returns non-zero values when something is wrong.
*/
int init_h5interface(void);

/*
    This function frees all the resources allocated by init_h5interface.
    It must be called at the end of ssdp main.
*/
void free_h5interface(void);

/*
    This function maps a string to a currently supported type.
    A type is supported if `type` matches to a `type_str` member of an object found in g_supported_h5types.
    If no matching type if found return a struct whose type_id member is H5I_INVALID_HID.
*/
struct supported_type map_supported_types_to_h5types(char* type);
