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

void init_h5interface(void);

void free_h5interface(void);


struct supported_type map_supported_types_to_h5types(char* type);
