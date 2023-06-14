#include "h5interface.h"
#include "string.h"
#include <stdlib.h>
struct H5FileIOHandlerPool *g_h5filepool = NULL;
# define N_SUPPORTED_TYPES 4
struct supported_type g_supported_h5types[N_SUPPORTED_TYPES];

struct supported_type build_supported_type(char *type_str, hid_t type_id, int is_custom){
    struct supported_type out = {type_str, type_id, is_custom};
    return out;
}
/*
    This function initializes the H5FileIOHandlerPool and all Datatypes via the supported_type struct.
    It must be called at the beginning of the ssdp main.
*/
int init_h5interface(){
    g_h5filepool = H5FileIOHandlerPool_init();
    // QTODO: I decided to opt for the catching situation when
    // something fails, because if a catastrophe happens you would like
    // to know that it happens, instead of getting some weird bug
    // somewhere in the middle.
    if (NULL == g_h5filepool) goto einit;

    if (NULL == g_supported_h5types) goto einittypes;

    g_supported_h5types [0] = build_supported_type("float16", H5T_define_16bit_float(), 1);
    g_supported_h5types [1] = build_supported_type("float64", H5T_NATIVE_DOUBLE, 0);
    g_supported_h5types [2] = build_supported_type("int32", H5T_NATIVE_INT32, 0);
    g_supported_h5types [3] = build_supported_type("int64", H5T_NATIVE_INT64, 0);

    return 0;
einittypes:
    free(g_h5filepool);
einit:
    return -1;
}

/*
    This function frees all the resources allocated by init_h5interface.
    It must be called at the end of ssdp main.
*/
void free_h5interface(){
    for(int i = 0; i < N_SUPPORTED_TYPES; i++){
        if(g_supported_h5types[i].is_custom){
            H5Tclose(g_supported_h5types[i].type_id);
        }
    }
    H5FileIOHandlerPool_free(&g_h5filepool);
}


/*
    This function maps a string to a currently supported type.
    A type is supported if `type` matches to a `type_str` member of an object found in g_supported_h5types.
    If no matching type if found return a struct whose type_id member is H5I_INVALID_HID.
*/
struct supported_type map_supported_types_to_h5types(char* type){
	int return_idx = -1;
    for(int i = 0; i < N_SUPPORTED_TYPES; i++){
		if (strncmp(type, g_supported_h5types[i].type_str, 255) == 0){
			return_idx = i;
		}
	}
	if(return_idx >=0){
		return g_supported_h5types[return_idx];
	}
	struct supported_type error_out = {type, H5I_INVALID_HID, 0};
	return error_out;
}
