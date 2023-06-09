#include "h5interface.h"
#include "string.h"
#include <stdlib.h>
struct H5FileIOHandlerPool *g_h5filepool = NULL;
# define N_SUPPORTED_TYPES 4
struct supported_type *g_supported_h5types;
// QTODO: if you choose to use arrays, then the following is a valid syntax
//
//         struct supported_type g_supported_h5types[N_SUPPORTED_TYPES];
//
// since you need the length of the array.

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
    // QTODO: I decided to opt for the catching sitation when
    // something fails, because if a catastroph happens you would like
    // to know that it happens, instead of getting some weird bug
    // somewhere in the middle.
    if (NULL == g_h5filepool) goto einit;

    // QTODO: your N_SUPPORTED_TYPES is hardcoded, why not just use
    // standard array? i.e.
    //    struct supported_type g_supported_h5types[4];
    g_supported_h5types = malloc(N_SUPPORTED_TYPES * sizeof(*g_supported_h5types));
    if (NULL == g_supported_h5types) goto einittypes;

    g_supported_h5types [0] = build_supported_type("float16", H5T_define_16bit_float(), 1);
    g_supported_h5types [1] = build_supported_type("float64", H5T_NATIVE_DOUBLE, 0);
    g_supported_h5types [2] = build_supported_type("int32", H5T_NATIVE_INT32, 0);
    g_supported_h5types [3] = build_supported_type("int64", H5T_NATIVE_INT64, 0);

    return 0;
    // free(g_supported_h5types);
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
                // QTODO: you check the ith is custom, but close the
                // first. It is something out of ordinary.
            H5Tclose(g_supported_h5types->type_id);
        }
    }
    // QTODO: you do malloc of this but I cannot see free. Free must
    // be somewhere here.
    free(g_supported_h5types);
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
