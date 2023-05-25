#include "h5interface.h"
#include "string.h"
#include <stdlib.h>
struct H5FileIOHandlerPool *g_h5filepool = NULL;
# define N_SUPPORTED_TYPES 4
struct supported_type *g_supported_h5types;

struct supported_type build_supported_type(char *type_str, hid_t type_id, int is_custom){
    struct supported_type out = {type_str, type_id, is_custom};
    return out;
}

void init_h5interface(){
    g_h5filepool = H5FileIOHandlerPool_init();
    g_supported_h5types = malloc(sizeof(struct supported_type)*N_SUPPORTED_TYPES);
    g_supported_h5types [0] = build_supported_type("float16", H5T_define_16bit_float(), 1);
    g_supported_h5types [1] = build_supported_type("float64", H5T_NATIVE_DOUBLE, 0);
    g_supported_h5types [2] = build_supported_type("int32", H5T_NATIVE_INT32, 0);
    g_supported_h5types [3] = build_supported_type("int64", H5T_NATIVE_INT64, 0);
}


void free_h5interface(){
    for(int i = 0; i < N_SUPPORTED_TYPES; i++){
        if(g_supported_h5types[i].is_custom){
            H5Tclose(g_supported_h5types->type_id);
        }
    }
    H5FileIOHandlerPool_free(&g_h5filepool);
}



struct supported_type map_supported_types_to_h5types(char* type){
	int return_idx = -1;
    for(int i = 0; i < N_SUPPORTED_TYPES; i++){
		if (strncmp(type, g_supported_h5types[i].type_str, 255) == 0){
			return_idx = i;
		}
	}
	// to avoid memory leaks close all custom types which are not returned
	// if valgrind says 1,848 bytes in 1 blocks are still reachable
	// caused by malloc used in H5E_get_stack this has nothing to do with this function
	// but with H5 as this leak also happens when H5T_define_16bit_float() is not called in this function
	//(I tested it by commenting it out) 
	/*
    for(int i = 0; i < N_SUPPORTED_TYPES; i++){
		if (i != return_idx && supported_types[i].is_custom){
			printf("Closed custom type called %s\n", supported_types[i].type_str);
			H5Tclose(supported_types[i].type_id);
		}
	}
    */
	if(return_idx >=0){
		return g_supported_h5types[return_idx];
	}
	struct supported_type error_out = {type, H5I_INVALID_HID, 0};
	return error_out;
}