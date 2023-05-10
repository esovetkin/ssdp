#include "H5Datatypes.h"
#include "H5Enums.h"

hid_t H5T_define_nbit_float(hid_t base_datatype, size_t sign_pos, size_t exponent_pos, size_t exponent_size, size_t mantissa_pos, size_t mantissa_size){
    hsize_t base_size = H5Tget_size(base_datatype);
    hid_t new_datatype;
    herr_t status;
    size_t ebias = 1;
    if((base_size == 0) || (base_size > 1 + exponent_size + mantissa_size)){
        return H5I_INVALID_HID;
    }
    new_datatype = H5Tcopy(base_datatype);

    status = H5Tset_fields(new_datatype, sign_pos, exponent_pos, exponent_size, mantissa_pos, mantissa_size);
    if(status<0){
        goto error;
    }

    status = H5Tset_offset(new_datatype, mantissa_pos);
    if(status<0){
        goto error;
    }

    status = H5Tset_precision(new_datatype, exponent_pos);
    if(status<0){
        goto error;
    }

    status = H5Tset_size(new_datatype, base_size);
    if(status<0){
        goto error;
    }
    
    for(size_t i = 0; i<exponent_size-1; i++){
        ebias*=2;
    }
    ebias-=1;
    
    status = H5Tset_ebias(new_datatype, ebias);
    if(status<0){
        goto error;
    }

    return new_datatype;

error:
    H5Tclose(new_datatype);
    return H5I_INVALID_HID;

}
