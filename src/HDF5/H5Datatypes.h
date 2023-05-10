#pragma once
#include <hdf5/serial/hdf5.h>
#include <stdlib.h>

/*
    This Module provides some helper functions to use custom HDF5 Datatypes
*/

/*
    Define a float datatype with an n-bit precision.
    The defined datatype must be closed with H5Tclose!

    In the following position referce to the position when counting from the least significant bit.
    args:
        base_datatype: H5 Datatype to base this type on. The new type will never be bigger than base_datatype
        sign_pos: position of the sign bit
        exponent_pos: position of the exponent
        exponent_size: number of bits in the exponent
        mantissa_pos: position of the mantissa
        mantissa_size: number of bits in the mantissa

    return:
        hid_t of the newly defined datatype or H5I_INVALID_HID if an error occurrs

    Example of a 20 bit float defined within 32 bits:
        size=4 byte, precision=20 bits, offset=7 bits,
        mantissa size=13 bits, mantissa position=7,
        exponent size=6 bits, exponent position=20,
        exponent bias=31.
        It can be illustrated in little-endian order as:
        (S - sign bit, E - exponent bit, M - mantissa bit, ? - padding bit)
        3       2        1        0
        ?????SEE EEEEMMMM MMMMMMMM M???????
*/
hid_t H5T_define_nbit_float(hid_t base_datatype, size_t sign_pos, size_t exponent_pos, size_t exponent_size, size_t mantissa_pos, size_t mantissa_size);

// TODO make this a macro for 16 bit
# define H5T_define_16bit_float() H5T_define_nbit_float(H5T_NATIVE_DOUBLE, 63, 58, 5, 48, 10)