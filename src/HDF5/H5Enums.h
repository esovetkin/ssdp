#pragma once
typedef enum {
    X,
    W,
    A,
    R,
} IOMode;

// Maybe use H5E_values ?
typedef enum {
    SUCCESS = 0,
    OUTOFMEMORY = 1,
    FAILURE = -1,
    // add others with different negative values
} ErrorCode;

const char* ErrorCode_to_string(ErrorCode err);