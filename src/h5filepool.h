#include "HDF5/hd5extention/src/H5FileIO.h"
#include <stdlib.h>

extern struct H5FileIOHandlerPool* g_h5filepool;

void init_h5filepool(void);

void free_h5filepool(void);
