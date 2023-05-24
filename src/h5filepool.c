#include "h5filepool.h"
/* TODO
 * For some reason the pool gets initialized with a poolsize of zero
 */
struct H5FileIOHandlerPool *g_h5filepool = NULL;;
void init_h5filepool(){
    g_h5filepool = H5FileIOHandlerPool_init();
}

void free_h5filepool(){
    H5FileIOHandlerPool_free(&g_h5filepool);
}

