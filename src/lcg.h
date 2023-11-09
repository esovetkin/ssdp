#ifndef _RUNIF_H
#define _RUNIF_H

#include <stdint.h>

/** linear conguent generator
 *
 *
 */
struct lcg {
		uint32_t* iu;
};

struct lcg* lcg_init(uint32_t seed);
void lcg_free(struct lcg* self);
double lcg_unif(struct lcg* self, int thread);

#endif
