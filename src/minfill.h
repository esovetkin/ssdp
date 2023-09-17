#ifndef _MINFILL_H
#define _MINFILL_H

#include "iset.h"

struct queue;

/** fill missing values based on maximum or minimum value available
   within the given connected components

   @z, @nx, @ny: **borrowed** 2d array

   @missing_value: upper limit, values below (or equal) are considered
   missing.

   @maxwalk: maximum number of walking step for search for the
   non-missing value (default: -1, i.e. check the whole connected
   component)

   @direction: either 1 (search maximum non-missing value), or -1
   (fill minimum non-missing value)

   @visited: hashmap to keep track of visited pixels

   @queue: queue of nodes to visit

   @component: current connected component
*/
struct minfill {
		double *z;
		int nx, ny;
		double missing_value;
		int maxwalk;
		int direction;
		struct iset *visited;
		struct queue *queue;
		struct queue *component;
		int n[9];
};

struct minfill* minfill_init(double *z, int nx, int ny, double missing_value);
void minfill_free(struct minfill *self);
int minfill_fill(struct minfill *self);

#endif
