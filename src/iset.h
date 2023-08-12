#ifndef _ISET_H
#define _ISET_H

/** simple hashmap to store a set of non-negative integers
 *
 *  values[hash(i)] = -1 if not set, i
 *  N allocated size of values
 *  n number of currently stored items
 *  n_max maximum value of n, when values are reallocated
 *
 */
struct iset {
		int* values;
		int N;
		int n;
		int n_max;
};


struct iset* iset_init(int size);
void iset_free(struct iset* self);

/** inserts a key
 *
 *  return a non-zero value if memory failed
 */
int iset_insert(struct iset* self, unsigned int key);


/** return 0 if key not in set, 1 if key in set
 */
int iset_isin(struct iset* self, unsigned int key);

#endif
