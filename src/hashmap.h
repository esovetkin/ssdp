#ifndef _HASHMAP_H
#define _HASHMAP_H


#include <stdint.h>


struct hash_v {
		uint64_t a,b;
};


struct hashmap {
		struct hash_v* keys;
		void **values;
		int N;
		int cap;

		// number of entries that triggers reallocate
		int n_realloc;
};


/** init hashmap
 *
 */
struct hashmap* hashmap_init(int cap);

typedef void value_free(void *value);
void hashmap_free(struct hashmap* self, value_free* vfree);


/** inserts a (key,value) pair
 *
 *  return a non-zero value if memory failed
 */
int hashmap_insert(struct hashmap* self, void *k, int len, void *data);


/** return 0 if key not in set, 1 if key in set
 */
int hashmap_isin(struct hashmap* self, void *k, int len);
void* hashmap_get(struct hashmap* self, void *k, int len);


#endif
