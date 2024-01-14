#include <stdlib.h>
#include "hashmap.h"
#include "murmurhash3.h"

#ifdef RUNTEST
#include <assert.h>
#include <stdio.h>
int COLLISIONS = 0;
int REALLOCATE = 0;
#endif


static struct hash_v hash(const void * key, int len)
{
		uint64_t h[2];
		MurmurHash3_x64_128(key, len, 0, &h);
		if (0 == h[0] && 0 == h[1]) ++h[1];
		return (struct hash_v){.a = h[0],.b=h[1]};
}


static int istaken(struct hash_v* key)
{
		if (0 == key->a && 0 == key->b) return 0;

		return 1;
}


static int reallocate(struct hashmap* self)
{
		if (self->N < self->n_realloc) return 0;

#ifdef RUNTEST
		++REALLOCATE;
#endif

		int i, n = 2*self->cap;
		struct hash_v *tmp;
		void **values=NULL;
		if (NULL==(tmp=calloc(n,sizeof(*tmp))))
				goto etmp;
		if (self->values)
				if (NULL==(values=calloc(n,sizeof(*values))))
						goto evls;

#ifdef RUNTEST
		for (i=0; i < n; ++i)
				assert(0 == istaken(tmp+i));
#endif

		for (i=0; i < self->cap; ++i) {
				if (0 == istaken(self->keys+i)) continue;
				int k = self->keys[i].b % n;
				while (istaken(tmp+k))
						k = (k + 1) % n;
				tmp[k] = self->keys[i];
				if (self->values)
						values[k] = self->values[i];
		}

		free(self->keys); self->keys=tmp;
		if (self->values) {free(self->values); self->values=values;}
		self->cap = n;
		self->n_realloc = 3*n/4;

		return 0;
evls:
etmp:
		return -1;
}


struct hashmap* hashmap_init(int cap)
{
		struct hashmap* self = malloc(sizeof(*self));
		if (NULL==self) goto eself;

		if (NULL==(self->keys=calloc(cap,sizeof(*self->keys)))) goto ekeys;
		self->values = NULL;
		self->cap = cap;
		self->N=0;
		self->n_realloc = 3*self->cap/4;

		return self;
		free(self->keys); self->keys=NULL;
ekeys:
		free(self);
eself:
		return NULL;
}


void hashmap_free(struct hashmap* self, value_free* vfree)
{
		if (NULL == self) return;
		if (self->keys) {free(self->keys); self->keys=NULL;}
		if (self->values && vfree) {
				int i;
				for (i=0; i < self->cap; ++i) {
						if (NULL == self->values[i])
								continue;

						vfree(self->values[i]);
						self->values[i] = NULL;
				}
				free(self->values);
				self->values=NULL;
		}
		free(self);
}


int hashmap_insert(struct hashmap* self, void *k, int len, void *data)
{
		if (reallocate(self)) goto erealloc;

		struct hash_v h = hash(k, len);
		int i = h.b % self->cap;

		while (istaken(self->keys+i)) {
#ifdef RUNTEST
				++COLLISIONS;
#endif
				i = (i + 1) % self->cap;
		}

		if (data && (NULL==self->values)) {
				self->values = calloc(self->cap,sizeof(*self->values));
				if (NULL == self->values) goto erealloc;
		}

		self->keys[i] = h;
		if (data) self->values[i] = data;
		++self->N;

		return 0;
erealloc:
		return -1;
}


static int hashmap_index(struct hashmap* self, void *k, int len)
{
		struct hash_v h = hash(k, len);
		int i = h.b % self->cap;
		int i0 = i;

		while (istaken(self->keys+i)) {
				if (h.a == self->keys[i].a &&
					h.b == self->keys[i].b)
						return i;

				if (i0 == (i=(i+1)%self->cap)) break;
		}

		return -1;
}


int hashmap_isin(struct hashmap* self, void *k, int len)
{
		int i = hashmap_index(self, k, len);
		if (-1 == i) return 0;

		return 1;
}


void* hashmap_get(struct hashmap* self, void *k, int len)
{
		if (NULL == self->values) return NULL;

		int i = hashmap_index(self, k, len);
		if (-1 == i) return NULL;
		return self->values[i];
}


#ifdef RUNTEST

#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>


void test(int n, int k)
{
		COLLISIONS=0;
		REALLOCATE=0;

		int i;
		struct hashmap* set;
		assert((set=hashmap_init(n)));

		for (i=0; i < k; ++i)
				assert(0==hashmap_insert(set, &i, sizeof(i), NULL));

		for (i=0; i < k; ++i)
				assert(hashmap_isin(set, &i, sizeof(i)));

		assert(k == set->N);

		for (i=k; i < 5*k; ++i)
				assert(0 == hashmap_isin(set, &i, sizeof(i)));

		printf("\nn = %d, k = %d\n", n, k);
		printf("        COLLISIONS = %d\n", COLLISIONS);
		printf("        REALLOCATE = %d\n", REALLOCATE);

		hashmap_free(set, NULL);
}


void test_s()
{
		struct hashmap* set;
		assert((set=hashmap_init(10)));

		char s[20] = "test";
		assert(0 == hashmap_insert(set, s, strlen(s)*sizeof(*s), NULL));
		assert(hashmap_isin(set, s, 4*sizeof(*s)));
		assert(0 == hashmap_isin(set, s, 5*sizeof(*s)));

		hashmap_free(set, NULL);
}


void test_data(int N)
{
		struct hashmap* map;
		assert((map=hashmap_init(10)));
		int i=0;
		for (i=0; i < N; ++i) {
				int *x = malloc(2*sizeof(*x));
				assert(x);
				x[0] = i; x[1] = i*i;
				assert(0 == hashmap_insert(map, &i, sizeof(i), x));
		}

		for (i=0; i < N; ++i) {
				assert(hashmap_isin(map, &i, sizeof(i)));
				int *x = hashmap_get(map, &i, sizeof(i));
				assert((i == x[0]));
				assert((i*i == x[1]));
		}

		hashmap_free(map, free);
}


void test_isin_speed(int n, int k)
{
		int i;
		double tic;
		struct hashmap* set;
		assert((set=hashmap_init(n)));

		for (i=0; i < k; ++i)
				assert(0 == hashmap_insert(set, &i, sizeof(i), NULL));

		tic = (double)clock();
		for (i=0; i < n; ++i)
				hashmap_isin(set, &i, sizeof(i));
		tic = (double)(clock()-tic)/CLOCKS_PER_SEC;
		printf("isin benchmark, n=%d, k=%d, %2.0f MB/s\n",
			   n, k, sizeof(n)*n*1e-6 / tic);

		hashmap_free(set, NULL);
}


int main(void)
{
		printf("testing iset ...");

		test_s();
		test(5, 1000);
		test(1000, 1000);
		test(3000, 1000);
		test_isin_speed(70000, 2000);
		test_isin_speed(700000, 20000);
		test_data(1000);
		test_data(10000);

		printf("PASSED\n");
		return 0;
}

#endif
