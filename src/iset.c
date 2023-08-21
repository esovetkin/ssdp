#include <stdlib.h>
#include "iset.h"

#ifdef RUNTEST
#include <assert.h>
int COLLISIONS = 0;
int REALLOCATE = 0;
#endif

#ifdef RUNMEMTEST
#include "random_fail_malloc.h"
#define malloc(x) random_fail_malloc(x)
#endif


static unsigned int hash(unsigned int x)
{
		// from: https://stackoverflow.com/a/12996028
		x = ((x >> 16) ^ x) * 0x45d9f3b;
		x = ((x >> 16) ^ x) * 0x45d9f3b;
		x = (x >> 16) ^ x;
		return x;
}


static int* reallocate(struct iset* self)
{
		if (self->n < self->n_max) return self->values;

		int i, n = 2*self->N;
		int *tmp = realloc(self->values, n*sizeof(*tmp));
		if (NULL == tmp) goto erealloc;

		for (i=self->N-1; i < n; ++i)
				tmp[i] = -1;

#ifdef RUNTEST
		++REALLOCATE;
		for (i=self->N-1; i < n; ++i)
				assert(-1 == tmp[i%n]);
#endif

		self->N = n;
		self->n_max = 3*n/4;

		return tmp;
erealloc:
		return NULL;
}


struct iset* iset_init(int n)
{
		struct iset* self = malloc(sizeof(*self));
		if (NULL==self) goto eself;

		if (NULL==(self->values=malloc(n*sizeof(*self->values)))) goto evalues;

		int i;
		for (i=0; i < n; ++i)
				self->values[i] = -1;

		self->N = n;
		self->n=0;
		self->n_max = 3*n/4;

		return self;
evalues:
		free(self);
eself:
		return NULL;
}


void iset_free(struct iset* self)
{
		free(self->values);
		free(self);
}


int iset_insert(struct iset* self, unsigned int key)
{
		if (NULL==(self->values=reallocate(self))) goto erealloc;

		unsigned int i = hash(key);

		while (-1 != self->values[i % self->N]) {
#ifdef RUNTEST
				++COLLISIONS;
#endif
				++i;
		}
		self->values[i % self->N] = (int) key;
		++self->n;

		return 0;
erealloc:
		return -1;
}


int iset_isin(struct iset* self, unsigned int key)
{
		unsigned int i=hash(key) % self->N;
		unsigned int i0 = i;

		while (-1 != self->values[i]) {
				if ((int)key == self->values[i]) return 1;
				if (i0 == (i=(i+1)%self->N)) break;
		}

		return 0;
}


#ifdef RUNTEST

#include <stdio.h>
#include <assert.h>

void test(int n, int k)
{
		COLLISIONS=0;
		REALLOCATE=0;

		int i;
		struct iset* set;
		assert((set=iset_init(n)));

		for (i=0; i < set->N; ++i)
				assert(-1 == set->values[i % set->N]);

		for (i=0; i < k; ++i) {
				assert(0==iset_insert(set, (unsigned int) i));
				assert(iset_isin(set, (unsigned int) i));
		}

		assert(k == set->n);

		for (i=k; i < 5*k; ++i)
				assert(0 == iset_isin(set, (unsigned int) i));

		printf("\nn = %d, k = %d\n", n, k);
		printf("        COLLISIONS = %d\n", COLLISIONS);
		printf("        REALLOCATE = %d\n", REALLOCATE);

		iset_free(set);
}


int main(void)
{
		printf("testing iset ...");

		test(1000, 1000);
		test(3000, 1000);

		printf("PASSED\n");
		return 0;
}

#endif
