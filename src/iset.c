#include <stdlib.h>
#include <string.h>
#include "iset.h"

#ifdef RUNTEST
#include <assert.h>
int COLLISIONS = 0;
int REALLOCATE = 0;
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
		int *tmp = malloc(n*sizeof(*tmp));
		if (NULL == tmp) goto erealloc;
		memset(tmp, -1, n*sizeof(*tmp));

#ifdef RUNTEST
		++REALLOCATE;
		for (i=0; i < n; ++i)
				assert(-1 == tmp[i]);
#endif

		for (i=0; i < self->N; ++i) {
				if (-1 == self->values[i]) continue;

				const unsigned int key = self->values[i];
				unsigned int k = hash(key);
				while (-1 != tmp[k % n]) ++k;
				tmp[k % n] = (int) key;
		}

		free(self->values);
		self->values=tmp;
		self->N = n;
		self->n_max = 3*n/4;

		return tmp;
erealloc:
		free(self->values); self->values=NULL;
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
		if (!self)
				return;
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
#include <time.h>


void test(int n, int k)
{
		COLLISIONS=0;
		REALLOCATE=0;

		int i;
		struct iset* set;
		assert((set=iset_init(n)));

		for (i=0; i < set->N; ++i)
				assert(-1 == set->values[i % set->N]);

		for (i=0; i < k; ++i)
				assert(0==iset_insert(set, (unsigned int) i));

		for (i=0; i < k; ++i)
				assert(iset_isin(set, (unsigned int) i));

		assert(k == set->n);

		for (i=k; i < 5*k; ++i)
				assert(0 == iset_isin(set, (unsigned int) i));

		printf("\nn = %d, k = %d\n", n, k);
		printf("        COLLISIONS = %d\n", COLLISIONS);
		printf("        REALLOCATE = %d\n", REALLOCATE);

		iset_free(set);
}


void test_isin_speed(int n, int k)
{
		int i;
		double tic;
		struct iset* set;
		assert((set=iset_init(n)));

		for (i=0; i < k; ++i)
				assert(0 == iset_insert(set, (unsigned int) i));

		tic = (double)clock();
		for (i=0; i < n; ++i)
				iset_isin(set, (unsigned int) i);
		tic = (double)(clock()-tic)/CLOCKS_PER_SEC;
		printf("isin benchmark, n=%d, k=%d, %2.0f MB/s\n",
			   n, k, sizeof(n)*n*1e-6 / tic);

		iset_free(set);
}


int main(void)
{
		printf("testing iset ...");

		test(1000, 1000);
		test(3000, 1000);
		test_isin_speed(70000, 2000);
		test_isin_speed(700000, 20000);

		printf("PASSED\n");
		return 0;
}

#endif
