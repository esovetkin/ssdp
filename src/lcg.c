#include "lcg.h"
#include <omp.h>
#include <stdlib.h>


struct lcg* lcg_init(uint32_t seed)
{
		struct lcg* self;
		if (NULL==(self=malloc(sizeof(*self)))) goto eself;

		if (NULL==(self->iu=malloc(omp_get_max_threads()*sizeof(*self->iu)))) goto eiu;

		for (int i=0; i < omp_get_max_threads(); ++i) {
				seed *= 2;
				self->iu[i] = ( ( seed%2 ) ? seed : seed + 1 );
		}

		return self;
eiu:
		free(self);
eself:
		return NULL;
}


void lcg_free(struct lcg* self)
{
		free(self->iu);
		self->iu = NULL;
		free(self);
}


double lcg_unif(struct lcg* self, int thread_id)
{
		self->iu[thread_id] *= 663608941l;
		return (self->iu[thread_id]*0.232830643654e-9);
}


#ifdef RUNTEST

#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>


double test(int n)
{
		int i;
		double x, tic;

		struct lcg* lcg;
		assert((lcg=lcg_init(1234)));
		
		tic = (double)clock();
		for (i=0; i < n; ++i)
				x += lcg_unif(lcg, 0);
		tic = (double)(clock()-tic)/CLOCKS_PER_SEC;
		printf("runif, n=%d, %2.0f MB/s\n", n, n*sizeof(x)*1e-6/tic);

		tic = (double)clock();
		for (i=0; i < n; ++i)
				x += rand();
		tic = (double)(clock()-tic)/CLOCKS_PER_SEC;
		printf("rand, n=%d, %2.0f MB/s\n", n, n*sizeof(x)*1e-6/tic);

		for (i=0; i < 10; ++i)
				printf("%f\n",lcg_unif(lcg,0));

		lcg_free(lcg);
		return x;
}


int main(void)
{
		printf("testing runif ...");

		printf("\n");
		test(50000000);

		printf("PASSED\n");
		return 0;
}

#endif
