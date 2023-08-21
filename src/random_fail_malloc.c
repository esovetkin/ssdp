#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "random_fail_malloc.h"


static int SRAND_INIT = 0;
static double MALLOCFAIL = 0.05;


static void init()
{
		if (SRAND_INIT)
				return;

		char *e;
		int i = (int) time(NULL);

		const char *v_init=getenv("SRAND");
		const char *v_fail=getenv("MALLOCFAIL");

		if (NULL != v_init) {
				i = strtol(v_init, &e, 10);
				if (*e != '\0')
						i = (int) time(NULL);
		}

		printf("SRAND = %d\n", i);
		srand(i);
		SRAND_INIT=1;

		if (NULL != v_fail) {
				double x = strtof(v_fail, &e);
				if (*e == '\0')
						MALLOCFAIL=x;
		}
		printf("MALLOCFAIL = %f\n", MALLOCFAIL);
}


void* random_fail_malloc(size_t n)
{
		init();
		double x = rand()/((double)RAND_MAX);

		if (x < MALLOCFAIL) {
				printf("random_fail_malloc: %f < %f\n", x, MALLOCFAIL);
				return NULL;
		}

		printf("random_fail_malloc: %f >= %f\n", x, MALLOCFAIL);
		return malloc(n);
}
