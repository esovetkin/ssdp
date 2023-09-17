#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "minfill.h"


struct ipair {
		int i;
		int s;
};


struct queue {
		struct ipair *x;
		int left;
		int right;
		int N;
};


struct queue* queue_init(int N)
{
		struct queue *self;
		if (NULL==(self=malloc(sizeof(*self)))) goto eself;
		self->N = N > 0 ? N : 1;
		self->left = 0;
		self->right = -1;
		if (NULL==(self->x=malloc(self->N*sizeof(*(self->x))))) goto ex;

		return self;
		free(self->x);
ex:
		free(self);
eself:
		return NULL;
}


void queue_free(struct queue *self)
{
		if (NULL==self) return;
		free(self->x);
		free(self);
}


int queue_isempty(struct queue *self)
{
		return self->left > self->right;
}


static struct ipair* reallocate(struct queue* self)
{
		if (queue_isempty(self)) {
				self->left = 0;
				self->right = -1;
		}

		if (self->right < self->N - 1)
				return self->x;

		int n = 2*self->N;
		struct ipair *tmp = realloc(self->x, n*sizeof(*tmp));
		if (NULL==tmp) goto erealloc;
		self->N = n;
		return tmp;
erealloc:
		return NULL;
}


int queue_append(struct queue *self, struct ipair item)
{
		if (NULL==(self->x=reallocate(self))) goto erealloc;

		self->x[++self->right] = item;
		return 0;
erealloc:
		return -1;
}


struct ipair queue_popleft(struct queue *self)
{
		return self->x[self->left++];
}


struct minfill* minfill_init(double *z, int nx, int ny, double missing_value)
{
		struct minfill *self;
		if (NULL==(self=malloc(sizeof(*self)))) goto eself;
		self->z = z;
		self->nx = nx;
		self->ny = ny;
		self->missing_value = nextafter(missing_value, DBL_MAX);
		self->maxwalk = -1;
		self->direction = -1;

		if (NULL==(self->visited=iset_init(nx*ny / 2))) goto evisited;
		if (NULL==(self->queue=queue_init(nx*ny / 10))) goto equeue;
		if (NULL==(self->component=queue_init(nx*ny / 10))) goto ecomponent;


		return self;
		queue_free(self->component);
ecomponent:
		queue_free(self->queue);
equeue:
		iset_free(self->visited);
evisited:
		free(self);
eself:
		return NULL;
}


void minfill_free(struct minfill* self)
{
		if (NULL==self) return;
		iset_free(self->visited);
		queue_free(self->queue);
		queue_free(self->component);
		free(self);
}


static void neighbours(int x, int *a, int nx, int ny)
{
		int i = 0, row = x % nx, col = x / nx;

		if (col - 1 > 0)
				a[i++] = (col-1)*nx + row;
		if (col + 1 < ny)
				a[i++] = (col+1)*nx + row;

		if (row - 1 > 0) {
				a[i++] = col*nx + (row-1);
				if (col - 1 > 0)
						a[i++] = (col-1)*nx + (row-1);
				if (col + 1 < ny)
						a[i++] = (col+1)*nx + (row-1);
		}
		if (row + 1 < nx) {
				a[i++] = col * nx + (row + 1);
				if (col - 1 > 0)
						a[i++] = (col-1)*nx + (row+1);
				if (col + 1 < ny)
						a[i++] = (col+1)*nx + (row+1);
		}
		a[i] = -1;
}


static int cmp(double x, double y, int d)
{
		if (d > 0)
				return x > y;

		return x <= y;
}


static void update_v(struct minfill* self, double* v, double x)
{
		if (*v < self->missing_value) {
				*v = x;
				return;
		}

		if (cmp(*v, x, self->direction)) return;

		*v = x;
}


static int walk(struct minfill* self, int q)
{
		int i;
		struct ipair x;
		double v = nextafter(self->missing_value, -DBL_MAX);
		if (queue_append(self->queue, (struct ipair){q, 0})) goto eappend;

		while (!queue_isempty(self->queue)) {
				x = queue_popleft(self->queue);

				if (iset_isin(self->visited, x.i)) continue;
				iset_insert(self->visited, x.i);

				if (queue_append(self->component, (struct ipair){x.i, 0}))
						goto eappend;
				if (self->maxwalk > 0 &&
					x.s >= self->maxwalk &&
					v >= self->missing_value)
						continue;

				neighbours(x.i, self->n, self->nx, self->ny);
				for (i=0; i < 9; ++i) {
						if (self->n[i] < 0) break;
						if (iset_isin(self->visited, self->n[i])) continue;

						if (self->z[self->n[i]] >= self->missing_value) {
								update_v(self, &v, self->z[self->n[i]]);
								continue;
						}

						if (queue_append(self->queue,
										 (struct ipair){self->n[i], x.s + 1}))
								goto eappend;
				}
		}

		while(!queue_isempty(self->component)) {
				x = queue_popleft(self->component);

				if (x.i < 0 || x.i >= self->nx*self->ny) goto ez;
				self->z[x.i] = v;
		}

		return 0;
ez:
		printf("Error: bug in minfill\n");
eappend:
		return -1;
}


int minfill_fill(struct minfill* self)
{
		int i;
		for (i=0; i < self->nx*self->ny; ++i) {
				if (self->z[i] >= self->missing_value)
						continue;

				if (walk(self, i)) goto efill;
		}

		return 0;
efill:
		return -1;
}


#ifdef RUNTEST

#include <stdio.h>
#include <math.h>
#include <assert.h>


void test_queue(int N)
{
		int i;
		struct queue *x;
		x = queue_init(N);

		for (i=0; i < 5*N; ++i)
				assert(0 == queue_append(x, (struct ipair){i,0}));

		for (i=0; i < 5*N; ++i)
				assert(i == queue_popleft(x).i);

		assert(queue_isempty(x));
		queue_free(x);
}


void print_z(double *z, int n, int m)
{
		int i;

		for (i=0; i < n*m; ++i) {
				if (0 == i % n)
						printf("\n");
				printf("%2.0f ",z[i]);
		}
		printf("\n");
}


double* init_z(int n, int m)
{
		int i;
		double *z = malloc(n*m*sizeof(*z));
		assert(z);

		for (i=0; i < n*m; ++i) {
				if (0 ==i%7 || 1==i%7 || 2==i%7) {
						z[i] = -1;
						continue;
				}
				z[i] = n*m - (double)i;
		}

		return z;
}


void test_minfill(int n, int m)
{
		int i;
		double *z = init_z(n, m);
		// print_z(z,n,m);

		struct minfill *mf = minfill_init(z, n, m, (double)0);
		assert(mf);
		mf->maxwalk = -3;
		mf->direction = -1;

		minfill_fill(mf);
		minfill_free(mf);
		// print_z(z,n,m);

		for (i=0; i < n*m; ++i)
				assert(z[i] > 0);

		free(z);
}


int main(void)
{
		printf("testing minfill ...");

		test_queue(0);
		test_queue(1);
		test_queue(2);
		test_queue(3);
		test_queue(4);
		test_queue(400);

		test_minfill(13, 17);
		test_minfill(10, 10);
		test_minfill(20, 20);
		test_minfill(1000, 1000);

		printf("PASSED\n");
}

#endif
