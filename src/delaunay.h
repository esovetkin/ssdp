
//BEGIN_SSDP_EXPORT
typedef struct nodetree {
	int *leafs;
	double bb[4];
	struct nodetree *N1;
	struct nodetree *N2;
	struct nodetree *N3;
	struct nodetree *N4;
} nodetree;
typedef struct triangles {
	int i, j, k;
	double ccx, ccy, ccz; 
} triangles;
//END_SSDP_EXPORT
typedef enum {SHULL_OK, SHULL_ERR} SHULL_STATE;
extern SHULL_STATE SSTATE;
triangles * Triangulate(double *x, double *y, int N, int *Nt);
int tsearch(triangles *T, double *X, double *Y, int N, double x, double y);
nodetree *InitTree(triangles *T, int Nt, double *X, double *Y, int Np);
void FreeTree(nodetree *P);
int treesearch(triangles *T, nodetree *P, double *X, double *Y, int N, double x, double y);
