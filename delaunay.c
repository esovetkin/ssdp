
#include "shull.h"
#include "ll.h"
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "delaunay.h"


SHULL_STATE SSTATE;

sh_point * MakeShullPoints(double *x, double *y, int N)
{
	int i;
	sh_point *p;
	p=malloc(N*sizeof(sh_point));
	for (i=0;i<N;i++)
	{
		p[i].x=x[i];
		p[i].y=y[i];
		p[i].index=i;
	}
	return p;
}

triangles * CollectTriangles(ll_node *n, int N, int *Nt)
{
	triangles *T;
	ll_node *o;
	sh_triangle *t;
	int i=0, Na;
	
	if (n == NULL)
		return NULL;
		
	o = n;	
	Na=2*N;
	T=malloc(Na*sizeof(triangles));
	do
	{
		t=(sh_triangle *)n->data;
		if (i==Na)
		{
			SSTATE=SHULL_ERR;
			free(T);
			return NULL;			
		}
		T[i].i=t->p[0]->index;
		T[i].j=t->p[1]->index;
		T[i].k=t->p[2]->index;
		T[i].ccx=t->cc.x;
		T[i].ccy=t->cc.y;
		T[i].ccr=sqrt(t->ccr2);
		i++;
		n = n->next;
	} while (n != o && n != NULL);
	(*Nt)=i;
	T=realloc(T,i*sizeof(triangles));
	return T;
}

triangles * Triangulate(double *x, double *y, int N, int *Nt)
{
	sh_point *P;
	sh_triangulation_data td;
	triangles *T;
	int result;
	P=MakeShullPoints(x, y, N);	
	result = delaunay(&td, P, N);
	if (result > 0)
	{
		SSTATE=SHULL_OK;
		T=CollectTriangles(td.triangles, N, Nt);

		ll_mapdestroy(td.triangles, free);
		ll_mapdestroy(td.hull_edges, free);
		ll_mapdestroy(td.internal_edges, free);
	}
	else
		SSTATE=SHULL_ERR;
	
	free(P);
	return T;
}

double dist2(double cx, double cy, double x, double y)
{
	return ((cx-x)*(cx-x)+(cy-y)*(cy-y));
}
int LineYCrossLine(double x0, double y0, double x1, double y1, double x, double y)
{
	double t, xx;
	// solve yt=y;
	// yt=y0*(1-t)+y1*t
	/*
	 * y0*(1-t)+y1*t=y
	 * y0+(y1-y0)*t=y
	 * t=(y-y0)/(y1-y0)
	 * where t between 0 and 1
	 * xx=x0*(1-t)+x1*t
	 * xx>x
	 */
	if (((y0<y)&&(y1>=y))||((y0>=y)&&(y1<y)))
	{
		if ((y==y0)&&(y==y1))
		{
			if ((x0>x)||(x1>x))
				return 1;
			return 0;
		}
		else
		{
			t=(y-y0)/(y1-y0);
			xx=x0*(1-t)+x1*t;
			if ((t>0)&&(t<1)&&(xx>=x))
				return 1;
		}
	}
	return 0;
}
int IsInTriangle(triangles t, double *X, double *Y, double x, double y, int *nearest, double *dmin)
{
	int res=0;
	double d;
	
	(*nearest)=0;
	if (LineYCrossLine(X[t.i], Y[t.i], X[t.j], Y[t.j], x, y))
		res=!res;
	if (LineYCrossLine(X[t.j], Y[t.j], X[t.k], Y[t.k], x, y))
		res=!res;
	if (LineYCrossLine(X[t.k], Y[t.k], X[t.i], Y[t.i], x, y))
		res=!res;
		
	d=dist2(t.ccx, t.ccy, x, y);
	if (d<(*dmin))
	{
		(*nearest)=1;
		(*dmin)=d;
	}
	return res;
}
	

int tsearch(triangles *T, double *X, double *Y, int N, double x, double y)
{ // find nearest by plain exhaustive search
	int i, im=0, imatch=0, n;
	double dm;
	if (N==0)
		return 0;
	dm=dist2(T[0].ccx, T[0].ccy, x, y);
	if (IsInTriangle(T[0], X, Y, x, y, &n, &dm))
		return 0;
	i=1;
	while ((imatch==0)&&(i<N))
	{
		if(IsInTriangle(T[i], X, Y, x, y, &n, &dm))
			imatch=i;
		if (n)
			im=i;
		i++;
	}
	if (imatch)
		return imatch;
	return im;	
}

#define TINY 1e-12
#define LEAFBLOCK 16
void NodeBB(triangles T,  double *X, double *Y, double bb[])
{
	bb[0]=X[T.i];
	bb[2]=X[T.i];
	bb[1]=Y[T.i];
	bb[3]=Y[T.i];
	if (X[T.j]<bb[0])
		bb[0]=X[T.j];
	if (X[T.j]>bb[2])
		bb[2]=X[T.j];
	if (Y[T.j]<bb[1])
		bb[1]=Y[T.j];
	if (Y[T.j]>bb[3])
		bb[3]=Y[T.j];
		
	if (X[T.k]<bb[0])
		bb[0]=X[T.k];
	if (X[T.k]>bb[2])
		bb[2]=X[T.k];
	if (Y[T.k]<bb[1])
		bb[1]=Y[T.k];
	if (Y[T.k]>bb[3])
		bb[3]=Y[T.k];
}

int Contained(double bb_parent[], double bb_child[])
{
	if ((bb_parent[0]-bb_child[0]<-TINY)&&(bb_parent[2]-bb_child[2]>TINY)&&(bb_parent[1]-bb_child[1]<-TINY)&&(bb_parent[3]-bb_child[3]>TINY))
		return 1;
	return 0;
}

nodetree *InitNode(double bb[])
{
	nodetree *T;
	T=malloc(sizeof(nodetree));
	T->bb[0]=bb[0];
	T->bb[1]=bb[1];
	T->bb[2]=bb[2];
	T->bb[3]=bb[3];
	T->leaves=NULL;
	T->N1=NULL;
	T->N2=NULL;
	T->N3=NULL;
	T->N4=NULL;
	return T;
}	

void AddNewNode(nodetree *P, int index, double bb[])
{
	double bbnew[4];
	/* first quarter */
	bbnew[0]=P->bb[0];
	bbnew[2]=(P->bb[0]+P->bb[2])/2;
	bbnew[1]=P->bb[1];
	bbnew[3]=(P->bb[1]+P->bb[3])/2;
	if (Contained(bbnew,bb))
	{
		P->N1=InitNode(bbnew);
		AddNewNode(P->N1, index, bb);
		return;
	}
	/* second quarter */
	bbnew[0]=(P->bb[0]+P->bb[2])/2;
	bbnew[2]=P->bb[2];
	if (Contained(bbnew,bb))
	{
		P->N2=InitNode(bbnew);
		AddNewNode(P->N2, index, bb);
		return;
	}
	/* third quarter */
	bbnew[1]=(P->bb[1]+P->bb[3])/2;
	bbnew[3]=P->bb[3];
	if (Contained(bbnew,bb))
	{
		P->N3=InitNode(bbnew);
		AddNewNode(P->N3, index, bb);
		return;
	}
	/* fourth quarter */
	bbnew[0]=P->bb[0];
	bbnew[2]=(P->bb[0]+P->bb[2])/2;
	if (Contained(bbnew,bb))
	{
		P->N4=InitNode(bbnew);
		AddNewNode(P->N4, index, bb);
		return;
	}
	/* still here? --> Add a leaf to parent and stop recursing
	printf("%e %e %e %e\n", P->bb[0],P->bb[1], bb[0],bb[1]);
	printf("%e %e %e %e\n", P->bb[2],P->bb[1], bb[2],bb[1]);
	printf("%e %e %e %e\n", P->bb[2],P->bb[3], bb[2],bb[3]);
	printf("%e %e %e %e\n", P->bb[0],P->bb[3], bb[0],bb[3]);
	printf("%e %e %e %e\n\n", P->bb[0],P->bb[1], bb[0],bb[1]);*/
	if (P->leaves)
	{
		if (P->leaves[0]%LEAFBLOCK==0)
			P->leaves=realloc(P->leaves, (P->leaves[0]+LEAFBLOCK+1)*sizeof(int));
		P->leaves[0]++;
		P->leaves[P->leaves[0]]=index;
	}
	else
	{
		P->leaves=malloc((LEAFBLOCK+1)*sizeof(int));
		P->leaves[1]=index;
		P->leaves[0]=1;
	}
	
	return;
}

void AddNode(nodetree *P, int index, double bb[])
{
	double bbnew[4];
	/* first quarter */
	bbnew[0]=P->bb[0];
	bbnew[2]=(P->bb[0]+P->bb[2])/2;
	bbnew[1]=P->bb[1];
	bbnew[3]=(P->bb[1]+P->bb[3])/2;
	if (Contained(bbnew,bb))
	{
		if (P->N1)
			AddNode(P->N1, index, bb); // data structures are already initialized
		else
		{
			P->N1=InitNode(bbnew); // first element in this sector
			AddNewNode(P->N1, index, bb);
		}		
		return;
	}
	
	/* second quarter */
	bbnew[0]=(P->bb[0]+P->bb[2])/2;
	bbnew[2]=P->bb[2];
	if (Contained(bbnew,bb))
	{
		if (P->N2)
			AddNode(P->N2, index, bb); // data structures are already initialized
		else
		{
			P->N2=InitNode(bbnew); // first element in this sector
			AddNewNode(P->N2, index, bb);
		}		
		return;
	}
	/* third quarter */
	bbnew[1]=(P->bb[1]+P->bb[3])/2;
	bbnew[3]=P->bb[3];
	if (Contained(bbnew,bb))
	{
		if (P->N3)
			AddNode(P->N3, index, bb); // data structures are already initialized
		else
		{
			P->N3=InitNode(bbnew); // first element in this sector
			AddNewNode(P->N3, index, bb);
		}		
		return;
	}
	/* fourth quarter */
	bbnew[0]=P->bb[0];
	bbnew[2]=(P->bb[0]+P->bb[2])/2;
	if (Contained(bbnew,bb))
	{
		if (P->N4)
			AddNode(P->N4, index, bb); // data structures are already initialized
		else
		{
			P->N4=InitNode(bbnew); // first element in this sector
			AddNewNode(P->N4, index, bb);
		}		
		return;
	}
	/* still here? --> Add a leaf to parent and stop recursing*/
	if (P->leaves)
	{
		if (P->leaves[0]%LEAFBLOCK==0)
			P->leaves=realloc(P->leaves, (P->leaves[0]+LEAFBLOCK+1)*sizeof(int));
		P->leaves[0]++;
		P->leaves[P->leaves[0]]=index;
	}
	else
	{
		P->leaves=malloc((LEAFBLOCK+1)*sizeof(int));
		P->leaves[1]=index;
		P->leaves[0]=1;
	}
	return;
}
nodetree *InitTree(triangles *T, int Nt, double *X, double *Y, int Np)
{
	nodetree *P;
	int i;
	double bb[4];
		
	if (Nt==0)
		return NULL;
	bb[0]=X[0];
	bb[2]=X[0];
	bb[1]=Y[0];
	bb[3]=Y[0];
	
	for (i=0;i<Np;i++)
	{
		if (X[i]<bb[0])
			bb[0]=X[i];
		if (X[i]>bb[2])
			bb[2]=X[i];
		if (Y[i]<bb[1])
			bb[1]=Y[i];
		if (Y[i]>bb[3])
			bb[3]=Y[i];
	}
	
	P=InitNode(bb);
	
	for (i=0;i<Nt;i++)
	{
		NodeBB(T[i], X, Y, bb);
		AddNode(P, i, bb);
	}
	return P;
}

void FreeTree(nodetree *P)
{
	if (P->leaves)
		free(P->leaves);
	if (P->N1)
	{
		 FreeTree(P->N1);
		 free(P->N1);
		 P->N1=NULL;
	}
	if (P->N2)
	{
		 FreeTree(P->N2);
		 free(P->N2);
		 P->N2=NULL;
	}
	if (P->N3)
	{
		 FreeTree(P->N3);
		 free(P->N3);
		 P->N3=NULL;
	}
	if (P->N4)
	{
		 FreeTree(P->N4);
		 free(P->N4);
		 P->N4=NULL;
	}
}
int PContained(double bb[], double x, double y)
{
	if ((bb[0]-x<-TINY)&&(bb[2]-x>TINY)&&(bb[1]-y<-TINY)&&(bb[3]-y>TINY))
		return 1;
	return 0;
}

int IsInTriangle2(triangles t, double *X, double *Y, double x, double y)
{
	int res=0;
	if ((X[t.i]<x)&&(X[t.j]<x)&&(X[t.k]<x))
		return 0;
	if ((X[t.i]>x)&&(X[t.j]>x)&&(X[t.k]>x))
		return 0;
	if ((Y[t.i]<y)&&(Y[t.j]<y)&&(Y[t.k]<y))
		return 0;
	if ((Y[t.i]>y)&&(Y[t.j]>y)&&(Y[t.k]>y))
		return 0;
	
	
	if (LineYCrossLine(X[t.i], Y[t.i], X[t.j], Y[t.j], x, y))
		res=!res;
	if (LineYCrossLine(X[t.j], Y[t.j], X[t.k], Y[t.k], x, y))
		res=!res;
	if (LineYCrossLine(X[t.k], Y[t.k], X[t.i], Y[t.i], x, y))
		res=!res;
		
	return res;
}
nodetree *SearchLeaf(nodetree *P, double x, double y, int restart)
{
	if ((P->leaves)&&(!restart))
		return P; // return to process leaves
	if (P->N1)
		if (PContained(P->N1->bb,x,y))
			return SearchLeaf(P->N1, x, y, 0);
	if (P->N2)
		if (PContained(P->N2->bb,x,y))
			return SearchLeaf(P->N2, x, y, 0);
	if (P->N3)
		if (PContained(P->N3->bb,x,y))
			return SearchLeaf(P->N3, x, y, 0);
	if (P->N4)
		if (PContained(P->N4->bb,x,y))
			return SearchLeaf(P->N4, x, y, 0);
	return NULL;
}	

int treesearch(triangles *T, nodetree *P, double *X, double *Y, int N, double x, double y)
{
	int iter=0;
	int i;
	if (!P)
		return tsearch(T, X, Y, N, x, y);
		
	while (P)
	{
		if (P->leaves)
		{
			for (i=1;i<=P->leaves[0];i++)
			{
				iter++;
				if(IsInTriangle2(T[P->leaves[i]], X, Y, x, y))
				{
					//printf("%d %d %e\n", iter, N, (double)iter/(double)N);
					return P->leaves[i];
				}
			}
		}
		P=SearchLeaf(P, x, y, 1);
	}
	printf("Failed! %e %e\n", x, y);	
	// failed!
	return tsearch(T, X, Y, N, x, y);
}
