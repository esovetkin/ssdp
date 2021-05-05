#include "shull.h"
#include "ll.h"
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "delaunay.h"
#include "print.h"
#include "error.h"

sh_point * MakeShullPoints(double *x, double *y, int N)
{
	int i;
	sh_point *p;
	
	// ERRORFLAG MALLOCFAILDELAUNAY  "Error memory allocation failed for Delaunay points"
	if ((p=malloc(N*sizeof(sh_point)))==NULL)
	{
		AddErr(MALLOCFAILDELAUNAY);
		return p;
	}
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
	
	// ERRORFLAG MALLOCFAILCOLLTRIANGLE  "Error memory allocation failed in collecting triangles"
	if ((T=malloc(Na*sizeof(triangles)))==NULL)
	{
		AddErr(MALLOCFAILCOLLTRIANGLE);
		return T;
	}
	do
	{
		t=(sh_triangle *)n->data;
		if (i==Na)
		{
			// ERRORFLAG TOOMANYTRIANGLES  "Error found more triangles than there should be, panic!"
			AddErr(TOOMANYTRIANGLES);
			free(T);
			return NULL;			
		}
		T[i].i=t->p[0]->index; // original indexes of the points
		T[i].j=t->p[1]->index;
		T[i].k=t->p[2]->index;
		T[i].ccx=(t->p[0]->x+t->p[1]->x+t->p[2]->x)/3;// average coordinate of the triangle
		T[i].ccy=(t->p[0]->y+t->p[1]->y+t->p[2]->y)/3;
		i++;
		n = n->next;
	} while (n != o && n != NULL);
	(*Nt)=i;
	if ((T=realloc(T,i*sizeof(triangles)))==NULL)
	{
		AddErr(MALLOCFAILCOLLTRIANGLE);
		return T;
	}
	return T;
}

#define DISTEPS 1e-12
void Distort(double *x, double *y, int N)
{
	// in case triangulation fails we can try distorting the arrays a bit
	int i;
	srand(time(NULL));
	for (i=0;i<N;i++)
	{
		x[i]*=(1.0+DISTEPS*(((double)rand())/RAND_MAX-0.5));
		y[i]*=(1.0+DISTEPS*(((double)rand())/RAND_MAX-0.5));
	}
}

triangles * Triangulate(double *x, double *y, int N, int *Nt)
{
	sh_point *P;
	sh_triangulation_data td;
	triangles *T;
	int result;
	P=MakeShullPoints(x, y, N);	
	if (ssdp_error_state)
		return NULL;
	result = delaunay(&td, P, N);
	if (result > 0)
	{
		T=CollectTriangles(td.triangles, N, Nt);

		ll_mapdestroy(td.triangles, free);
		ll_mapdestroy(td.hull_edges, free);
		ll_mapdestroy(td.internal_edges, free);
	}
	else
	{
		// give it one more try after a random distortion of the mesh
		Print(WARNING,"Warning triangulation failed, retry with random pertubations\n");
		free(P);
		Distort(x,y, N);
		P=MakeShullPoints(x, y, N);	
		result = delaunay(&td, P, N);
		if (result > 0)
		{
			T=CollectTriangles(td.triangles, N, Nt);
	
			ll_mapdestroy(td.triangles, free);
			ll_mapdestroy(td.hull_edges, free);
			ll_mapdestroy(td.internal_edges, free);
		}
		else
		{
			// ERRORFLAG TRIANGULATIONFAILED  "Error triangulation failed!"
			AddErr(TRIANGULATIONFAILED);
			T=NULL;
		}
	}
	free(P);
	return T;
}

// measure quadratic distance
double dist2(double cx, double cy, double x, double y)
{
	return ((cx-x)*(cx-x)+(cy-y)*(cy-y));
}

// check whether a horizontal line from (x,y)->(inf,y)
// passes through a line segment (x0,y0)->(x1,y1)
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

/* here follows some very crude and simple quaternary search tree to locate triangles
 * The concept is simple:
 * - build a tree of bounding boxes
 * - each bounding box may be divided in 4, equal sized, sub-bounding boxes
 * - each branch may have leafs. Leafs are those triangles that do fit in the 
 *   current bounding box but in none of the sub bounding boxes
 * With this we can search for a triangle given a coordinate as we just traverse 
 * down the tree and check every leaf on the way
 * The downside is that there are many leafs at the bottom of the tree as any triangle
 * that is cut by the any of the sub-bounding boxes will become a leaf. In practice I see
 * that I need to check a few percent of the triangles for any search (I saw 4%). I suppose
 * that should be fairly size independent, or? */

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
	// ERRORFLAG MALLOCFAILSEARTREENODE  "Error memory allocation failed in creating search tree branches"
	if ((T=malloc(sizeof(nodetree)))==NULL)
	{
		AddErr(MALLOCFAILSEARTREENODE);
		return T;
	}
	T->bb[0]=bb[0];
	T->bb[1]=bb[1];
	T->bb[2]=bb[2];
	T->bb[3]=bb[3];
	T->leafs=NULL;
	T->N1=NULL;
	T->N2=NULL;
	T->N3=NULL;
	T->N4=NULL;
	return T;
}	

/* when initializing go up the tree creating new branches
 * untill we can hang our triangle to the highest branch
 * that still fits the triangle */
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
	/* still here? --> Add a leaf to parent and stop recursing */
	// ERRORFLAG MALLOCFAILLEAFS  "Error memory allocation failed in creating search tree leafs"
	if (P->leafs)
	{
		if (P->leafs[0]%LEAFBLOCK==0)
			P->leafs=realloc(P->leafs, (P->leafs[0]+LEAFBLOCK+1)*sizeof(int));
		if (P->leafs==NULL)
		{
			AddErr(MALLOCFAILLEAFS);
			return;
		}
		
		P->leafs[0]++;
		P->leafs[P->leafs[0]]=index;
	}
	else
	{
		P->leafs=malloc((LEAFBLOCK+1)*sizeof(int));
		if (P->leafs==NULL)
		{
			AddErr(MALLOCFAILLEAFS);
			return;
		}
		P->leafs[1]=index;
		P->leafs[0]=1;
	}
	
	return;
}

/* when initializing go up the tree. As long as branches are there
 * that still fit the triangle we go up. If we find we are at the top
 * of the current tree and find we need to grow further we switche to 
 * the AddNewNode routine. */
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
	if (P->leafs)
	{
		if (P->leafs[0]%LEAFBLOCK==0)
			P->leafs=realloc(P->leafs, (P->leafs[0]+LEAFBLOCK+1)*sizeof(int));
		if (P->leafs==NULL)
		{
			AddErr(MALLOCFAILLEAFS);
			return;
		}
		P->leafs[0]++;
		P->leafs[P->leafs[0]]=index;
	}
	else
	{
		P->leafs=malloc((LEAFBLOCK+1)*sizeof(int));
		if (P->leafs==NULL)
		{
			AddErr(MALLOCFAILLEAFS);
			return;
		}
		P->leafs[1]=index;
		P->leafs[0]=1;
	}
	return;
}
/* grow a tree */
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

/* kill a tree */
void FreeTree(nodetree *P)
{
	if (P->leafs)
	{
		free(P->leafs);
		P->leafs=NULL;
	}
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
/* check point within a bounding box */
int PContained(double bb[], double x, double y)
{
	if ((bb[0]-x<-TINY)&&(bb[2]-x>TINY)&&(bb[1]-y<-TINY)&&(bb[3]-y>TINY))
		return 1;
	return 0;
}

/* check point within a triangle */
int IsInTriangle2(triangles t, double *X, double *Y, double x, double y)
{
	int res=0;
	// check bounding box
	if ((X[t.i]<x)&&(X[t.j]<x)&&(X[t.k]<x))
		return 0;
	if ((X[t.i]>x)&&(X[t.j]>x)&&(X[t.k]>x))
		return 0;
	if ((Y[t.i]<y)&&(Y[t.j]<y)&&(Y[t.k]<y))
		return 0;
	if ((Y[t.i]>y)&&(Y[t.j]>y)&&(Y[t.k]>y))
		return 0;
	
	
	// still here, do a thorough check
	if (LineYCrossLine(X[t.i], Y[t.i], X[t.j], Y[t.j], x, y))
		res=!res;
	if (LineYCrossLine(X[t.j], Y[t.j], X[t.k], Y[t.k], x, y))
		res=!res;
	if (LineYCrossLine(X[t.k], Y[t.k], X[t.i], Y[t.i], x, y))
		res=!res;
		
	return res;
}

/* climb up a tree till we find new leafs */
nodetree *SearchLeaf(nodetree *P, double x, double y, int restart)
{
	if ((P->leafs)&&(!restart))
		return P; // return to process leafs
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

/* find that one particular leaf */
int treesearch(triangles *T, nodetree *P, double *X, double *Y, int N, double x, double y)
{
	int i;
	int im=0;
	double dm=-1, d;
	if (!P)
		return tsearch(T, X, Y, N, x, y);
		
	dm=dist2(T[P->leafs[im]].ccx, T[P->leafs[im]].ccy, x, y);
	while (P)
	{
		if (P->leafs)
		{
			for (i=1;i<=P->leafs[0];i++)
			{
				if(IsInTriangle2(T[P->leafs[i]], X, Y, x, y))
					return P->leafs[i];
				else
				{
					d=dist2(T[P->leafs[i]].ccx, T[P->leafs[i]].ccy, x, y);
					if (dm>d)
					{
						dm=d;
						im=P->leafs[i];
					}
				}
			}
				
		}
		P=SearchLeaf(P, x, y, 1);
	}
	// failed! Hopefully the point just lies out of the convex hull
	// return the closest triangle we encountered on the way up
	return im;
}
