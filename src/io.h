#ifndef _IO_H
#define _IO_H
topology LoadTopo(char *fn);
void WriteTopo(char *fn, topology *T);
void WriteTriangles(char *fn, topology *T);
void WriteDome4D(char *fn, sky_grid *sky, location *l);
double **ReadArrays(char *fn, int Narr, int *N);
void WriteArrays(char *fn, double **data, int Narr, int N);
#endif
