#ifndef _IO_H
#define _IO_H
topology LoadTopo(char *fn);
void WriteTopo(char *fn, topology T);
void WriteTriangles(char *fn, topology T);
void WriteDome3D(char *fn, sky_grid sky, int sun, int mask);
void WriteDome4D(char *fn, sky_grid sky, int sun, int mask);
void RasterPOA(char *fn, sky_grid sky, topology T, double albedo, double dz, double a, double tilt, double x1, double y1, double x2, double y2, int Nx, int Ny);
#endif
