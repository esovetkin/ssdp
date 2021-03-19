#ifndef _IO_H
#define _IO_H
topology LoadTopo(char *fn);

void WriteDome3D(char *fn, sky_grid sky, int sun, int mask);
void WriteDome4D(char *fn, sky_grid sky, int sun, int mask);
#endif
