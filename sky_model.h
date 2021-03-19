#ifndef _SKY_MODEL_H
#define _SKY_MODEL_H
void UniformSky(sky_grid *sky, sky_pos sun, double GHI, double DHI);
void PerezSky(sky_grid * sky, sky_pos sun, double GHI, double DHI, double dayofyear);
#endif
