#ifndef _PNGOUT_H
#define _PNGOUT_H


/**
 * NORM_NONE just casts double as (uint8_t)
 *
 * NORM_MAXMIN scales array 255*(z - z.min()) / (z.max() - z.min());
 */
enum normalisation {
        NORM_NONE,
        NORM_MAXMIN
};


/**
 * normalise and write double array to a file
 *
 * @ofn: output filename
 *
 * @z, @nx, @ny: array of size nx*ny
 *
 * @norm: how to normalise array
 *
 * returns -1 if somethings goes wrong
 */
int write_png(const char* ofn, double *z, int nx, int ny, enum normalisation norm);

#endif
