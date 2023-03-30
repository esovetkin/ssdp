#ifndef _EDT_H
#define _EDT_H

/**
 * implements eucledian distance transform from the Meijster's paper
 *
 * @article{meijster2000general,
 *         title={A general algorithm for computing distance transforms in linear time},
 *         author={Meijster, Arnold and Roerdink, Jos BTM and Hesselink, Wim H},
 *         journal={Mathematical Morphology and its applications to image and signal processing},
 *         pages={331--340},
 *         year={2000},
 *         publisher={Springer}
 * }
 *
 * @x: stores index of the closes non-missing element
 *
 * @xg: stores index of the closes non-missing in one axis (first stage result)
 *
 * @G: stores the distance in one axis
 *
 * @z: **borrowed** array. memory management of z is done outside.
 *
 * @missing_value: upper limit, values below (or equal) are considered
 * missing.
 *
 * @s, @t: s and t arrays from the second stage
 */
struct edt {
        double *z;
        int *x, *xg, *G;
        int n, nx, ny;
        double missing_value;
        int *s, *t;
};


struct edt* edt_init(double *z, int n, int nx, int ny, double missing_value);
void edt_free(struct edt*);

void edt_compute(struct edt*);
void edt_fill(struct edt*);

#endif
