/* Copyright (c) 2007 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


#ifndef _SOBOLSEQ_H
#define _SOBOLSEQ_H


struct sobolseq {
		unsigned sdim;              /* dimension of sequence being generated */
		uint32_t *mdata;            /* array of length 32 * sdim */
		uint32_t *m[32];            /* more convenient pointers to mdata, of direction #s */
		uint32_t *x;                /* previous x = x_n, array of length sdim */
		unsigned *b;                /* position of fixed point in x[i] is after bit b[i] */
		uint32_t n;                 /* number of x's generated so far */
};


struct sobolseq* sobolseq_init(unsigned sdim);
void sobolseq_free(struct sobolseq* sd);
int sobolseq_gen(struct sobolseq * sd, double *x);

#endif
