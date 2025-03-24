#ifndef SHORT_PREAMBLE_H
#define SHORT_PREAMBLE_H

#include <complex.h>  // Include necessary libraries


void create_short_preamble(complex double *Short_preamble);
void ifft(complex double *X, complex double *x, int N);
void fftshift(complex double *X, int N);
void fft(complex double *x, complex double *X, int N);
#endif // SHORT_PREAMBLE_H
