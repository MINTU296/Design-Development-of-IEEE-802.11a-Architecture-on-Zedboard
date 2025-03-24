#ifndef CORRELATION_H
#define CORRELATION_H
#include "data.h"
#include <complex.h>  // Include necessary libraries


void compute_normalized_correlation (const double _Complex Rx_signal[],
    int len_Rx,
    int delay_param,
    int window_length,
    double corr_out[],
    int out_length);

#endif