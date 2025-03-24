#include "long_preamble.h"
#include "short_preamble.h"
#include <stdio.h>
#include <math.h>
#include "data.h"
#include <complex.h>

void compute_normalized_correlation(const double complex Rx_signal[],
    int len_Rx,
    int delay_param,
    int window_length,
    double corr_out[],
    int out_length)
    
    {
    int n, k;
    // Temporary arrays for storing correlation sums and power sums
    double complex corr_arr[out_length];
    double peak_arr[out_length];
    
    // For each window starting at index n:
    for (n = 0; n < out_length; n++) {
        double complex sum_corr = 0.0 + 0.0 * I;
        double sum_peak = 0.0;
        for (k = 0; k < window_length; k++) {
            double complex sample1 = Rx_signal[n + k];
            double complex sample2 = Rx_signal[n + k + delay_param];
            // Cross-correlation: product of sample1 and sample2.
            sum_corr += sample1 * sample2;
            // Compute power (squared magnitude) of the delayed sample.
            sum_peak += pow(cabs(sample2), 2);
        }
        corr_arr[n] = sum_corr;
        peak_arr[n] = sum_peak;
    }

       // Compute normalized correlation: (|corr_arr|^2) / (peak_arr^2)
       for (n = 0; n < out_length; n++) {
        double mag_corr = cabs(corr_arr[n]);
        double mag_corr_sq = mag_corr * mag_corr;
        double peak_sq = peak_arr[n] * peak_arr[n];
        if (peak_sq != 0)
            corr_out[n] = mag_corr_sq / peak_sq;
        else
            corr_out[n] = 0.0;
    }
}