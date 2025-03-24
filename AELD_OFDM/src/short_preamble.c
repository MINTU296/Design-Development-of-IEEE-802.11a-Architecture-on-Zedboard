#include "short_preamble.h"
#include <stdio.h>
#include <math.h>
#include "data.h"

// Function to perform IFFT using naive DFT approach
void ifft(complex double *X, complex double *x, int N) {
    for (int n = 0; n < N; n++) {
        x[n] = 0.0 + 0.0 * I;
        for (int k = 0; k < N; k++) {
            double angle = 2 * PI * k * n / N;
            x[n] += X[k] * cexp(I * angle);
        }
        x[n] /= N;
    }
}

// Function to perform FFT shift (reordering the frequency bins)
void fftshift(complex double *X, int N) {
    int mid = N / 2;
    for (int i = 0; i < mid; i++) {
        complex double temp = X[i];
        X[i] = X[i + mid];
        X[i + mid] = temp;
    }
}

// Function to create the short preamble
void create_short_preamble(complex double *Short_preamble) {
    // Virtual subcarriers (padding zeros)
    complex double virtual_subcarrier[11] = {0};

    // Constructing the frequency-domain preamble (64 elements)
    complex double Short_preamble_slot_Frequency[N_FFT];
    
    for (int i = 0; i < 6; i++)
        Short_preamble_slot_Frequency[i] = virtual_subcarrier[i];
    
    for (int i = 0; i < 53; i++)
        Short_preamble_slot_Frequency[i + 6] = S_k[i] * sqrt(13.0 / 6);
    
    for (int i = 0; i < 5; i++)
        Short_preamble_slot_Frequency[i + 59] = virtual_subcarrier[i];

    // Perform IFFT shift
    fftshift(Short_preamble_slot_Frequency, N_FFT);

    // Time-domain representation of the short preamble
    complex double Short_preamble_slot_Time[N_FFT];
    ifft(Short_preamble_slot_Frequency, Short_preamble_slot_Time, N_FFT);

    // Generate the final Short Preamble sequence (repeat first 16 values 10 times)
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 16; j++)
            Short_preamble[i * 16 + j] = Short_preamble_slot_Time[j];
}

// Function to perform FFT using naive DFT approach
void fft(complex double *x, complex double *X, int N) {
    for (int k = 0; k < N; k++) {
        X[k] = 0.0 + 0.0 * I;
        for (int n = 0; n < N; n++) {
            double angle = -2 * PI * k * n / N;
            X[k] += x[n] * cexp(I * angle);
        }
    }
}
