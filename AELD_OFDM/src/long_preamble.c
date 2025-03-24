#include "long_preamble.h"
#include "short_preamble.h"
#include <stdio.h>
#include <math.h>
#include "data.h"

// Function to create the long preamble
void create_long_preamble(complex double *Long_preamble) {
    // Long Preamble sequence (53 elements)
    // double L_k[53] = {
    //     1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 
    //     0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1
    // };

    // Virtual subcarriers (padding zeros)
    complex double virtual_subcarrier[11] = {0};

    // Constructing the frequency-domain preamble (64 elements)
    complex double Long_preamble_slot_Frequency[N_FFT];
    
    for (int i = 0; i < 6; i++)
        Long_preamble_slot_Frequency[i] = virtual_subcarrier[i];
    
    for (int i = 0; i < 53; i++)
        Long_preamble_slot_Frequency[i + 6] = L_k[i];

    for (int i = 0; i < 5; i++)
        Long_preamble_slot_Frequency[i + 59] = virtual_subcarrier[i];

    // Perform IFFT shift
    fftshift(Long_preamble_slot_Frequency, N_FFT);

    // Time-domain representation of the long preamble
    complex double Long_preamble_slot_Time[N_FFT];
    ifft(Long_preamble_slot_Frequency, Long_preamble_slot_Time, N_FFT);

    // Generate the final Long Preamble sequence
    for (int i = 0; i < 32; i++)
        Long_preamble[i] = Long_preamble_slot_Time[i + 32];

    for (int i = 0; i < 64; i++)
        Long_preamble[i + 32] = Long_preamble_slot_Time[i];

    for (int i = 0; i < 64; i++)
        Long_preamble[i + 96] = Long_preamble_slot_Time[i];
}