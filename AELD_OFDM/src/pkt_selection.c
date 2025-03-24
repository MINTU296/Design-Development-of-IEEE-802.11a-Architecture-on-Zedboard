#include "long_preamble.h"
#include "short_preamble.h"
#include <stdio.h>
#include <math.h>
#include "data.h"

// Packet Selection Function
// Uses the normalized correlation output (corr_out) to identify the packet start index.
// Returns the packet start index (adjusted by LENGTH_RRC_RX) if found; otherwise returns -1.

int packetSelection(const double corr_out[], int corr_out_length)
{
    // Array to store indices where corr_out > PACKET_THRESHOLD.
    int packet_idx_arr[OUT_LENGTH];
    int idx_count = 0;
    
    // Loop through corr_out to save indices that meet the threshold.
    for (int i = 0; i < corr_out_length; i++) {
        if (corr_out[i] > PACKET_THRESHOLD) {
            packet_idx_arr[idx_count] = i;
            idx_count++;
            if (idx_count >= OUT_LENGTH)
                break;
        }
    }

    // Create a new static array to store the trimmed elements
    int trim_idx_arr[idx_count];

    // Copy only the first idx_count elements
    for (int i = 0; i < idx_count; i++) {
        trim_idx_arr[i] = packet_idx_arr[i];
    }
   
    int trimm_arr_size = sizeof(trim_idx_arr)/sizeof(trim_idx_arr[0]);
    //printf("\n Trimmed_Array_length: %d \n",sizeof(trim_idx_arr)/sizeof(trim_idx_arr[0]));
    int temp_1[trimm_arr_size + 1]; // Shifted right
    int temp_2[trimm_arr_size + 1]; // Shifted left
    int temp_3[trimm_arr_size + 1]; // Difference array

        // Create temp_1 (shift right, last element is 0)
        for (int i = 0; i < trimm_arr_size; i++) {
            temp_1[i] = packet_idx_arr[i];
        }
        temp_1[trimm_arr_size] = 0; // Last element is 0
    
        // Create temp_2 (shift left, first element is 0)
        temp_2[0] = 0; // First element is 0
        for (int i = 1; i <= trimm_arr_size; i++) {
            temp_2[i] = packet_idx_arr[i - 1];
        }
    
        // Compute temp_3 (element-wise subtraction)
        for (int i = 0; i <= trimm_arr_size; i++) {
            temp_3[i] = temp_1[i] - temp_2[i];
        }

        // printf("\n Temp3: %d \n",sizeof(temp_3)/sizeof(temp_3[0]));

int temp_3_size = sizeof(temp_3)/sizeof(temp_3[0]);


// Determine packet front positions by checking the difference between consecutive indices.
int packet_front[temp_3_size];
int packet_front_count = 0;

for (int i = 0; i < temp_3_size; i++) {
     if (temp_3[i] > 300) {
        packet_front[packet_front_count] = i;  // Save the position in packet_idx_arr.
        packet_front_count++;
    }
}
     // Create a new static array to store the trimmed elements
     int trim_packet_front_arr[packet_front_count];

    //printf("\n Trimmed_packet_front: %d \n",sizeof(trim_packet_front_arr)/sizeof(trim_packet_front_arr[0]));

    // Copy only the first idx_count elements
    for (int i = 0; i < packet_front_count; i++) {
        trim_packet_front_arr[i] = packet_front[i];
    }
    

    //Printing the trimmed_packet_front 3 variables
    // for(int i=0;i<packet_front_count;i++){
    //     printf("\n Trimmed_packet_front[%d]: %d \n",i,trim_packet_front_arr[i]);
    // }
   

    // Extract the actual packet front indices from packet_idx_arr.
    int packet_front_idx[packet_front_count];

    for (int i = 0; i < packet_front_count; i++) {
        packet_front_idx[i] = packet_idx_arr[ packet_front[i] ];
    }
    

    //Printing the packet_front_idx
    // for (int i = 0; i < packet_front_count; i++) {
    //     printf("\n packet_front_idx[%d]: %d \n",i,packet_front_idx[i]);
    // }
    


    int packet_idx = -1;
    // Loop through the packet front indices and check if the value at (index + THRESHOLD_LENGTH)
    // exceeds the threshold. Adjust the found index by LENGTH_RRC_RX.
    for (int x = 0; x < packet_front_count - 1; x++) {
        int index_check = packet_front_idx[x] + THRESHOLD_LENGTH;
        if (corr_out[index_check] > PACKET_THRESHOLD) {
            packet_idx = packet_front_idx[x] + 10 + 1;
            break;
        }
    }
    return packet_idx;
}

// Function for Coarse CFO Estimation and Correction
void coarseCFOEstimation(complex double *rx_frame, complex double *rx_frame_after_coarse, int rx_frame_len) {
    // Calculate the complex conjugate product
    complex double prod_consq_frame_coarse = 0.0;
    for (int i = 0; i < SHORT_PREAMBLE_SLOT_LENGTH; i++) {
        int idx1 = SHORT_PREAMBLE_SLOT_LENGTH * 5 + i;
        int idx2 = SHORT_PREAMBLE_SLOT_LENGTH * 6 + i;
        prod_consq_frame_coarse += rx_frame[idx1] * conj(rx_frame[idx2]);
    }
    
    // Estimate the coarse frequency offset
    double phase = atan2(cimag(prod_consq_frame_coarse), creal(prod_consq_frame_coarse));
    double freq_coarse_est = (-1.0 / (2.0 * PI * SHORT_PREAMBLE_SLOT_LENGTH * TS_SEC)) * phase;
    
    // Apply the coarse frequency offset to the received frame
    for (int n = 0; n < rx_frame_len; n++) {
        complex double correction = cexp(-I * 2.0 * PI * freq_coarse_est * TS_SEC * n);
        rx_frame_after_coarse[n] = rx_frame[n] * correction;
    }
}


// Function for Fine CFO Estimation and Correction
void fineCFOEstimation(complex double *rx_frame_after_coarse, complex double *rx_frame_after_fine, int rx_frame_len) {
    // Calculate the complex conjugate product for fine CFO estimation
    complex double prod_consq_frame_fine = 0.0;
    for (int i = 0; i < SHORT_PREAMBLE_SLOT_LENGTH * 4; i++) {
        int idx1 = SHORT_PREAMBLE_SLOT_LENGTH * 12 + i;
        int idx2 = SHORT_PREAMBLE_SLOT_LENGTH * 16 + i;
        prod_consq_frame_fine += rx_frame_after_coarse[idx1] * conj(rx_frame_after_coarse[idx2]);
    }
    
    // Estimate the fine frequency offset
    double phase = atan2(cimag(prod_consq_frame_fine), creal(prod_consq_frame_fine));
    double freq_fine_est = (-1.0 / (2.0 * PI * 64 * TS_SEC)) * phase;
    
    // Apply the fine frequency offset to the received frame
    for (int n = 0; n < rx_frame_len; n++) {
        complex double correction = cexp(-I * 2.0 * PI * freq_fine_est * TS_SEC * n);
        rx_frame_after_fine[n] = rx_frame_after_coarse[n] * correction;
    }
}

void channelEstimation(complex double *rx_frame_after_fine, complex double *H_est_used_for_fft, complex double *H_est, complex double *H_est_time){
    double complex Long_preamble_1[N_FFT]; // Output for Long_preamble_1
    double complex Long_preamble_2[N_FFT]; // Output for Long_preamble_2
   
    for (int i = 0; i < N_FFT; i++) {
        Long_preamble_1[i] = rx_frame_after_fine[SHORT_PREAMBLE_SLOT_LENGTH * 12 + i];
        Long_preamble_2[i] = rx_frame_after_fine[SHORT_PREAMBLE_SLOT_LENGTH * 16 + i];
    }
    
   
    double complex Long_preamble_1_After_FFT[N_FFT];
    double complex Long_preamble_2_After_FFT[N_FFT];
   
        // Perform FFT on Long_preamble_1
    fft(Long_preamble_1, Long_preamble_1_After_FFT, N_FFT);
    fftshift(Long_preamble_1_After_FFT, N_FFT);
    // Perform FFT on Long_preamble_2
    fft(Long_preamble_2, Long_preamble_2_After_FFT, N_FFT);
    fftshift(Long_preamble_2_After_FFT, N_FFT);
   
   //VALUES OF SLOT FREQUENCY 
   complex double Long_preamble_slot_Frequency[N_FFT];
   // Virtual subcarriers (padding zeros)
   complex double virtual_subcarrier[11] = {0};
   
   
   for (int i = 0; i < 6; i++)
   Long_preamble_slot_Frequency[i] = virtual_subcarrier[i];
   
   for (int i = 0; i < 53; i++)
   Long_preamble_slot_Frequency[i + 6] = L_k[i];
   
   for (int i = 0; i < 5; i++)
   Long_preamble_slot_Frequency[i + 59] = virtual_subcarrier[i];
   //DEFINED IN LONG_PREABLE, BUT NOT ABLE TO CALL IT 

// Perform the computation
for (int i = 0; i < 64; i++) {
H_est[i] = 0.5 * (Long_preamble_1_After_FFT[i] + Long_preamble_2_After_FFT[i]) * conj(Long_preamble_slot_Frequency[i]);
}




for(int i=0;i<N_FFT;i++){
H_est_used_for_fft[i] = H_est[i];
}


fftshift(H_est_used_for_fft, N_FFT);

ifft(H_est_used_for_fft, H_est_time, N_FFT);
}


void oneTapEqualizer(complex double *rx_frame_after_fine, complex double *H_est, complex double *RX_Payload_1_Frequency, complex double *RX_Payload_2_Frequency, complex double *RX_Payload_1_Frequency_Equalizer, complex double *RX_Payload_2_Frequency_Equalizer){
        // Step 1: Extract RX_Payload_1_time (elements 321 to 400)
        complex double RX_Payload_1_time[80];
        for (int i = 0; i < 80; i++) {
            RX_Payload_1_time[i] = rx_frame_after_fine[320 + i];
        }
        
        // Step 2: Remove cyclic prefix (first 16 elements)
        complex double RX_Payload_1_no_CP[N_FFT];
        for (int i = 0; i < N_FFT; i++) {
            RX_Payload_1_no_CP[i] = RX_Payload_1_time[16 + i];
        }
        
        // Step 3: Perform FFT and FFT shift
        
        fft(RX_Payload_1_no_CP, RX_Payload_1_Frequency, N_FFT);
        fftshift(RX_Payload_1_Frequency, N_FFT);
        
        // Step 4: Equalize the frequency-domain data
        for (int i = 0; i < N_FFT; i++) {
            RX_Payload_1_Frequency_Equalizer[i] = RX_Payload_1_Frequency[i] / H_est[i];
        }
        
        // Step 1: Extract RX_Payload_2_time (elements 401 to 480)
        complex double RX_Payload_2_time[80];
        for (int i = 0; i < 80; i++) {
            RX_Payload_2_time[i] = rx_frame_after_fine[400 + i];
        }
        
        // Step 2: Remove cyclic prefix (first 16 elements)
        complex double RX_Payload_2_no_CP[N_FFT];
        for (int i = 0; i < N_FFT; i++) {
            RX_Payload_2_no_CP[i] = RX_Payload_2_time[16 + i];
        }
        
        // Step 3: Perform FFT and FFT shift
        
        fft(RX_Payload_2_no_CP, RX_Payload_2_Frequency, N_FFT);
        fftshift(RX_Payload_2_Frequency, N_FFT);
        
        // Step 4: Equalize the frequency-domain data
        for (int i = 0; i < N_FFT; i++) {
            RX_Payload_2_Frequency_Equalizer[i] = RX_Payload_2_Frequency[i] / H_est[i];
        }
}

void demapping_RX_Payload(complex double *RX_Payload_Frequency, complex double *RX_Payload_Frequency_Equalizer,complex double *RX_Payload_no_Equalizer,complex double *RX_Payload_no_pilot){
    // Extract the required elements for (RX_Payload_1_no_Equalizer)
int idx = 0;
// Copy elements 7:11 (indices 6 to 10 in C)
for (int i = 6; i <= 10; i++) {
    RX_Payload_no_Equalizer[idx++] = RX_Payload_Frequency[i];
}
// Copy elements 13:25 (indices 12 to 24 in C)
for (int i = 12; i <= 24; i++) {
    RX_Payload_no_Equalizer[idx++] = RX_Payload_Frequency[i];
}
// Copy elements 27:32 (indices 26 to 31 in C)
for (int i = 26; i <= 31; i++) {
    RX_Payload_no_Equalizer[idx++] = RX_Payload_Frequency[i];
}
// Copy elements 34:39 (indices 33 to 38 in C)
for (int i = 33; i <= 38; i++) {
    RX_Payload_no_Equalizer[idx++] = RX_Payload_Frequency[i];
}
// Copy elements 41:53 (indices 40 to 52 in C)
for (int i = 40; i <= 52; i++) {
    RX_Payload_no_Equalizer[idx++] = RX_Payload_Frequency[i];
}
// Copy elements 55:59 (indices 54 to 58 in C)
for (int i = 54; i <= 58; i++) {
    RX_Payload_no_Equalizer[idx++] = RX_Payload_Frequency[i];
}


// Extract the required elements for (RX_Payload_1_no_pilot)
 idx = 0;

// Copy elements 7:11 (indices 6 to 10 in C)
for (int i = 6; i <= 10; i++) {
    RX_Payload_no_pilot[idx++] = RX_Payload_Frequency_Equalizer[i];
}

// Copy elements 13:25 (indices 12 to 24 in C)
for (int i = 12; i <= 24; i++) {
    RX_Payload_no_pilot[idx++] = RX_Payload_Frequency_Equalizer[i];
}

// Copy elements 27:32 (indices 26 to 31 in C)
for (int i = 26; i <= 31; i++) {
    RX_Payload_no_pilot[idx++] = RX_Payload_Frequency_Equalizer[i];
}

// Copy elements 34:39 (indices 33 to 38 in C)
for (int i = 33; i <= 38; i++) {
    RX_Payload_no_pilot[idx++] = RX_Payload_Frequency_Equalizer[i];
}

// Copy elements 41:53 (indices 40 to 52 in C)
for (int i = 40; i <= 52; i++) {
    RX_Payload_no_pilot[idx++] = RX_Payload_Frequency_Equalizer[i];
}

// Copy elements 55:59 (indices 54 to 58 in C)
for (int i = 54; i <= 58; i++) {
    RX_Payload_no_pilot[idx++] = RX_Payload_Frequency_Equalizer[i];
}
}

void Rx_Payload_AGC(complex double *RX_Payload_Final,complex double *RX_Payload_no_pilot){
// Initialize RX_Payload_1_Final with zeros
for (int i = 0; i < 48; i++) {
    RX_Payload_Final[i] = 0 + 0 * I; // Initialize to zero
}

// Process each element in RX_Payload_1_no_pilot
for (int idx = 0; idx < 48; idx++) {
    // Get the current point
    double complex RX_Payload_1_no_pilot_current = RX_Payload_no_pilot[idx];

    // Initialize Mapped Real and Imaginary Parts
    double RX_Payload_1_mapped_real = 0.0; // Real part of the mapped point
    double RX_Payload_1_mapped_imag = 0.0; // Imaginary part of the mapped point

    // Decision Logic for Real Part
    if (creal(RX_Payload_1_no_pilot_current) > 0) {
        RX_Payload_1_mapped_real = 1.0 / sqrt(2.0); // Map to +1/sqrt(2) for positive real part
    } else if (creal(RX_Payload_1_no_pilot_current) < 0) {
        RX_Payload_1_mapped_real = -1.0 / sqrt(2.0); // Map to -1/sqrt(2) for negative real part
    }

    // Decision Logic for Imaginary Part
    if (cimag(RX_Payload_1_no_pilot_current) > 0) {
        RX_Payload_1_mapped_imag = 1.0 / sqrt(2.0); // Map to +1/sqrt(2) for positive imaginary part
    } else if (cimag(RX_Payload_1_no_pilot_current) < 0) {
        RX_Payload_1_mapped_imag = -1.0 / sqrt(2.0); // Map to -1/sqrt(2) for negative imaginary part
    }

    // Combine Mapped Real and Imaginary Parts
    RX_Payload_Final[idx] = RX_Payload_1_mapped_real + I * RX_Payload_1_mapped_imag;
}
}

void demodulation_Rx_Payload(int RX_Payload_Final_len, complex double *RX_Payload_Final, double *RX_Payload_demod,int RX_Payload_demod_len){
    for (int i = 0; i < RX_Payload_demod_len; i++) {
        RX_Payload_demod[i] = 0 + 0 * I; // Initialize to zero (Preallocate the demodulated bits array for better performance)
    }
    // Loop through each symbol
    for (int i = 0; i < RX_Payload_Final_len; i++) {
        // Extract real and imaginary parts of the current symbol
        double RX_Payload_1_demod_real = creal(RX_Payload_Final[i]);
        double RX_Payload_1_demod_imag = cimag(RX_Payload_Final[i]);

        // Determine bits based on QPSK mapping
        if (RX_Payload_1_demod_real > 0 && RX_Payload_1_demod_imag > 0) {
            // Mapping for 00 -> 0.707 + i*0.707
            RX_Payload_demod[2 * i] = 0;
            RX_Payload_demod[2 * i + 1] = 0;
        } else if (RX_Payload_1_demod_real < 0 && RX_Payload_1_demod_imag > 0) {
            // Mapping for 01 -> -0.707 + i*0.707
            RX_Payload_demod[2 * i] = 0;
            RX_Payload_demod[2 * i + 1] = 1;
        } else if (RX_Payload_1_demod_real < 0 && RX_Payload_1_demod_imag < 0) {
            // Mapping for 10 -> -0.707 - i*0.707
            RX_Payload_demod[2 * i] = 1;
            RX_Payload_demod[2 * i + 1] = 0;
        } else if (RX_Payload_1_demod_real > 0 && RX_Payload_1_demod_imag < 0) {
            // Mapping for 11 -> 0.707 - i*0.707
            RX_Payload_demod[2 * i] = 1;
            RX_Payload_demod[2 * i + 1] = 1;
        }
    }
}