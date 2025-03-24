#include "data_payload.h"
#include "short_preamble.h"
#include <stdint.h>
#include "data.h"
#include <stdio.h>
#include <math.h>

// Function to generate random binary data
void generate_random_bits(int *data, int size) {
    for (int i = 0; i < size; i++) {
        data[i] = rand() % 2;
    }
}

//Function to generate bits from text data
void text_to_binary(const char *text, int *binary_data) {
    for (int i = 0; i < CHAR_COUNT; i++) {
        uint8_t ch = text[i];
        for (int j = 7; j >= 0; j--) {  // Extract bits from MSB to LSB
            binary_data[i * 8 + (7 - j)] = (ch >> j) & 1;
        }
    }
}

// Function to perform QPSK modulation
void qpsk_modulation(int *data_bits, complex double *modulated_symbols, int size) {
    for (int i = 0; i < size / 2; i++) {
        int bit1 = data_bits[2 * i];
        int bit2 = data_bits[2 * i + 1];

        if (bit1 == 0 && bit2 == 0)
            modulated_symbols[i] = (1 + I) / sqrt(2);
        else if (bit1 == 0 && bit2 == 1)
            modulated_symbols[i] = (-1 + I) / sqrt(2);
        else if (bit1 == 1 && bit2 == 0)
            modulated_symbols[i] = (-1 - I) / sqrt(2);
        else
            modulated_symbols[i] = (1 - I) / sqrt(2);
    }
}

// Function to construct the frame
void construct_frame(complex double *payload, complex double *frame, complex double *virtual_subcarrier, complex double *pilot) {
    for (int i = 0; i < 6; i++) frame[i] = virtual_subcarrier[i];
    for (int i = 0; i < 5; i++) frame[i + 6] = payload[i];
    frame[11] = pilot[0];
    for (int i = 0; i < 13; i++) frame[i + 12] = payload[i + 5];
    frame[25] = pilot[1];
    for (int i = 0; i < 6; i++) frame[i + 26] = payload[i + 18];
    frame[32] = 0;
    for (int i = 0; i < 6; i++) frame[i + 33] = payload[i + 24];
    frame[39] = pilot[2];
    for (int i = 0; i < 13; i++) frame[i + 40] = payload[i + 30];
    frame[53] = pilot[3];
    for (int i = 0; i < 5; i++) frame[i + 54] = payload[i + 43];
    for (int i = 0; i < 5; i++) frame[i + 59] = virtual_subcarrier[i + 6];
}

// Function to create the payload
void create_payload(int *data_bits_1, int *data_bits_2, complex double *data_payload_1, complex double *data_payload_2, complex double *data_1_TX_payload, complex double *data_2_TX_payload) {
    //int data_bits_1[N_BITS], data_bits_2[N_BITS];
    
    //complex double data_payload_1[48], data_payload_2[48];
    complex double virtual_subcarrier[11] = {0};
    complex double pilot[4] = {1, 1, 1, -1};

   
    
   //Converting the text to binary data
    text_to_binary(text_1, data_bits_1);
    text_to_binary(text_2,data_bits_2);

//    printf("Printing data_payload_1: \n");
//    for(int i=0;i<N_BITS;i++){
//     printf("data_bit_1[%d]: %d \n",i,data_bits_1[i]);
//    }
   
//    printf("Printing data_payload_2: \n");
//    for(int i=0;i<N_BITS;i++){
//     printf("data_bit_2[%d]: %d \n",i,data_bits_2[i]);
//    }


    // Perform QPSK modulation
    qpsk_modulation(data_bits_1, data_payload_1, N_BITS);
    qpsk_modulation(data_bits_2, data_payload_2, N_BITS);

    // Construct frames
    complex double data_frame_1[N_FFT], data_frame_2[N_FFT];
    construct_frame(data_payload_1, data_frame_1, virtual_subcarrier, pilot);
    construct_frame(data_payload_2, data_frame_2, virtual_subcarrier, pilot);

    // Apply FFT shift
    fftshift(data_frame_1, N_FFT);
    fftshift(data_frame_2, N_FFT);

    // Apply IFFT
    complex double data_frame_1_ifft[N_FFT], data_frame_2_ifft[N_FFT];
    ifft(data_frame_1, data_frame_1_ifft, N_FFT);
    ifft(data_frame_2, data_frame_2_ifft, N_FFT);

    // Add cyclic prefix
    for (int i = 0; i < 16; i++) {
        data_1_TX_payload[i] = data_frame_1_ifft[i + 48];
        data_2_TX_payload[i] = data_frame_2_ifft[i + 48];
    }
    for (int i = 0; i < 64; i++) {
        data_1_TX_payload[i + 16] = data_frame_1_ifft[i];
        data_2_TX_payload[i + 16] = data_frame_2_ifft[i];
    }
}

// Function to convert 96-bit double array into a string
void decode_bits_to_text(double bit_array[], char text[]) {
    for (int i = 0; i < 12; i++) {
        int ascii_value = 0;

        // Convert 8-bit segment to an ASCII character
        for (int j = 0; j < 8; j++) {
            ascii_value = (ascii_value << 1) | (int)bit_array[i * 8 + j]; // Cast double to int
        }

        text[i] = (char)ascii_value; // Store character in output string

        // // Debugging: Print each ASCII value
        // printf("Char %d: ASCII %d -> '%c'\n", i, ascii_value, text[i]);
    }

    text[12] = '\0'; // Null-terminate the string
}

void create_frame_tx(double complex *Short_preamble,double complex *Long_preamble, double complex *data_1_TX_payload , double complex *data_2_TX_payload, double complex *Frame_Tx){
    	   // Copy Short_preamble to Frame_Tx
           for (int i = 0; i < 160; i++) {
	        Frame_Tx[i] = Short_preamble[i];
	    }

	    // Copy Long_preamble to Frame_Tx
	    for (int i = 0; i < 160; i++) {
	        Frame_Tx[160 + i] = Long_preamble[i];
	    }

	    // Copy data_1_TX_payload to Frame_Tx
	    for (int i = 0; i < 80; i++) {
	        Frame_Tx[160 + 160 + i] = data_1_TX_payload[i];
	    }

	    // Copy data_2_TX_payload to Frame_Tx
	    for (int i = 0; i < 80; i++) {
	        Frame_Tx[160 + 160 + 80 + i] = data_2_TX_payload[i];
	    }
}

// Function to oversample Frame_Tx
void oversampleFrameTx(complex double *Frame_Tx_Oversampled, complex double *Frame_Tx) {
    for (int i = 0; i < FRAME_TX_OVERSAMP_LENGTH; i++) {
        if (i % OVERSAMPLINGFACTOR == 0) {
            Frame_Tx_Oversampled[i] = Frame_Tx[i / OVERSAMPLINGFACTOR]; // Copy original element
        } else {
            Frame_Tx_Oversampled[i] = 0 + 0 * I; // Insert zero (complex zero)
        }
    }
}

void convolveWithRRC(complex double *Tx_signal, complex double *Frame_Tx_Oversampled, const double *rrc_filter_tx, int filter_length, int output_signal_len){
     //Creating the padding signal 
   //padding length for each side,ie the number of zeros to be added on each end side.
   int pad_length = filter_length - 1; 
   int len_padded = FRAME_TX_OVERSAMP_LENGTH + 2 * pad_length; // Total length of padded signal 1000

   //signal for padding 
   complex double inp_sig_padded[len_padded];
   // Fill the padded array with zeros then Frame_Tx_Oversamp values and again zeros
   for (int i = 0; i < filter_length - 1; i++) {
   inp_sig_padded[i] = 0;
   }
   for (int i = 0; i < FRAME_TX_OVERSAMP_LENGTH; i++) {
   inp_sig_padded[filter_length - 1 + i] = Frame_Tx_Oversampled[i];
   }
   for (int i = 0; i < filter_length - 1; i++) {
   inp_sig_padded[filter_length - 1 + FRAME_TX_OVERSAMP_LENGTH + i] = 0;
   }

   //Initializing empty Tx_signal
for (int i = 0; i < output_signal_len; i++) {
    Tx_signal[i] = 0; 
}
   //perform convolution between rcc coefff and padded signal 
   for (int n = 0; n < output_signal_len; n++) {
    for (int k = 0; k < filter_length; k++) {
        Tx_signal[n] += rrc_filter_tx[k] * inp_sig_padded[n + k]; // Multiply and accumulate
         }
}
//end of convolution. 
//Now we have signal which is the result of convolution of coefficient and signal wiith padded zeroes. 980 - Tx_signals

}

void repeat10times(int output_signal_len, complex double *Tx_signal, complex double *Tx_signal_9800){
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < output_signal_len; j++) {
            Tx_signal_9800[i * output_signal_len + j] = Tx_signal[j];
        }
     }
}