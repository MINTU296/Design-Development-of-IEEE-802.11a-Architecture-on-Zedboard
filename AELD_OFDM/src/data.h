/*
 * data.h
 *
 *  Created on: 11-Feb-2025
 *      Author: shree
 */

 #ifndef SRC_DATA_H_
 #define SRC_DATA_H_
 #include <complex.h>
 #define N_FFT 64 // FFT Size
 #define PI 3.14159265358979323846
 #define FC_HZ 5000000000 // Carrier frequency (5 GHz)
 #define FS_HZ 20000000   // Sampling frequency (20 MHz)
 #define TS_SEC (1.0 / FS_HZ) // Time period (50 ns)
 #define N_FFT 64        // FFT Size
 #define N_BITS 96       // Number of bits in one frame
 #define CHAR_COUNT (N_BITS / 8) //Number of characters that fit
 #define M 4             // QPSK modulation
 #define FRAME_TX_LEN 480    // Define the total length of Frame_Tx
 #define OVERSAMPLINGFACTOR 2 //Oversampling factor
 #define FRAME_TX_OVERSAMP_LENGTH (FRAME_TX_LEN  * OVERSAMPLINGFACTOR)
 #define NUM_SNR_VALUES 10
 #define NP_PACKETS_CAPTURE 3000 // Number of packets to capture

// Inputs for Packet Detection:
//   Rx_signal: input array of complex samples (length = LEN_RX)
//   len_Rx: length of Rx_signal (fixed as LEN_RX)
//   delay_param: delay applied for correlation
//   window_length: length of the window over which correlation is computed
//   out_length: number of windows (should be LEN_RX - delay_param + 1 - window_length)
#define LEN_RX 3000         // Length of Rx_signal array
#define DELAY_PARAM 16      // Delay parameter for correlation
#define WINDOW_LENGTH 32    // Correlation window length
#define OUT_LENGTH (LEN_RX - DELAY_PARAM + 1 - WINDOW_LENGTH) // Number of output windows 

#define SHORT_PREAMBLE_SLOT_LENGTH 16
#define PACKET_THRESHOLD 0.75
#define THRESHOLD_LENGTH 230



//**Defined the array elements in data.c
//// Short Preamble sequence (53 elements) **//
extern const complex double S_k[53];  // Declare the array using extern;
///Defined the array elements in data.c
//// Long Preamble sequence (53 elements)
extern const complex double L_k[53]; //Declare the array using extern;

 
extern const char text_1[12];  // Use 'extern' to declare but not define
extern const char text_2[12];

// Defining  the Tx RRC filter coefficients
extern const double rrc_filter_tx[21];
// Defining  the Rx RRC filter coefficients
extern const double rrc_filter_rx[21];



 #endif /* SRC_DATA_H_ */