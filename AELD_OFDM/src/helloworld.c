/******************************************************************************
*
* Copyright (C) 2009 - 2014 Xilinx, Inc.  All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* Use of the Software is limited solely to applications:
* (a) running on a Xilinx device, or
* (b) that interact with a Xilinx device through a bus or interconnect.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
* XILINX  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
* OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
* Except as contained in this notice, the name of the Xilinx shall not be used
* in advertising or otherwise to promote the sale, use or other dealings in
* this Software without prior written authorization from Xilinx.
*
******************************************************************************/

/*
 * helloworld.c: simple test application
 *
 * This application configures UART 16550 to baud rate 9600.
 * PS7 UART (Zynq) is not initialized by this application, since
 * bootrom/bsp configures it to baud rate 115200
 *
 * ------------------------------------------------
 * | UART TYPE   BAUD RATE                        |
 * ------------------------------------------------
 *   uartns550   9600
 *   uartlite    Configurable only in HW design
 *   ps7_uart    115200 (configured by bootrom/bsp)
 */

#include <stdio.h>
#include "platform.h"
#include "xil_printf.h"
#include <XTime_l.h>
#include "data.h"
#include "short_preamble.h"
#include "long_preamble.h"
#include "data_payload.h"
#include "correlation.h"
#include "pkt_selection.h"
#include "math.h"
int sign_func(double x)
{
    if (x > 0)
        return 1;
    else if (x < 0)
        return -1;
    else
        return 0;
}

// Structure to store results for each SNR index.
typedef struct
{
    double evm;
    double evm_AGC;
    double ber;
} Results;
void do_main_function(float *frame_tx_time,float *convolution_tx_time,float *total_tx_time,float *detect_time,float *select_time,float *estimation_time,float *total_rx_time, float *TOTAL_EXECUTION_TIME){
	  complex double Short_preamble[160];
	  complex double Long_preamble[160];
	  int data_bits_1[N_BITS];
	  int data_bits_2[N_BITS];
	  complex double data_payload_1_mod[48], data_payload_2_mod[48];
	  complex double data_1_TX_payload[80], data_2_TX_payload[80];
	  complex double Frame_Tx[FRAME_TX_LEN];
	  XTime Total_TX_Start, Total_TX_Stop;
	  XTime Frame_Start,Frame_Stop;
      XTime TOTAL_START, TOTAL_STOP;


//****Creating short_preamble, long_preamble, data_payloads and Frame_Tx*************//
    XTime_GetTime(&TOTAL_START);
    XTime_GetTime(&Total_TX_Start);
	XTime_GetTime(&Frame_Start);
	create_short_preamble(Short_preamble);
	create_long_preamble(Long_preamble);
	create_payload(data_bits_1, data_bits_2, data_payload_1_mod, data_payload_2_mod, data_1_TX_payload, data_2_TX_payload);
	create_frame_tx(Short_preamble,Long_preamble,data_1_TX_payload,data_2_TX_payload,Frame_Tx);
	XTime_GetTime(&Frame_Stop);
	*frame_tx_time = (float)1.0 * (Frame_Stop - Frame_Start) / (COUNTS_PER_SECOND/1000000);
//******Oversampling Frame_Tx  (480*2)****************************************//
	complex double Frame_Tx_Oversampled[FRAME_TX_OVERSAMP_LENGTH];
	oversampleFrameTx(Frame_Tx_Oversampled,Frame_Tx);

//********Convolution with RRC************************************************//
    XTime Convolution_TX_Start,Convolution_TX_Stop;
    XTime_GetTime(&Convolution_TX_Start);
	int rrc_filter_tx_length = sizeof(rrc_filter_tx) / sizeof(rrc_filter_tx[0]); //filterlength 21
	int output_signal_len= (FRAME_TX_OVERSAMP_LENGTH + rrc_filter_tx_length - 1); //length of output signal ,980
	complex double Tx_signal[output_signal_len]; //980 Tx_signal samples
	convolveWithRRC(Tx_signal,Frame_Tx_Oversampled,rrc_filter_tx,rrc_filter_tx_length,output_signal_len);
	XTime_GetTime(&Convolution_TX_Stop);
	*convolution_tx_time = (float)1.0 * (Convolution_TX_Stop - Convolution_TX_Start) / (COUNTS_PER_SECOND/1000000);
//********Repeating the convoluted Tx_signal(980 samples) 10 times to generate 9800 samples************************************************//
    int Tx_signal_9800len= output_signal_len*10;
	complex double Tx_signal_9800[Tx_signal_9800len];
	repeat10times(output_signal_len, Tx_signal, Tx_signal_9800);
    XTime_GetTime(&Total_TX_Stop);
    *total_tx_time  =  (float)1.0 * (Total_TX_Stop - Total_TX_Start) / (COUNTS_PER_SECOND/1000000);
//**********************Over The Air Transmission (Channel)*********************************************//
	 complex double Tx_ota_signal[9800];
     double snr_dB_values[NUM_SNR_VALUES] = {5,6,7,8,9,10,11,12,13,14};
	 // Calculate the average power of the transmitted signal
	double Tx_signal_power = 0.0;
	for (int i = 0; i < Tx_signal_9800len; i++) {
	   Tx_signal_power += creal(Tx_signal_9800[i] * conj(Tx_signal_9800[i])); // |Tx_signal|^2
	}
	Tx_signal_power /= output_signal_len;

	Results results[NUM_SNR_VALUES];
	for (int idx_snr_j = 0; idx_snr_j < NUM_SNR_VALUES; idx_snr_j++) {
	   double snr_dB = snr_dB_values[idx_snr_j];
	   double snr_linear = pow(10, snr_dB / 10.0); //Converting SNR from dB to linear scale
	   double noise_power = Tx_signal_power / snr_linear;
	   double noise_stddev = sqrt(noise_power/2);
	   for (int i = 0; i < Tx_signal_9800len; i++) {
	       complex double noise = noise_stddev * (((double)rand()/RAND_MAX) + ((double)rand()/RAND_MAX) * I); // generating noise
	       Tx_ota_signal[i] = Tx_signal_9800[i] + noise;  //adding noise
	   }

	   printf("Printing the output for SNR: %.1lf dB \n ", snr_dB_values[idx_snr_j]);
	   	//************************RECEIVER STARTS**************************/
	   XTime Total_Rx_Start, Total_Rx_Stop;
	   XTime_GetTime(&Total_Rx_Start);
	   	srand(2); // Equivalent to MATLAB's rng('default')
	   	int rx_start = rand() % (Tx_signal_9800len - NP_PACKETS_CAPTURE + 1); //generating a random variable between(0,9800-3000)
        //printf("rx_start_value_randomly generated: %d \n",rx_start);

	   	//Rx_signal variable which will hold the captured received signal
	   	complex double Rx_signal[NP_PACKETS_CAPTURE];
	   	//printf("\n Randomly starting to capture packets from Channel ...\n");
	   	//creating rx signal of 3000 samples from Tx_ota_signal
	   	for (int i = 0; i < NP_PACKETS_CAPTURE; i++) {
	   	    Rx_signal[i] = Tx_ota_signal[rx_start + i];
	   	}

	   	int Rx_signal_len = sizeof(Rx_signal) / sizeof(Rx_signal[0]);  //calculating length of rx_signal 3000

	   	//printf("Rx_Signal_length: %d \n", Rx_signal_len);


	   	 //Convolve the Rx frame with the RRC filter:-
	   	 int rrc_filter_rx_length = sizeof(rrc_filter_rx) / sizeof(rrc_filter_rx[0]); //filterlength  21
	   	 //printf("rrc_filter_rx_length: %d \n", rrc_filter_rx_length);

	   	 int Rx_filtered_signal_len= (Rx_signal_len + rrc_filter_rx_length - 1); //length of output signal ,3020
	   	 //printf("Rx_filtered_signal_len: %d \n", Rx_filtered_signal_len);

	   	  //Creating empty Rx_filtered_signal, for output
	   	  complex double Rx_filtered_signal[Rx_filtered_signal_len];
	   	  for (int i = 0; i < Rx_filtered_signal_len; i++) {
	   	      Rx_filtered_signal[i] = 0 + 0 * I;
	   	  }


	   	  //Creating the padding signal
	   	  //padding length for each side,ie the number of zeros to be added on each end side.
	   	   int pad_length_rx = rrc_filter_rx_length - 1;


	   	   int len_padded_rx = Rx_signal_len + 2 * pad_length_rx; // Total length of padded signal
	   	 //  printf("len_padded_rx: %d \n", len_padded_rx);


	   	   //Signal for padding rx
	   	   complex double inp_sig_padded_rx[len_padded_rx];


	   	   // Fill the padded array with zeros then Rx_signal values and again zeros
	   	   for (int i = 0; i < rrc_filter_rx_length - 1; i++) {
	   	   inp_sig_padded_rx[i] = 0;
	   	   }
	   	   for (int i = 0; i < Rx_signal_len; i++) {
	   	   inp_sig_padded_rx[rrc_filter_tx_length - 1 + i] = Rx_signal[i];
	   	   }
	   	   for (int i = 0; i < rrc_filter_tx_length - 1; i++) {
	   	   inp_sig_padded_rx[rrc_filter_tx_length - 1 + FRAME_TX_OVERSAMP_LENGTH + i] = 0;
	   	   }



	   	     //Perform convolution between rcc coefff rx and padded signal rx
	   	    for (int n = 0; n < Rx_filtered_signal_len; n++) {
	   	         for (int k = 0; k < rrc_filter_rx_length; k++) {
	   	             Rx_filtered_signal[n] += rrc_filter_rx[k] * inp_sig_padded_rx[n + k]; // Multiply and accumulate
	   	              }
	   	    }
	   	     //end of convolution for receiver part.
	   	     //Now we have Rx_filtered_signal which is the result of convolution of coefficient and signal with padded zeroes.
	   	     //We have the Rx_filtered_signal with us
//******************************************************************************************************************************** *//
// Prepare output array for normalized correlation.
//Packet Detection
	   	XTime Detect_Start, Detect_Stop;
	    XTime_GetTime(&Detect_Start);
	   	double corr_out[OUT_LENGTH];

	   	//   corr_out: output array of normalized correlation values (length = out_length)
	   	// Compute normalized correlation using the function.
	   	compute_normalized_correlation(Rx_signal, LEN_RX, DELAY_PARAM, WINDOW_LENGTH, corr_out, OUT_LENGTH);
	    XTime_GetTime(&Detect_Stop);
	    *detect_time = (float)1.0 * (Detect_Stop - Detect_Start) / (COUNTS_PER_SECOND/1000000);
//******************************************************************************************************************************** */
// Packet selection.
	   	XTime Select_Start, Select_Stop;
	    XTime_GetTime(&Select_Start);
	   	int packet_idx = packetSelection(corr_out, OUT_LENGTH);

	   	if (packet_idx != -1){
	   	    printf("Packet found!\n");
	   	    //printf("Packet start index: %d\n", packet_idx);
	   	}
	   	else
	   	    printf("No valid packet detected.\n");



	   	XTime_GetTime(&Select_Stop);
	   	*select_time = (float)1.0 * (Select_Stop - Select_Start) / (COUNTS_PER_SECOND/1000000);
//******************************************************************************************************************************** */
//Downsampling
//Creating a variable to store the rx_frame
	   	double complex rx_frame[480];

	   	 // Calculate end index
	   	 int end_index = OVERSAMPLINGFACTOR * FRAME_TX_LEN + packet_idx - 1;
	   	 //printf("\n end_index = %d \n",end_index);
	   	// Extract values with step size oversampling_rate_rx
	   	int j = 0;
	   	for (int i = packet_idx; i <= end_index; i += OVERSAMPLINGFACTOR) {
	   	    rx_frame[j++] = Rx_filtered_signal[i];
	   	}

//******************************************************************************************************************************** */
//Coarse CFO Estimation
	   	XTime Estimation_Start, Estimation_Stop;
	    XTime_GetTime(&Estimation_Start);
	   	int rx_frame_len = sizeof(rx_frame) / sizeof(rx_frame[0]); //filterlength
	   	double complex rx_frame_after_coarse[rx_frame_len];

	   	coarseCFOEstimation(rx_frame,rx_frame_after_coarse,rx_frame_len);

//******************************************************************************************************************************** */
// Fine CFO Estimation
	   	 double complex rx_frame_after_fine[rx_frame_len];
	   	fineCFOEstimation(rx_frame_after_coarse,rx_frame_after_fine,rx_frame_len);

//******************************************************************************************************************************** */
// Channel Estimation
	   	complex double H_est_time[N_FFT];
	   	double complex H_est_used_for_fft[N_FFT];
	   	double complex H_est[64];

	   	channelEstimation(rx_frame_after_fine, H_est_used_for_fft, H_est, H_est_time);
	    XTime_GetTime(&Estimation_Stop);
	    *estimation_time =  (float)1.0 * (Estimation_Stop - Estimation_Start) / (COUNTS_PER_SECOND/1000000);
//******************************************************************************************************************************** */
//ONE_TAP_EQUALIZER

	   	 complex double RX_Payload_1_Frequency[N_FFT];
	   	 complex double RX_Payload_2_Frequency[N_FFT];
	   	 complex double RX_Payload_1_Frequency_Equalizer[N_FFT];
	   	 complex double RX_Payload_2_Frequency_Equalizer[N_FFT];

	   	 oneTapEqualizer(rx_frame_after_fine,H_est,RX_Payload_1_Frequency,RX_Payload_2_Frequency, RX_Payload_1_Frequency_Equalizer,RX_Payload_2_Frequency_Equalizer);

//******************************************************************************************************************************** */
//DE_MAPPING

	   	double complex RX_Payload_1_no_Equalizer[48];
	   	double complex RX_Payload_1_no_pilot[48];
	   	double complex RX_Payload_2_no_Equalizer[48];
	   	double complex RX_Payload_2_no_pilot[48];
	   	demapping_RX_Payload(RX_Payload_1_Frequency, RX_Payload_1_Frequency_Equalizer,RX_Payload_1_no_Equalizer, RX_Payload_1_no_pilot);
	   	demapping_RX_Payload(RX_Payload_2_Frequency, RX_Payload_2_Frequency_Equalizer,RX_Payload_2_no_Equalizer, RX_Payload_2_no_pilot);

//END OF DE_MAPPING
//******************************************************************************************************************************** */
//AGC For Rx Data Payload 1
	   	double complex RX_Payload_1_Final[48];
	   	double complex RX_Payload_2_Final[48];
	   	Rx_Payload_AGC(RX_Payload_1_Final,RX_Payload_1_no_pilot);
	   	Rx_Payload_AGC(RX_Payload_2_Final,RX_Payload_2_no_pilot);
//******************************************************************************************************************************** */
//QPSK Demdoulation For Rx Data Payload 1
	   	 int RX_Payload_1_Final_len = sizeof(RX_Payload_1_Final) / sizeof(RX_Payload_1_Final[0]);
	   	 double  RX_Payload_1_demod[2 * RX_Payload_1_Final_len]; //96
	   	 int RX_Payload_1_demod_len = sizeof(RX_Payload_1_demod) / sizeof(RX_Payload_1_demod[0]);

	   	 int RX_Payload_2_Final_len = sizeof(RX_Payload_2_Final) / sizeof(RX_Payload_2_Final[0]);
	   	 double  RX_Payload_2_demod[2 * RX_Payload_2_Final_len]; //96
	   	 int RX_Payload_2_demod_len = sizeof(RX_Payload_2_demod) / sizeof(RX_Payload_2_demod[0]);

	   	  demodulation_Rx_Payload(RX_Payload_1_Final_len,RX_Payload_1_Final, RX_Payload_1_demod,RX_Payload_1_demod_len);
	   	  demodulation_Rx_Payload(RX_Payload_2_Final_len,RX_Payload_2_Final,RX_Payload_2_demod, RX_Payload_2_demod_len);
	   	XTime_GetTime(&Total_Rx_Stop);
	   	*total_rx_time = (float)1.0 * (Total_Rx_Stop - Total_Rx_Start) / (COUNTS_PER_SECOND/1000000);

	   	char text_1_output[12];
	   	char text_2_output[12];

	   	decode_bits_to_text(RX_Payload_1_demod, text_1_output);
	   	printf("\n Decoded text_1: %s\n", text_1_output);

	   	decode_bits_to_text(RX_Payload_2_demod, text_2_output);
	   	 printf("\n Decoded text_2: %s\n", text_2_output);
//******************************************************************************************************************************** */
// EVM Calculation Before AGC
// Allocate an array for the concatenated error vector (size = 2*N)

	   	double complex error_vector[N_BITS];

// Calculate error_vector = [RX_Payload_1_no_pilot, RX_Payload_2_no_pilot] - [data_payload_1_mod, data_payload_2_mod]
	    for (int i = 0; i < N_BITS/2; i++)
	    {
	   	    error_vector[i] = RX_Payload_1_no_pilot[i] - data_payload_1_mod[i];
	   	    error_vector[i + N_BITS/2] = RX_Payload_2_no_pilot[i] - data_payload_2_mod[i];
	    }

	   	//     printf("Error Vector:\n");
	   // for (int i = 0; i < 10; i++) {
	   //     printf("Error_vector[%d]:(%f, %f)\n",i, creal(error_vector[i]), cimag(error_vector[i]));
	   // }
// Calculate sum of squared magnitudes for the error vector and for the transmitted symbols
	   double sum_error = 0.0;
	   double sum_ref = 0.0;
	   for (int i = 0; i < N_BITS/2; i++)
	   {
	   	   sum_error += pow(cabs(error_vector[i]), 2);
	   	   sum_ref += pow(cabs(data_payload_1_mod[i]), 2);
	   }
	   for (int i = 0; i < N_BITS/2; i++)
	   {
	   	   sum_error += pow(cabs(error_vector[i + N_BITS/2]), 2);
	   	   sum_ref += pow(cabs(data_payload_2_mod[i]), 2);
	   }

// Calculate EVM as the RMS value of error magnitude normalized by RMS value of transmitted symbols
	   double evm = sqrt(sum_error / N_BITS) / sqrt(sum_ref / N_BITS);

	   // Convert EVM to dB
	   double evm_dB = 20 * log10(evm);
	   // Print results
	   //printf("EVM (dB)    = %f\n", evm_dB);
//******************************************************************************************************************************** */
	   // EVM Calculation After AGC

	   double complex error_vector_AGC[N_BITS];
// Calculate error_vector_AGC = [RX_Payload_1_Final, RX_Payload_2_Final] - [data_payload_1_mod, data_payload_2_mod]
	   for (int i = 0; i < N_BITS/2; i++)
	   {
	   	   error_vector_AGC[i] = RX_Payload_1_Final[i] - data_payload_1_mod[i];
	   	   error_vector_AGC[i + N_BITS/2] = RX_Payload_2_Final[i] - data_payload_2_mod[i];
	   }
// Calculate the sum of squared magnitudes for the error vector and the transmitted symbols
	   	    sum_error = 0.0;
	   	    sum_ref = 0.0;

	   	    for (int i = 0; i < N_BITS/2; i++)
	   	    {
	   	           sum_error += pow(cabs(error_vector_AGC[i]), 2);
	   	           sum_ref += pow(cabs(data_payload_1_mod[i]), 2);
	   	    }
	   	    for (int i = 0; i < N_BITS/2; i++)
	   	    {
	   	           sum_error += pow(cabs(error_vector_AGC[i + N_BITS/2]), 2);
	   	           sum_ref += pow(cabs(data_payload_2_mod[i]), 2);
	   	    }

// Calculate EVM as the RMS error normalized by the RMS of the transmitted symbols
	   	    double evm_AGC = sqrt(sum_error / N_BITS) / sqrt(sum_ref / N_BITS);
// Convert EVM to dB
	   	    double evm_AGC_dB = 20 * log10(evm_AGC);

// Print results
	   	   // printf("EVM_AGC (dB)    = %f\n", evm_AGC_dB);

//******************************************************************************************************************************** */
// BER Calculation Before AGC

	   	        // Calculate the number of error bits.
	   	        int Error_bits = 0;
	   	        for (int i = 0; i < 96; i++)
	   	        {
	   	            // Using sign_func to check if there is a difference.
	   	            if (abs(sign_func(data_bits_1[i] - RX_Payload_1_demod[i])) == 1)
	   	            {
	   	                Error_bits++;
	   	            }
	   	        }
	   	        for (int i = 0; i < 96; i++)
	   	        {
	   	            if (abs(sign_func(data_bits_2[i] - RX_Payload_2_demod[i])) == 1)
	   	            {
	   	                Error_bits++;
	   	            }
	   	        }

	   	        // Calculate the Bit Error Rate (BER)
	   	        double BER = (double)Error_bits / (96 + 96);

	   	        // Adjust the size as needed.
	   	        results[idx_snr_j].evm = evm_dB;         // Store EVM (in dB) value
	   	        results[idx_snr_j].evm_AGC = evm_AGC_dB; // Store EVM after AGC (in dB) value
	   	        results[idx_snr_j].ber = BER;            // Store the Bit Error Rate
	   	        // Display the results.
	   	       // printf("Results[%d]: evm = %lf, evm_AGC = %lf, ber = %lf\n", idx_snr_j, results[idx_snr_j].evm, results[idx_snr_j].evm_AGC, results[idx_snr_j].ber);
	   	}
	XTime_GetTime(&TOTAL_STOP);
	*TOTAL_EXECUTION_TIME = (float)1.0 * (TOTAL_STOP - TOTAL_START) / (COUNTS_PER_SECOND/1000000);
//******************************************************************************************************************************** */
//Creating the variable array for plotting
	   	double evm_dB_values[NUM_SNR_VALUES];      // Equivalent to results.evm
	   	double evm_AGC_dB_values[NUM_SNR_VALUES];  // Equivalent to results.evm_AGC
	   	double ber_values[NUM_SNR_VALUES];         // Equivalent to results.ber

	   	// Initialize with some values (Replace with actual data)
	   	for (int i = 0; i < NUM_SNR_VALUES; i++) {
	   	    evm_dB_values[i] = results[i].evm;
	   	    evm_AGC_dB_values[i] = results[i].evm_AGC;
	   	    ber_values[i] = results[i].ber;
	   	}

	   	// Process the arrays
	   	for (int i = 0; i < NUM_SNR_VALUES; i++) {
	   	    // Replace -inf with -40
	   	    if (isinf(evm_AGC_dB_values[i]) && evm_AGC_dB_values[i] < 0) {
	   	        evm_AGC_dB_values[i] = -40;
	   	    }

	   	    // Replace 0 with 1e-6
	   	    if (ber_values[i] == 0) {
	   	        ber_values[i] = 1e-6;
	   	    }
	   	}

	   	//  // Print results for verification
	   	//  for (int i = 0; i < NUM_SNR_VALUES; i++) {
	   	//     printf("EVM: %.2f, EVM_AGC: %.2f, BER: %.6f\n",
	   	//            evm_dB_values[i], evm_AGC_dB_values[i], ber_values[i]);
	   	// }

	   	//****************************************************************** */
	   	//To copy to matlab for verification  - THE OUTPUT VALUES for EVM, BER
	   	printf("EVM_in_dB: \t");
	   	for (int i = 0; i < NUM_SNR_VALUES; i++) {
	   	    printf("%.2f", evm_dB_values[i]);  // Print the value

	   	    if (i < NUM_SNR_VALUES - 1) {
	   	        printf(" , ");  // Print comma and space only if it's not the last element
	   	    }
	   	}
	   	printf("\n");  // New line at the end

	   	printf("EVM_AGC_in_dB: \t");
	   	for (int i = 0; i < NUM_SNR_VALUES; i++) {
	   	    printf("%.2f", evm_AGC_dB_values[i]);  // Print the value

	   	    if (i < NUM_SNR_VALUES - 1) {
	   	        printf(" , ");  // Print comma and space only if it's not the last element
	   	    }
	   	}
	   	printf("\n");  // New line at the end

	   	printf("BER: \t");
	   	for (int i = 0; i < NUM_SNR_VALUES; i++) {
	   	    printf("%.2f", ber_values[i]);  // Print the value

	   	    if (i < NUM_SNR_VALUES - 1) {
	   	        printf(" , ");  // Print comma and space only if it's not the last element
	   	    }
	   	}
	   	printf("\n");  // New line at the end

}






int main()
{
   init_platform();
    XTime PS_Start, PS_Stop;
    float TOTAL_EXECUTION_TIME = 0;
    float total_tx_time = 0;
    float frame_tx_time = 0;
    float convolution_tx_time = 0;
    float detect_time = 0;
    float select_time = 0;
    float estimation_time = 0;
    float total_rx_time = 0;
    XTime_SetTime(0);
    XTime_GetTime(&PS_Start);

    do_main_function(&frame_tx_time,&convolution_tx_time,&total_tx_time,&detect_time,&select_time,&estimation_time,&total_rx_time, &TOTAL_EXECUTION_TIME);

    XTime_GetTime(&PS_Stop);

    float time_processor = 0;

    time_processor = (float)1.0 * (PS_Stop - PS_Start) / (COUNTS_PER_SECOND/1000000);
    printf("Execution Time for Creating the Frame_TX  : %f ms \n" , frame_tx_time/1000);
    printf("Execution Time for Convolution in Tx  : %f ms \n" , convolution_tx_time/1000);
    printf("Execution Time for Transmitter_time : %f ms \n" , total_tx_time/1000);
    printf("Execution Time for Rx_Packet_Detection  : %f ms \n" , detect_time/(1000*NUM_SNR_VALUES));
    printf("Execution Time for Rx_Packet_Selection  : %f ms \n" , select_time/(1000*NUM_SNR_VALUES));
    printf("Execution Time for RX_Channel_Estimation  : %f ms \n" , estimation_time/(1000*NUM_SNR_VALUES));
    printf("Execution Time for Receiver_Time  : %f ms \n" , total_rx_time/(1000*NUM_SNR_VALUES));
    printf("Execution Time for Total_Execution_till_EVM_BER_Calculation  : %f ms \n" , TOTAL_EXECUTION_TIME/1000);



    printf("Execution Time for PS  : %f s \n" , time_processor/1000000);

    cleanup_platform();
    return 0;
}
