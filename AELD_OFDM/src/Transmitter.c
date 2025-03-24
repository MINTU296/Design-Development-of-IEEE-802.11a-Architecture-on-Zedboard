#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "data.h"
#include "short_preamble.h"
#include "long_preamble.h"
#include "data_payload.h"
#include "correlation.h"
#include "pkt_selection.h"



//int sign_func(double x)
//{
//    if (x > 0)
//        return 1;
//    else if (x < 0)
//        return -1;
//    else
//        return 0;
//}
//
//// Structure to store results for each SNR index.
//typedef struct
//{
//    double evm;
//    double evm_AGC;
//    double ber;
//} Results;
//
//
//int main() {
//
//    complex double Short_preamble[160];
//    complex double Long_preamble[160];
//    int data_bits_1[N_BITS];
//    int data_bits_2[N_BITS];
//    complex double data_payload_1_mod[48], data_payload_2_mod[48];
//    complex double data_1_TX_payload[80], data_2_TX_payload[80];
//    complex double Frame_Tx[FRAME_TX_LEN];
//
// //************************************************************* */
//    // Call the function to create the short preamble
//    printf("\n Creating short preamble...\n");
//    create_short_preamble(Short_preamble);
////************************************************************* */
//
//    // Call the function to create the long preamble
//    printf("\n Creating long preamble...\n");
//    create_long_preamble(Long_preamble);
//
////************************************************************* */
//    // Create payload
//    printf("\n Creating two data_payloads for texts ...\n");
//    printf("\n  Text_1: %s\n", text_1);
//    printf("\n  Text_2: %s\n", text_2);
//
//    create_payload(data_bits_1, data_bits_2, data_payload_1_mod, data_payload_2_mod, data_1_TX_payload, data_2_TX_payload);
//
////************************************************************* */
//
//   // Concatenate the arrays into Frame_Tx using memcpy - a string function
//   printf("\n Creating Frame_Tx of length 480  ...\n");
//   create_frame_tx(Short_preamble,Long_preamble,data_1_TX_payload,data_2_TX_payload,Frame_Tx);
//
//
////************************************************************* */
//
//   //Oversampling:- Adding zeros in between each element.The number of zeros added is given by oversamplingFactor
//   printf("\n Creating Oversampled Frame_Tx of length 960  ...\n");
//   complex double Frame_Tx_Oversampled[FRAME_TX_OVERSAMP_LENGTH];
//
//   oversampleFrameTx(Frame_Tx_Oversampled,Frame_Tx);
////************************************************************* */
//
////Convolve the Tx frame with the RRC filter:-
//    printf("\n Convolving Frame_Tx with RRC filters  ...\n");
//   int rrc_filter_tx_length = sizeof(rrc_filter_tx) / sizeof(rrc_filter_tx[0]); //filterlength 21
//   int output_signal_len= (FRAME_TX_OVERSAMP_LENGTH + rrc_filter_tx_length - 1); //length of output signal ,980
//   complex double Tx_signal[output_signal_len]; //980 Tx_signal samples
//
//   convolveWithRRC(Tx_signal,Frame_Tx_Oversampled,rrc_filter_tx,rrc_filter_tx_length,output_signal_len);
//
//// ************************************************************* */
//
//   printf("\n Repeating the Frame_Tx of 980 samples  10 times  ...\n");
//   // Tx signal :- repeating 10 times
//    // Replicate the original Tx_signal 10 times
//    int Tx_signal_9800len= output_signal_len*10;
//
//    complex double Tx_signal_9800[Tx_signal_9800len];
//
//    repeat10times(output_signal_len, Tx_signal, Tx_signal_9800);
//  // ***********************************************************************************************/
//
//// Over The Air Transmission (Channel)
//printf("\n Sending it for Transmission over Air (Channel) with added gaussian noise  ...\n");
// complex double Tx_ota_signal[9800];
//
// double snr_dB_values[NUM_SNR_VALUES] = {8,9,10,11,12,13,14,15,16,17};
//
// // Calculate the average power of the transmitted signal
//double Tx_signal_power = 0.0;
//for (int i = 0; i < Tx_signal_9800len; i++) {
//   Tx_signal_power += creal(Tx_signal_9800[i] * conj(Tx_signal_9800[i])); // |Tx_signal|^2
//}
//Tx_signal_power /= output_signal_len;
//
//Results results[NUM_SNR_VALUES];
//
//for (int idx_snr_j = 0; idx_snr_j < NUM_SNR_VALUES; idx_snr_j++) {
//   double snr_dB = snr_dB_values[idx_snr_j];
//   double snr_linear = pow(10, snr_dB / 10.0); //Converting SNR from dB to linear scale
//   double noise_power = Tx_signal_power / snr_linear;
//   double noise_stddev = sqrt(noise_power/2);
//   for (int i = 0; i < Tx_signal_9800len; i++) {
//       double noise = noise_stddev * ((double)rand() / RAND_MAX ); // generating noise
//       Tx_ota_signal[i] = Tx_signal_9800[i] + noise;  //adding noise
//   }
//
//printf("Printing the output for SNR: %lf dB", snr_dB_values[idx_snr_j]);
//
////*****RECEIVER STARTS**************************/
//srand(0); // Equivalent to MATLAB's rng('default')
//int rx_start = rand() % (Tx_signal_9800len - NP_PACKETS_CAPTURE + 1); //generating a random variable between(0,9800-3000)
//
//
////Rx_signal variable which will hold the captured received signal
//complex double Rx_signal[NP_PACKETS_CAPTURE];
//printf("\n Randomly starting to capture packets from Channel ...\n");
////creating rx signal of 3000 samples from Tx_ota_signal
//for (int i = 0; i < NP_PACKETS_CAPTURE; i++) {
//    Rx_signal[i] = Tx_ota_signal[rx_start + i];
//}
//
//int Rx_signal_len = sizeof(Rx_signal) / sizeof(Rx_signal[0]);  //calculating length of rx_signal 3000
//
////printf("Rx_Signal_length: %d \n", Rx_signal_len);
//
//
// //Convolve the Rx frame with the RRC filter:-
// int rrc_filter_rx_length = sizeof(rrc_filter_rx) / sizeof(rrc_filter_rx[0]); //filterlength  21
// //printf("rrc_filter_rx_length: %d \n", rrc_filter_rx_length);
//
// int Rx_filtered_signal_len= (Rx_signal_len + rrc_filter_rx_length - 1); //length of output signal ,3020
// //printf("Rx_filtered_signal_len: %d \n", Rx_filtered_signal_len);
//
//  //Creating empty Rx_filtered_signal, for output
//  complex double Rx_filtered_signal[Rx_filtered_signal_len];
//  for (int i = 0; i < Rx_filtered_signal_len; i++) {
//      Rx_filtered_signal[i] = 0 + 0 * I;
//  }
//
//
//  //Creating the padding signal
//  //padding length for each side,ie the number of zeros to be added on each end side.
//   int pad_length_rx = rrc_filter_rx_length - 1;
//
//
//   int len_padded_rx = Rx_signal_len + 2 * pad_length_rx; // Total length of padded signal
// //  printf("len_padded_rx: %d \n", len_padded_rx);
//
//
//   //Signal for padding rx
//   complex double inp_sig_padded_rx[len_padded_rx];
//
//
//   // Fill the padded array with zeros then Rx_signal values and again zeros
//   for (int i = 0; i < rrc_filter_rx_length - 1; i++) {
//   inp_sig_padded_rx[i] = 0;
//   }
//   for (int i = 0; i < Rx_signal_len; i++) {
//   inp_sig_padded_rx[rrc_filter_tx_length - 1 + i] = Rx_signal[i];
//   }
//   for (int i = 0; i < rrc_filter_tx_length - 1; i++) {
//   inp_sig_padded_rx[rrc_filter_tx_length - 1 + FRAME_TX_OVERSAMP_LENGTH + i] = 0;
//   }
//
//
//
//     //Perform convolution between rcc coefff rx and padded signal rx
//    for (int n = 0; n < Rx_filtered_signal_len; n++) {
//         for (int k = 0; k < rrc_filter_rx_length; k++) {
//             Rx_filtered_signal[n] += rrc_filter_rx[k] * inp_sig_padded_rx[n + k]; // Multiply and accumulate
//              }
//    }
//     //end of convolution for receiver part.
//     //Now we have Rx_filtered_signal which is the result of convolution of coefficient and signal wiith padded zeroes.
//     //We have the Rx_filtered_signal with us
//
//
//
//
//
////******************************************************************************************************************************** */
//
////Variables needed are defined in header files
//// Prepare output array for normalized correlation.
//double corr_out[OUT_LENGTH];
//
////   corr_out: output array of normalized correlation values (length = out_length)
//// Compute normalized correlation using the function.
//compute_normalized_correlation(Rx_signal, LEN_RX, DELAY_PARAM, WINDOW_LENGTH, corr_out, OUT_LENGTH);
//
////******************************************************************************************************************************** */
//
//// Run packet selection.
//int packet_idx = packetSelection(corr_out, OUT_LENGTH);
//
//if (packet_idx != -1){
//    printf("Packet found!\n");
//    //printf("Packet start index: %d\n", packet_idx);
//}
//else
//    printf("No valid packet detected.\n");
//
//
//
//
//
////******************************************************************************************************************************** */
////Downsampling
//
////Creating a variable to store the rx_frame
//double complex rx_frame[480];
//
// // Calculate end index
// int end_index = OVERSAMPLINGFACTOR * FRAME_TX_LEN + packet_idx - 1;
// //printf("\n end_index = %d \n",end_index);
//// Extract values with step size oversampling_rate_rx
//int j = 0;
//for (int i = packet_idx; i <= end_index; i += OVERSAMPLINGFACTOR) {
//    rx_frame[j++] = Rx_filtered_signal[i];
//}
//
////******************************************************************************************************************************** */
//
//
//
////Coarse CFO Estimation
//int rx_frame_len = sizeof(rx_frame) / sizeof(rx_frame[0]); //filterlength
//double complex rx_frame_after_coarse[rx_frame_len];
//
//coarseCFOEstimation(rx_frame,rx_frame_after_coarse,rx_frame_len);
//
//  //******************************************************************************************************************************** */
//  //Fine Estimation
//  // Fine CFO Estimation
//  double complex rx_frame_after_fine[rx_frame_len];
//fineCFOEstimation(rx_frame_after_coarse,rx_frame_after_fine,rx_frame_len);
//
// //******************************************************************************************************************************** */
// // Channel Estimation
// complex double H_est_time[N_FFT];
//double complex H_est_used_for_fft[N_FFT];
//double complex H_est[64];
//
//channelEstimation(rx_frame_after_fine, H_est_used_for_fft, H_est, H_est_time);
//
////******************************************************************************************************************************** */
////ONE_TAP_EQUALIZER
// //One tap Equalizer
// complex double RX_Payload_1_Frequency[N_FFT];
// complex double RX_Payload_2_Frequency[N_FFT];
// complex double RX_Payload_1_Frequency_Equalizer[N_FFT];
// complex double RX_Payload_2_Frequency_Equalizer[N_FFT];
//
// oneTapEqualizer(rx_frame_after_fine,H_est,RX_Payload_1_Frequency,RX_Payload_2_Frequency, RX_Payload_1_Frequency_Equalizer,RX_Payload_2_Frequency_Equalizer);
//
////******************************************************************************************************************************** */
////DE_MAPPING
////De-Mapping
//double complex RX_Payload_1_no_Equalizer[48];
//double complex RX_Payload_1_no_pilot[48];
//double complex RX_Payload_2_no_Equalizer[48];
//double complex RX_Payload_2_no_pilot[48];
//demapping_RX_Payload(RX_Payload_1_Frequency, RX_Payload_1_Frequency_Equalizer,RX_Payload_1_no_Equalizer, RX_Payload_1_no_pilot);
//demapping_RX_Payload(RX_Payload_2_Frequency, RX_Payload_2_Frequency_Equalizer,RX_Payload_2_no_Equalizer, RX_Payload_2_no_pilot);
//
////END OF DE_MAPPING
////******************************************************************************************************************************** */
////AGC for Data_payload_1
////AGC For Rx Data Payload 1
//double complex RX_Payload_1_Final[48];
//double complex RX_Payload_2_Final[48];
//Rx_Payload_AGC(RX_Payload_1_Final,RX_Payload_1_no_pilot);
//Rx_Payload_AGC(RX_Payload_2_Final,RX_Payload_2_no_pilot);
//
////******************************************************************************************************************************** */
//
// //QPSK Demdoulation For Rx Data Payload 1
// int RX_Payload_1_Final_len = sizeof(RX_Payload_1_Final) / sizeof(RX_Payload_1_Final[0]);
// double  RX_Payload_1_demod[2 * RX_Payload_1_Final_len]; //96
// int RX_Payload_1_demod_len = sizeof(RX_Payload_1_demod) / sizeof(RX_Payload_1_demod[0]);
//
// int RX_Payload_2_Final_len = sizeof(RX_Payload_2_Final) / sizeof(RX_Payload_2_Final[0]);
// double  RX_Payload_2_demod[2 * RX_Payload_2_Final_len]; //96
// int RX_Payload_2_demod_len = sizeof(RX_Payload_2_demod) / sizeof(RX_Payload_2_demod[0]);
//
//  demodulation_Rx_Payload(RX_Payload_1_Final_len,RX_Payload_1_Final, RX_Payload_1_demod,RX_Payload_1_demod_len);
//  demodulation_Rx_Payload(RX_Payload_2_Final_len,RX_Payload_2_Final,RX_Payload_2_demod, RX_Payload_2_demod_len);
//
//char text_1_output[12];
//char text_2_output[12];
//
//decode_bits_to_text(RX_Payload_1_demod, text_1_output);
//printf("\n Decoded text_1: %s\n", text_1_output);
//
//decode_bits_to_text(RX_Payload_2_demod, text_2_output);
// printf("\n Decoded text_2: %s\n", text_2_output);
//
//
//
//    //END OF QPSK Demdoulation For Rx Data Payload 2
////******************************************************************************************************************************** */
//// EVM Calculation Before AGC
//
//        // Allocate an array for the concatenated error vector (size = 2*N)
//
//        double complex error_vector[N_BITS];
//
//        // Calculate error_vector = [RX_Payload_1_no_pilot, RX_Payload_2_no_pilot] - [data_payload_1_mod, data_payload_2_mod]
//        for (int i = 0; i < N_BITS/2; i++)
//        {
//            error_vector[i] = RX_Payload_1_no_pilot[i] - data_payload_1_mod[i];
//            error_vector[i + N_BITS/2] = RX_Payload_2_no_pilot[i] - data_payload_2_mod[i];
//        }
//
//    //     printf("Error Vector:\n");
//    // for (int i = 0; i < 10; i++) {
//    //     printf("Error_vector[%d]:(%f, %f)\n",i, creal(error_vector[i]), cimag(error_vector[i]));
//    // }
//
//        // Calculate sum of squared magnitudes for the error vector and for the transmitted symbols
//        double sum_error = 0.0;
//        double sum_ref = 0.0;
//        for (int i = 0; i < N_BITS/2; i++)
//        {
//            sum_error += pow(cabs(error_vector[i]), 2);
//            sum_ref += pow(cabs(data_payload_1_mod[i]), 2);
//        }
//        for (int i = 0; i < N_BITS/2; i++)
//        {
//            sum_error += pow(cabs(error_vector[i + N_BITS/2]), 2);
//            sum_ref += pow(cabs(data_payload_2_mod[i]), 2);
//        }
//
//        // Calculate EVM as the RMS value of error magnitude normalized by RMS value of transmitted symbols
//        double evm = sqrt(sum_error / N_BITS) / sqrt(sum_ref / N_BITS);
//
//        // Convert EVM to dB
//        double evm_dB = 20 * log10(evm);
//
//        // Print results
//        printf("EVM (linear) = %f\n", evm);
//        printf("EVM (dB)    = %f\n", evm_dB);
//
//
////******************************************************************************************************************************** */
//        // EVM Calculation After AGC
//
//        double complex error_vector_AGC[N_BITS];
//
//        // Calculate error_vector_AGC = [RX_Payload_1_Final, RX_Payload_2_Final] - [data_payload_1_mod, data_payload_2_mod]
//        for (int i = 0; i < N_BITS/2; i++)
//        {
//            error_vector_AGC[i] = RX_Payload_1_Final[i] - data_payload_1_mod[i];
//            error_vector_AGC[i + N_BITS/2] = RX_Payload_2_Final[i] - data_payload_2_mod[i];
//        }
//
//        // Calculate the sum of squared magnitudes for the error vector and the transmitted symbols
//        sum_error = 0.0;
//        sum_ref = 0.0;
//
//        for (int i = 0; i < N_BITS/2; i++)
//        {
//            sum_error += pow(cabs(error_vector_AGC[i]), 2);
//            sum_ref += pow(cabs(data_payload_1_mod[i]), 2);
//        }
//        for (int i = 0; i < N_BITS/2; i++)
//        {
//            sum_error += pow(cabs(error_vector_AGC[i + N_BITS/2]), 2);
//            sum_ref += pow(cabs(data_payload_2_mod[i]), 2);
//        }
//
//        // Calculate EVM as the RMS error normalized by the RMS of the transmitted symbols
//        double evm_AGC = sqrt(sum_error / N_BITS) / sqrt(sum_ref / N_BITS);
//
//        // Convert EVM to dB
//        double evm_AGC_dB = 20 * log10(evm_AGC);
//
//        // Print results
//        printf("EVM_AGC (linear) = %f\n", evm_AGC);
//        printf("EVM_AGC (dB)    = %f\n", evm_AGC_dB);
//
//// //******************************************************************************************************************************** */
//        // BER Calculation Before AGC
//
//        // Calculate the number of error bits.
//        int Error_bits = 0;
//        for (int i = 0; i < 96; i++)
//        {
//            // Using sign_func to check if there is a difference.
//            if (abs(sign_func(data_bits_1[i] - RX_Payload_1_demod[i])) == 1)
//            {
//                Error_bits++;
//            }
//        }
//        for (int i = 0; i < 96; i++)
//        {
//            if (abs(sign_func(data_bits_2[i] - RX_Payload_2_demod[i])) == 1)
//            {
//                Error_bits++;
//            }
//        }
//
//        // Calculate the Bit Error Rate (BER)
//        double BER = (double)Error_bits / (96 + 96);
//
//        // Adjust the size as needed.
//        results[idx_snr_j].evm = evm_dB;         // Store EVM (in dB) value
//        results[idx_snr_j].evm_AGC = evm_AGC_dB; // Store EVM after AGC (in dB) value
//        results[idx_snr_j].ber = BER;            // Store the Bit Error Rate
//        // Display the results.
//       // printf("Results[%d]: evm = %lf, evm_AGC = %lf, ber = %lf\n", idx_snr_j, results[idx_snr_j].evm, results[idx_snr_j].evm_AGC, results[idx_snr_j].ber);
//}
// //******************************************************************************************************************************** */
////Creating the variable array for plotting
//double evm_dB_values[NUM_SNR_VALUES];      // Equivalent to results.evm
//double evm_AGC_dB_values[NUM_SNR_VALUES];  // Equivalent to results.evm_AGC
//double ber_values[NUM_SNR_VALUES];         // Equivalent to results.ber
//
//// Initialize with some values (Replace with actual data)
//for (int i = 0; i < NUM_SNR_VALUES; i++) {
//    evm_dB_values[i] = results[i].evm;
//    evm_AGC_dB_values[i] = results[i].evm_AGC;
//    ber_values[i] = results[i].ber;
//}
//
//// Process the arrays
//for (int i = 0; i < NUM_SNR_VALUES; i++) {
//    // Replace -inf with -40
//    if (isinf(evm_AGC_dB_values[i]) && evm_AGC_dB_values[i] < 0) {
//        evm_AGC_dB_values[i] = -40;
//    }
//
//    // Replace 0 with 1e-6
//    if (ber_values[i] == 0) {
//        ber_values[i] = 1e-6;
//    }
//}
//
////  // Print results for verification
////  for (int i = 0; i < NUM_SNR_VALUES; i++) {
////     printf("EVM: %.2f, EVM_AGC: %.2f, BER: %.6f\n",
////            evm_dB_values[i], evm_AGC_dB_values[i], ber_values[i]);
//// }
//
////****************************************************************** */
////To copy to matlab for verification  - THE OUTPUT VALUES for EVM, BER
//printf("EVM_in_dB: \t");
//for (int i = 0; i < NUM_SNR_VALUES; i++) {
//    printf("%.2f", evm_dB_values[i]);  // Print the value
//
//    if (i < NUM_SNR_VALUES - 1) {
//        printf(" , ");  // Print comma and space only if it's not the last element
//    }
//}
//printf("\n");  // New line at the end
//
//printf("EVM_AGC_in_dB: \t");
//for (int i = 0; i < NUM_SNR_VALUES; i++) {
//    printf("%.2f", evm_AGC_dB_values[i]);  // Print the value
//
//    if (i < NUM_SNR_VALUES - 1) {
//        printf(" , ");  // Print comma and space only if it's not the last element
//    }
//}
//printf("\n");  // New line at the end
//
//printf("BER: \t");
//for (int i = 0; i < NUM_SNR_VALUES; i++) {
//    printf("%.2f", ber_values[i]);  // Print the value
//
//    if (i < NUM_SNR_VALUES - 1) {
//        printf(" , ");  // Print comma and space only if it's not the last element
//    }
//}
//printf("\n");  // New line at the end
//
//  //******************************************************************************************************************************** */
//
////Printing Output values one by one ****************************************************************************************************
//    // //Print the Short Preamble
//    // printf("Short Preamble:\n");
//    // for (int i = 0; i < 10; i++) {
//    //     printf("(%f, %f)\n", creal(Short_preamble[i]), cimag(Short_preamble[i]));
//    // }
//
//    // printf("\n");
//    // printf("\n");
//
//    //  // Print the Long Preamble
//    //  printf("Long Preamble:\n");
//    //  for (int i = 0; i < 10; i++) {
//    //      printf("(%f, %f)\n", creal(Long_preamble[i]), cimag(Long_preamble[i]));
//    //  }
//
//    // printf("\n");
//    // printf("\n");
//
//    //  // Print the Payload
//    // printf("Data 1 TX Payload:\n");
//    // for (int i = 0; i < 80; i++) {
//    //     printf("(%f, %f)\n", creal(data_1_TX_payload[i]), cimag(data_1_TX_payload[i]));
//    // }
//
//    // printf("\n");
//
//
//    // printf("\nData 2 TX Payload:\n");
//    // for (int i = 0; i < 80; i++) {
//    //     printf("(%f, %f)\n", creal(data_2_TX_payload[i]), cimag(data_2_TX_payload[i]));
//    // }
//
//    // // Print the final concatenated Frame_Tx
//    // printf("\nConcatenated Frame (Frame_Tx):\n");
//    // for (int i = 0; i < 100; i++) {
//    //     printf("Frame_Tx[%d] = (%f, %f)\n", i, creal(Frame_Tx[i]), cimag(Frame_Tx[i]));
//    // }
//
//    // printf("\n Oversampled Frame (Frame_Tx):\n");
//    // for (int i = 0; i < 100; i++) {
//    //     printf("Frame_Tx_Oversampled[%d] = (%f, %f)\n", i, creal(Frame_Tx_Oversampled[i]), cimag(Frame_Tx_Oversampled[i]));
//    // }
//
//    //  //Printing Input_padded signals
//    // printf("\nInput signals padded:\n");
//    // for (int i = 0; i < len_padded; i++) {
//    //     printf("Input_signals_padded[%d]: (%f, %f)\n", i, creal(inp_sig_padded[i]), cimag(inp_sig_padded[i]));
//    // }
//
//    //  //Printing Tx_signals after convolution
//    //  printf("\nTx_signals after convolution:\n");
//    //  for (int i = 0; i < 10; i++) {
//    //      printf("Tx_signal[%d]: (%f, %f)\n", i, creal(Tx_signal[i]), cimag(Tx_signal[i]));
//    //  }
//
//
////    //   *****Caution Prints out 9800 samples**********
//    // //Printing After repeating 10 times Tx_signals_9800 (It prints out only 10 )
//    //     printf("\nTx signals 9800 after repeating for 10 times:\n");
//    // for (int i = 0; i < 10; i++) {
//    //     printf("Tx_signal_9800[%d]: (%f, %f)\n", i, creal(Tx_signal_9800[i]), cimag(Tx_signal_9800[i]));
//    // }
//
//
//    // //Printing Tx_OTA_signals (Only first 100 signals)
//    //     printf("\nTx_OTA_signals:\n");
//    // for (int i = 0; i < 10; i++) {
//    //     printf("Tx_OTA_signals[%d]: (%f, %f)\n", i, creal(Tx_ota_signal[i]), cimag(Tx_ota_signal[i]));
//    // }
//
//
//    //     //Printing 3000 randomly captured signals from Channel
//    //     printf("\nReceiver captured signals_3000:\n");
//    // for (int i = 0; i < 10; i++) {
//    //     printf("Rx_captured_signals[%d]: (%f, %f)\n", i, creal(Rx_signal[i]), cimag(Rx_signal[i]));
//    // }
//
//// //Printing Rx_padded_signals 3040
////     printf("\nRx_padded_signals:\n");
//// for (int i = 0; i < 100; i++) {
////     printf("Rx_padded_signals[%d]: (%f, %f)\n", i, creal(inp_sig_padded_rx[i]), cimag(inp_sig_padded_rx[i]));
//// }
//
//// //Printing Rx_filtered_signals 3020
////     printf("\nRx_filtered_signals:\n");
//// for (int i = 0; i < 100; i++) {
////     printf("Rx_filtered_signals[%d]: (%f, %f)\n", i, creal(Rx_filtered_signal[i]), cimag(Rx_filtered_signal[i]));
//// }
//
//// //Printing Correlation output
////     printf("\nCorrelated_Outputs:\n"); //First few outputs
//// for (int i = 0; i < 10; i++) {
////     printf("Correlated_Output_signals[%d]: %lf \n", i, corr_out[i]);
//// }
//
//
//// //Printing Rx_frame after coarse
////     printf("\nRx frame after coarse:\n");
//// for (int i = 0; i < rx_frame_len; i++) {
////     printf("Rx_after_coarse[%d]: (%f, %f)\n", i, creal(rx_frame_after_coarse[i]), cimag(rx_frame_after_coarse[i]));
//// }
//
//// //Printing Rx_frame after fine estimation
//// printf("\nRx frame after fine estimation:\n");
//// for (int i = 0; i < rx_frame_len; i++) {
////     printf("Rx_after_fine[%d]: (%f, %f)\n", i, creal(rx_frame_after_fine[i]), cimag(rx_frame_after_fine[i]));
//// }
//
//// //Printing Long preamble_1 in Rx
//// printf("\nLong preamble_1 in Rx:\n");
//// for (int i = 0; i < N_FFT; i++) {
////     printf("Long_preamble_1[%d]: (%f, %f)\n", i, creal(Long_preamble_1[i]), cimag(Long_preamble_1[i]));
//// }
//
//// //Printing Long preamble_2 in Rx
//// printf("\nLong preamble_2 in Rx:\n");
//// for (int i = 0; i < N_FFT; i++) {
////     printf("Long_preamble_2[%d]: (%f, %f)\n", i, creal(Long_preamble_2[i]), cimag(Long_preamble_2[i]));
//// }
//
//
//// //Printing Long preamble_1_after_fft in Rx
//// printf("\nLong preamble_1_after_fft in Rx:\n");
//// for (int i = 0; i < N_FFT; i++) {
////     printf("Long_preamble_1_after_fft[%d]: (%f, %f)\n", i, creal(Long_preamble_1_After_FFT[i]), cimag(Long_preamble_1_After_FFT[i]));
//// }
//
//// //Printing Long preamble_2_after_fft in Rx
//// printf("\nLong preamble_2_after_fft in Rx:\n");
//// for (int i = 0; i < N_FFT; i++) {
////     printf("Long_preamble_2_after_fft[%d]: (%f, %f)\n", i, creal(Long_preamble_2_After_FFT[i]), cimag(Long_preamble_2_After_FFT[i]));
//// }
//
//// //Printing H_est_time after doing ifft
//// printf("\nH_est_time after doing ifft in Rx:\n");
//// for (int i = 0; i < N_FFT; i++) {
////     printf("H_est_time after ifft[%d]: (%f, %f)\n", i, creal(H_est_time[i]), cimag(H_est_time[i]));
//// }
//
//// //Printing RX_Payload_2_no_Equalizer
//// printf("\n RX_Payload_1_no_Equalizer in Rx:\n");
//// for (int i = 0; i < N_FFT; i++) {
////     printf("RX_Payload_1_no_Equalizer[%d]: (%f, %f)\n", i, creal(RX_Payload_1_no_Equalizer[i]), cimag(RX_Payload_1_no_Equalizer[i]));
//// }
//
//
//// //Printing RX_Payload_2_no_Equalizer
//// printf("\n RX_Payload_2_no_Equalizer in Rx:\n");
//// for (int i = 0; i < N_FFT; i++) {
////     printf("RX_Payload_2_no_Equalizer[%d]: (%f, %f)\n", i, creal(RX_Payload_2_no_Equalizer[i]), cimag(RX_Payload_2_no_Equalizer[i]));
//// }
//
//// //Printing Rx_payload_1_Frequency_Equalizer
//// printf("\nRx_payload_1_Frequency_Equalizer in Rx:\n");
//// for (int i = 0; i < N_FFT; i++) {
////     printf("Rx_payload_1_Frequency_Equalizer[%d]: (%f, %f)\n", i, creal(RX_Payload_1_Frequency_Equalizer[i]), cimag(RX_Payload_1_Frequency_Equalizer[i]));
//// }
//
//
//// //Printing Rx_payload_2_Frequency_Equalizer
//// printf("\nRx_payload_2_Frequency_Equalizer in Rx:\n");
//// for (int i = 0; i < N_FFT; i++) {
////     printf("Rx_payload_2_Frequency_Equalizer[%d]: (%f, %f)\n", i, creal(RX_Payload_2_Frequency_Equalizer[i]), cimag(RX_Payload_2_Frequency_Equalizer[i]));
//// }
//
//// //Printing Rx_payload_1_no_pilot
//// printf("\n Rx_payload_1_no_pilot in Rx:\n");
//// for (int i = 0; i < 48; i++) {
////     printf("Rx_payload_1_no_pilot[%d]: (%f, %f)\n", i, creal(RX_Payload_1_no_pilot[i]), cimag(RX_Payload_1_no_pilot[i]));
//// }
//
//
//// //Printing Rx_payload_2_no_pilot
//// printf("\n Rx_payload_2_no_pilot in Rx:\n");
//// for (int i = 0; i < 48; i++) {
////     printf("Rx_payload_2_no_pilot[%d]: (%f, %f)\n", i, creal(RX_Payload_2_no_pilot[i]), cimag(RX_Payload_2_no_pilot[i]));
//// }
//
//
//// //Printing Rx_payload_1_Final
//// printf("\n Rx_payload_1_Final in Rx:\n");
//// for (int i = 0; i < 48; i++) {
////     printf("Rx_payload_1_Final[%d]: (%f, %f)\n", i, creal(RX_Payload_1_Final[i]), cimag(RX_Payload_1_Final[i]));
//// }
//
//
//// //Printing Rx_payload_2_Final
//// printf("\n Rx_payload_2_Final in Rx:\n");
//// for (int i = 0; i < 48; i++) {
////     printf("Rx_payload_2_Final[%d]: (%f, %f)\n", i, creal(RX_Payload_2_Final[i]), cimag(RX_Payload_2_Final[i]));
//// }
//
//
//// //Printing Rx_payload_1_demod
//// printf("\n Rx_payload_1_demod in Rx:\n");
//// for (int i = 0; i < N_BITS; i++) {
////     printf("Rx_payload_1_demod[%d]: %lf\n", i, RX_Payload_1_demod[i]);
//// }
//
//
//// //Printing Rx_payload_2_demod
//// printf("\n Rx_payload_1_demod in Rx:\n");
//// for (int i = 0; i < N_BITS; i++) {
////     printf("Rx_payload_1_demod[%d]: %lf\n", i, RX_Payload_2_demod[i] );
//// }
//
//    return 0;
//}
