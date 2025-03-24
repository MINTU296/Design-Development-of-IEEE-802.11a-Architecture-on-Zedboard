#ifndef PKT_SELECTION_H
#define PKT_SELECTION_H
#include "data.h"
#include <complex.h>
#include <stdlib.h>

int packetSelection(const double corr_out[], int corr_out_length);
void coarseCFOEstimation(complex double *rx_frame, complex double *rx_frame_after_coarse, int rx_frame_len);
void fineCFOEstimation(complex double *rx_frame_after_coarse, complex double *rx_frame_after_fine, int rx_frame_len);
void channelEstimation(complex double *rx_frame_after_fine, complex double *H_est_used_for_fft, complex double *H_est, complex double *H_est_time);
void oneTapEqualizer(complex double *rx_frame_after_fine, complex double *H_est, complex double *RX_Payload_1_Frequency, complex double *RX_Payload_2_Frequency, complex double *RX_Payload_1_Frequency_Equalizer, complex double *RX_Payload_2_Frequency_Equalizer);
void demapping_RX_Payload(complex double *RX_Payload_Frequency, complex double *RX_Payload_Frequency_Equalizer,complex double *RX_Payload_no_Equalizer,complex double *RX_Payload_no_pilot);
void Rx_Payload_AGC(complex double *RX_Payload_Final,complex double *RX_Payload_no_pilot);
void demodulation_Rx_Payload(int RX_Payload_Final_len, complex double *RX_Payload_Final, double *RX_Payload_demod,int RX_Payload_demod_len);

#endif // PKT_SELECTION_H