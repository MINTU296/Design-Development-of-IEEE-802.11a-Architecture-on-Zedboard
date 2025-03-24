#ifndef DATA_PAYLOAD_H
#define DATA_PAYLOAD_H
#include "data.h"
#include <complex.h>
#include <stdlib.h>

void generate_random_bits(int *data, int size);
void text_to_binary(const char *text, int *binary_data);
void qpsk_modulation(int *data_bits, complex double *modulated_symbols, int size);
void create_payload(int *data_bits_1, int *data_bits_2, complex double *data_payload_1, complex double *data_payload_2, complex double *data_1_TX_payload, complex double *data_2_TX_payload);
void construct_frame(complex double *payload, complex double *frame, complex double *virtual_subcarrier, complex double *pilot);
void decode_bits_to_text(double bit_array[], char text[]);
void create_frame_tx(double complex *Short_preamble,double complex *Long_preamble, double complex *data_1_TX_payload , double complex *data_2_TX_payload, double complex *Frame_Tx);
void oversampleFrameTx(complex double *Frame_Tx_Oversampled, complex double *Frame_Tx) ;
void convolveWithRRC(complex double *Tx_signal, complex double *Frame_Tx_Oversampled,const double *rrc_filter_tx, int filter_length,int output_signal_len);
void repeat10times(int output_signal_len, complex double *Tx_signal, complex double *Tx_signal_9800);
#endif // DATA_PAYLOAD_H