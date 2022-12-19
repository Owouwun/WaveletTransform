#include <iostream>
#include "WaveletTrans.h"

using namespace wavelet_trans;

int main() {
	const double t_input_step = 16E-9;
	const int t_input_size = 1E5;

	int t_start_index = 0.125E5;
	int t_end_index = 0.625E5;
	//int t_start_index = 0;
	//int t_end_index = 1E5;
	int t_sizing = 100;

	int f_start = 0;
	int f_step = 1E3;
	int f_size = 5E2;

	wavelets wavelet = hermitian3;

	string transformingSignal_fileName = "input1.txt";
	string transformedSignal_fileName = "output1.txt";
	string backTransformedSignal_fileName = "output4.txt";

	threads_number = 8;

	double t_start = t_start_index * t_input_step;
	double t_end = t_end_index * t_input_step;
	double t_step = t_input_step * t_sizing;
	int t_size = (t_end_index - t_start_index) / t_sizing;

	complex<double>** wavelet_function = new complex<double>*[f_size];
	for (int i = 0; i < f_size; i++)
		wavelet_function[i] = new complex<double>[2 * t_size];
	make_waveletFunction_equalStep(t_start, t_step, t_size, f_start, f_step, f_size, wavelet, wavelet_function);
	
	double** wavelet_transform = new double*[f_size];
	for (int i = 0; i < f_size; i++)
		wavelet_transform[i] = new double[t_size];
	waveletTransform_fromIndexToIndex(transformingSignal_fileName, t_input_step, t_input_size, t_start_index, t_end_index, t_sizing, f_start, f_step, f_size, wavelet_function, transformedSignal_fileName);
	//waveletTransform_fromTimeToTime(transformingSignal_fileName, t_input_step, t_input_size, t_start, t_end, t_sizing, f_start, f_step, f_size, wavelet_function, transformedSignal_fileName);
	backWavelet(transformedSignal_fileName, t_start, t_step, t_size, f_start, f_step, f_size, wavelet_function, backTransformedSignal_fileName);
	
	double* axis_tau = new double[t_size];
	double* axis_s = new double[f_size];
	for (int i = 0; i < t_size; i++)
		axis_tau[i] = t_start + i * t_step;
	for (int i = 0; i < f_size; i++)
		axis_s[i] = f_start + i * f_step;

	printArray(axis_tau, t_size, "output2.txt");
	printArray(axis_s, f_size, "output3.txt");
}