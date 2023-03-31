#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_DEPRECATE)
/* Functions like strcpy are technically not secure because they do */
/* not contain a 'length'. But we disable this warning for the VISA */
/* examples since we never copy more than the actual buffer size.   */
#define _CRT_SECURE_NO_DEPRECATE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <Windows.h>
#include <map>
//#include <omp.h>
#include <time.h> 
#include <chrono>

#include <C:\Users\user1\source\repos\owon\owon\WaveletTransform.h>
#include <C:\Users\user1\source\repos\owon\owon\OWON.h>

#include "C:\Program Files\IVI Foundation\VISA\Win64\Include\visa.h"

#pragma comment(lib, "C:/Program Files/IVI Foundation/VISA/Win64/Lib_x64/msc/visa64.lib")

using namespace std;
using namespace std::chrono;

//Program example
int main(void) {
	//Variables for getting data from OWON
	int sample_size = 4E3;
	double sample_step = 100E-9; //0 to min
	int sample_perSec_max = 125E6; //0 to default (125E6)
	double voltage_max = 0.1;
	int reading_number = 100;
	int trigger_level = 0; //from -5 to 5, where 5 equals voltage_max
	double* read_arr = new double[sample_size];
	string read_fileName = "read.txt"; //"" to skip writing

	for (int i = 0; i < sample_size; i++)
		read_arr[i] = 0;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	ViStatus status = OWON::getOWONData2(sample_size, sample_step, sample_perSec_max, voltage_max, trigger_level, reading_number, read_arr, read_fileName);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = time_span + duration_cast<duration<double>>(t2 - t1);
	std::cout << "It took me " << time_span.count() << " seconds.\n";
	
	//Variables for wavelet transform
	int t_sizing = 1; //If you need to use only every t_sizing element of samples

	//Frequency splitting
	int f_start = 0;
	int f_step = 1E3;
	int f_size = 5E2;

	//The wavelet that will be used in the transform
	wavelet_trans::wavelets wavelet = wavelet_trans::wavelets(6);

	//Let it be "" if you don't need to transform the signal
	//string transformedSignal_fileName = "";
	string transformedSignal_fileName = "signal_wavelet.txt";
	string backTransformedSignal_fileName = "signal_backWavelet.txt";


	if (transformedSignal_fileName!="") {
		double t_start, t_end, t_step;
		int t_size;

		//It is used in backWavelet() and in waveletTransform_fromTimeToTime()
		wavelet_trans::sizing_toDouble(0, sample_size, sample_step, sample_size, t_sizing, t_start, t_end, t_step, t_size);

		//Wavelet function generation
		complex<double>** wavelet_function = new complex<double>*[f_size];
		for (int i = 0; i < f_size; i++)
			wavelet_function[i] = new complex<double>[2 * t_size];
		make_waveletFunction_equalStep(t_start, t_step, t_size, f_start, f_step, f_size, wavelet, wavelet_function);

		//Wavelet transforms
		double** wavelet_transform = new double* [f_size];
		for (int i = 0; i < f_size; i++)
			wavelet_transform[i] = new double[t_size];
		
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		wavelet_trans::waveletTransform_fromIndexToIndex(read_fileName, sample_step, sample_size, 0, sample_size, t_sizing, f_start, f_step, f_size, wavelet_function, transformedSignal_fileName);
		//wavelet_trans::waveletTransform_fromTimeToTime(read_fileName, sample_step, sample_size, t_start, t_end, t_sizing, f_start, f_step, f_size, wavelet_function, transformedSignal_fileName);
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		printf("The signal has been transformed.\n");
		cout << time_span.count() << "s\n";

		if (backTransformedSignal_fileName!="") {
			t1 = high_resolution_clock::now();
			wavelet_trans::backWavelet(transformedSignal_fileName, t_start, t_step, t_size, f_start, f_step, f_size, wavelet_function, backTransformedSignal_fileName);
			t2 = high_resolution_clock::now();
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
			printf("The wavelet has been retransformed.\n");
			cout << time_span.count() << "s\n";

			if (transformedSignal_fileName!="" && backTransformedSignal_fileName!="") {
				double* axis_tau = new double[t_size];
				double* axis_s = new double[f_size];
				for (int i = 0; i < t_size; i++)
					axis_tau[i] = t_start + i * t_step;
				for (int i = 0; i < f_size; i++)
					axis_s[i] = f_start + i * f_step;

				wavelet_trans::printArray(axis_tau, t_size, "time_arr.txt");
				wavelet_trans::printArray(axis_s, f_size, "freq_arr.txt");
			}
		}
	}
	system("pause");
}