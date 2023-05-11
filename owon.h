#pragma once

#ifdef OWON_EXPORTS
#define OWON_API __declspec(dllexport)
#else
#define OWON_API __declspec(dllimport)
#endif

#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_DEPRECATE)
/* Functions like strcpy are technically not secure because they do */
/* not contain a 'length'. But we disable this warning for the VISA */
/* examples since we never copy more than the actual buffer size.   */
#define _CRT_SECURE_NO_DEPRECATE
#endif

#include "pch.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <Windows.h>
#include <map>
#include <chrono>
#include <omp.h>
#include <complex>

#include "common.h"

#include "C:\Program Files\IVI Foundation\VISA\Win64\Include\visa.h"

#pragma comment(lib, "C:/Program Files/IVI Foundation/VISA/Win64/Lib_x64/msc/visa64.lib")

using namespace std;
using namespace std::chrono;

namespace OWON {
	static ViSession defaultRM;
	static ViSession instr;
	static ViUInt32 numInstrs;
	static ViFindList findList;
	static ViUInt32 retCount;
	static ViUInt32 writeCount;
	static ViStatus status;
	static ViStatus status1;
	static ViStatus status2;
	static ViChar instrResourceString[VI_FIND_BUFLEN];

	static char* buffer;
	static char stringinput[512];

	//����� ������ ������
	enum COUPLING {
		COUPLING_AC,
		COUPLING_DC
	};

	//������������ ���������������
	enum VERTICAL_SCALE {
		VS2mV = 2,
		VS5mV = 5,
		VS10mV = 10,
		VS20mV = 20,
		VS50mV = 50,
		VS100mV = 100,
		VS200mV = 200,
		VS500mV = 500,
		VS1V = 1000,
		VS2V = 2000,
		VS5V = 5000
	};

	enum HORIZONTAL_SCALE {
		//HS1ns,
		//HS2ns,
		HS5ns = 5,
		HS10ns = 10,
		HS20ns = 20,
		HS50ns = 50,
		HS100ns = 100,
		HS200ns = 200,
		HS500ns = 500,
		HS1us = 1000,
		HS2us = 2000,
		HS5us = 5000,
		HS10us = 10000,
		HS20us = 20000,
		HS50us = 50000,
		HS100us = 100000,
		HS200us = 200000,
		HS500us = 500000,
		HS1ms = 1000000,
		HS2ms = 2000000,
		HS5ms = 5000000,
		HS10ms = 10000000,
		HS20ms = 20000000,
		HS50ms = 50000000,
		HS100ms = 100000000,
		HS200ms = 200000000,
		HS500ms = 500000000,
		HS1s = 1000000000,
		HS2s = 2000000000,
		HS5s = 5000000000,
		HS10s = 10000000000,
		HS20s = 20000000000,
		HS50s = 50000000000,
		HS100s = 100000000000
	};

	enum AQUIRE_DEPMEM {
		AD1K,
		AD10K,
		AD100K,
		AD1M,
		AD10M,
		AD100M
	};

	void vertScale(VERTICAL_SCALE sc, float& scCoef, const char*& chSc);

	string vertScale_to_string(VERTICAL_SCALE sc);

	string horiScale_to_string(HORIZONTAL_SCALE hs);
	double horiScale_to_number(HORIZONTAL_SCALE hs);
	HORIZONTAL_SCALE number_to_horiScale(double& hs);

	const char* acqDepmem(AQUIRE_DEPMEM ad);

	const char* coup(COUPLING c);

	VERTICAL_SCALE voltageMax_to_vertScale(double& vm);

	AQUIRE_DEPMEM size_to_acquire(int& sample_size, int& sample_acquire);

	void doLog(string msg);

	extern "C" __declspec(dllexport) int getOWONData2(int sample_size, double sample_step, int sample_perSec_max, double voltage_max, int trigger_level, int reading_number, double* result_arr, const char* result_fileName);
}

namespace Wavelet {
	// ���������� �������, ���������� �� ������������ ����������
	const double h3cst = 1. / (2. / 15. * sqrt(30.) * pow(PI, -0.25));

	// ������ �������������� �������-�������
	enum wavelets
	{
		hermitian1 = 0, // ����� ������� �� ����������� ���������
		hermitian2 = 1, // ������������� �������� "������������ �����".
		hermitian3 = 2, // ������� ���� ����� ���� �� ����� ������.
		hermitian4 = 3,
		poisson2 = 4, // ������� ���� ���� ���� �� ����� ������.
		modified_morlet = 6,
		morlet = 5,
	};

	enum window
	{
		none,
		square,
		triangle,
		RCF // Raised Cosine Filter
	};

	// ��������� �������������� �������� ������� f � ����� step ������� size.
	// ����� ��������������: ��������� ������� �������� (������� ������)
	// ��������! ���� �������� ��������� ������ �� ������� ����������� ��������� (����� ����� inf), ������� 0.
	complex<double> kotes(complex<double>* f, double step, int size);

	double kotes(double* f, double step, int size);

	// ������ ��������� �������-��������������.
	// � ��������������� ����������� ������� ����� �� �������� ���������. ��������� ������ ������������ ������������ ����, �������� ������� ������������ ��������� �������������� ��������� ���.
	// ��� ������� ������ �������� ������ ��� "�������������" ��� "���������" ������������ ��� ������.
	// ��� ������ ������������� �������� ������������� ��������� ������� ������� (3) � ���� ��������� ������������� �� �������. �� ���� ��� ����� ������ ����� ��������� ������� ��� ��������� ��������.
	// (���������� ��������� ������� ������ ���������, �.�. ����� ������ ��������� ������� �� ���������� ��������� ������� �����������, � �� ����� ��� �� ��������� - ���� �������)
	complex<double> morlet_wavelet(double t);
	double modified_morlet_wavelet(double t);
	double hermitian1_wavelet(double t);
	double hermitian2_wavelet(double t);
	double hermitian3_wavelet(double t);
	double hermitian4_wavelet(double t);
	double poisson2_wavelet(double t);

	extern "C" __declspec(dllexport) void sizing_toDouble(
		int t_input_startIndex,
		int t_input_endIndex,
		double t_input_step,
		int t_input_size,
		int t_sizing,
		double* t_sizing_start,
		double* t_sizing_end,
		double* t_sizing_step,
		int* t_sizing_size
	);

	// ����������� �������-������� (������������ ��������) res ���� wavelet ��� ������� �� t_start-t_step*t_size �� t_start+t_step*t_size � ����� t_step � ������ �� f_start �� f_start+f_step*f_size c ����� f_step.
	// ������� f ������� � ��������� s, ������������ � ������������ �������� ������������ �������-��������������, ������������ f=1/s.
	// ��������! � �������� �������-������� ������� ��������� 1/sqrt(f). � ������������ �������� ������� �� �� �����������, ������ ��� ��� �� ������ ������������ � ������ � �������� �������-��������������, �� ��� ������� � ���� �������-�������.
	// ��������! ����������� ������� "�������" ������������ t_start �� ������� � ����� ����������� ����������: ������ ����, ����� ����� ������� ������ ���� [f_size x t_size x t_size], ���������� ����� ������ [f_size x 2*t_size], �.�. tau (������������ �������) ����� ��� �� ��� �� ���������, ��� � �����.
	void make_waveletFunction_equalStep(double t_start, double t_step, int t_size, double f_start, double f_step, int f_size, wavelets wavelet, complex<double>* res);
	extern "C" __declspec(dllexport) void make_waveletFunction_equalStep(double t_start, double t_step, int t_size, double f_start, double f_step, int f_size, wavelets wavelet, double* res);

	// ������ �������-�������������� �������, �������� ������� ����� �� ����� inputFileName, ����� ��������� ��� t_input_step � ���������� �������� t_input_size;
	// � �������������� ������������ ������ t_sizing ��������, ������� � t_start_index � ���������� t_end_index;
	// �� ������ ���������� ������� ����������� [f_size x t_end_index-t_start_index)/sizing], ���������������
	// �������� � f_start �� f_start+f_step*f_size � ����� f_step
	// � �������� � t_input_step*t_start_index �� t_input_step*t_end_index � ����� t_input_step*t_sizing.
	// ����������� ������� ����� �������� psitaus.
	// �������� �������������� ������������ � ���� outputFile_name.
	// ������� f ������� � ��������� s, ������������ � ������������ �������� ������������ �������-��������������, ������������ f=1/s.
	// �������� �������� �� ����������� �������, ��������� ����������� �������! ��������� � �������� ������� ���������� ������������ �������� make_waveletFunction_equalStep.

	// Double
	extern "C" __declspec(dllexport) void waveletTransform_fromIndexToIndex(
		char* inputFileName,
		double t_input_step, double t_input_size,
		int t_start_index, int t_end_index, int t_sizing,
		double f_start, double f_step, int f_size,
		double* psitaus,
		char* outputFile_name,
		window win,
		double* win_top_t,
		double* win_top_f,
		int win_top_n,
		double* win_bot_t,
		double* win_bot_f,
		int win_bot_n);

	// �������, ����������� waveletTransform_fromIndexToIndex, �� ����� ������������ �� ��������� � �������� ���������, � ���������� � ��������.
	extern "C" __declspec(dllexport) void waveletTransform_fromTimeToTime(
		char* inputFileName,
		double t_input_step, double t_input_size,
		int t_start, int t_end, int t_sizing,
		double f_start, double f_step, int f_size,
		double* psitaus,
		char* outputFile_name,
		window win,
		double* win_top_t,
		double* win_top_f,
		int win_top_n,
		double* win_bot_t,
		double* win_bot_f,
		int win_bot_n);

	// �������� �������������� �������-��������������, �������� ������� ���������� � ����� transformedSignal_fileName ��� ��������
	// �������: �� t_start �� t_start+t_step*t_size � ����� t_step;
	// ������: �� f_start �� f_start+f_step*f_size � ����� f_step.
	// ����������� ������� ����� �������� waveletFunction. �������� �������� �� ����������� ����������� ��������! ��������� � �������� make_waveletFunction_equalStep.
	// ������� ��������� �������������� ������������ � ���� originalSignal_fileName.
	// ������� f ������� � ��������� s, ������������ � ������������ �������� ������������ �������-��������������, ������������ f=1/s.
	// ��������! ���������� ������� ����� ����� ������ ������� ������� ���������, ��� ������. ��� ����� ���� ������� � ��������� ��������������� �����������.
	// � ����� ������ ������������� ���� ����������� ������������ ������ ����������� �������, ���� ��������� ���������� ��������� ������� � �������������� (��������� t_sizing).
	extern "C" __declspec(dllexport) void backWavelet(char* transformedSignal_fileName,
		double t_start, double t_step, int t_size,
		double f_start, double f_step, int f_size,
		double* waveletFunction,
		char* originalSignal_fileName);
}
/*
int main(void){
	int sample_size = 4E3;
	double sample_step = 100E-9; //0 to min
	int sample_perSec_max = 125E6; //0 to default (125E6)
	double voltage_max = 5;
	int reading_number = 100;
	int trigger_level = 0; //from -5 to 5, where 5 equals voltage_max
	double* result_arr = new double[sample_size];
	string result_fileName = "result.txt"; //"" to skip writing
	//#pragma omp parallel for num_threads(THREADS_NUM)
	for (int i = 0; i < sample_size; i++)
		result_arr[i] = 0;

	ViStatus status = OWON::getOWONData2(sample_size, sample_step, sample_perSec_max, voltage_max, trigger_level, reading_number, result_arr, result_fileName);
	if (status == VI_SUCCESS) {
		OWON::doLog("Success!\n");
		printf("Success!\n");
	}
}
*/