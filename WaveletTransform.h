#pragma once
#ifndef WAVELET_TRANS
#define WAVELET_TRANS

#include <C:\Users\user1\source\repos\owon\owon\common.h>
#include <ctime>
#include <omp.h>

using namespace std;
namespace wavelet_trans
{
	// Количество потоков, выделяемых на параллельное вычисление
	int threads_number = 8;
	const double h3cst = 1. / (2. / 15. * sqrt(30.) * pow(PI, -0.25));

	// Список поддерживаемых вейвлет-функций
	enum wavelets
	{
		hermitian1 = 0, // Самый простой из описываемых вейвлетов
		hermitian2 = 1, // Соответствует вейвлету "мексиканская шляпа".
		hermitian3 = 2, // Показал себя лучше всех во время тестов.
		hermitian4 = 3,
		poisson2 = 4, // Показал себя хуже всех во время тестов.
		morlet = 5,
		modified_morlet = 6
	};

	// Численное интегрирование значений функции f с шагом step размера size.
	// Метод интегрирования: Составная формула Симпсона (формула Котеса)
	// Внимание! Если значение интеграла выйдет за пределы допустимого диапазона (будет равен inf), вернётся 0.
	complex<double> kotes(complex<double>* f, double step, int size) {
		complex<double> res = f[0] + f[size - 1];

		for (int k = 1; k < size - 2; k += 2) {
			res += 2. * f[k];
			res += 4. * f[k + 1];
		}

		if (res == res)
			return res * step * 0.33333333333333333;
		else
			return 0;
	}

	// Список доступных вейвлет-преобразований.
	// В преобразованиях применяется кратный сдвиг по значению аргумента. Кратность сдвига определялась эмперическим путём, выполняя условие соответствия конечного преобразования частотной оси.
	// Без наличия сдвига итоговый график мог "растягиваться" или "сжиматься" относительно оси частот.
	// При выборе используемого вейвлета рекомендуется применять эрмитов вейвлет (3) в виду наилучшей разрешенности по времени. То есть ему нужно меньше всего временных семплов для схождения значений.
	// (Количество временных семплов важнее частотных, т.к. время работы алгоритма зависит от количества временных семплов квадратично, в то время как от частотных - лишь линейно)
	complex<double> morlet_wavelet(double t) {
		const double sigma = 2 * PI;
		const double kappa = exp(-0.5 * sigma * sigma);
		return (exp(1.0i * sigma * t) - kappa) * exp(-0.5 * t * t);
	}
	double modified_morlet_wavelet(double t) {
		const double sigma = 2 * PI;
		return cos(sigma * t) / cosh(t);
	}
	double hermitian1_wavelet(double t) {
		t *= 4;
		return t * exp(-t * t * 0.5);
	}
	double hermitian2_wavelet(double t) {
		t *= 4;
		return (1 - t * t) * exp(-t * t * 0.5);
	}
	double hermitian3_wavelet(double t) {
		t *= 4;
		return h3cst*(t * t * t - 3 * t) * exp(-t * t * 0.5);
	}
	double hermitian4_wavelet(double t) {
		t *= 4;
		return (t * t * t * t - 6 * t * t + 3) * exp(-t * t * 0.5);
	}
	double poisson2_wavelet(double t) {
		t *= 2;
		if (t < 0)
			return 0;
		else
			return (t - 2) * t * exp(-t) * 0.5;
	}

	void sizing_toDouble(int t_input_startIndex, int t_input_endIndex, double t_input_step, int t_input_size, int t_sizing, double &t_sizing_start, double &t_sizing_end, double &t_sizing_step, int &t_sizing_size) {
		t_sizing_start = t_input_startIndex * t_input_step;
		t_sizing_end = t_input_endIndex * t_input_step;
		t_sizing_step = t_input_step * t_sizing;
		t_sizing_size = (t_input_endIndex - t_input_startIndex) / t_sizing;
		return;
	}

	// Определение вейвлет-функции (материнского вейвлета) res вида wavelet для времени от t_start-t_step*t_size до t_start+t_step*t_size с шагом t_step и частот от f_start до f_start+f_step*f_size c шагом f_step.
	// Частота f связана с масштабом s, используемом в оригинальных формулах непрерывного вейвлет-преобразования, соотношением f=1/s.
	// Внимание! В значение вейвлет-функции включен множитель 1/sqrt(f). В оригинальных формулах функций он не представлен, однако так как он всегда используется в прямом и обратном вейвлет-преобразовании, он был включен в саму вейвлет-функцию.
	// Внимание! Материнский вейвлет "отражен" относительно t_start по времени с целью оптимищзации вычислений: вместо того, чтобы иметь тройной массив типа [f_size x t_size x t_size], достаточно иметь массив [f_size x 2*t_size], т.к. tau (параллельный перенос) имеет тот же шаг по значениям, что и время.
	void make_waveletFunction_equalStep(double t_start, double t_step, int t_size, double f_start, double f_step, int f_size, wavelets wavelet, complex<double>** res) {
		switch (wavelet) {
		hermitian1:
			for (int i = 0; i < f_size; i++) {
				double divsqrt = 1 / sqrt(f_start + i * f_step);
				for (int j = 0; j < 2 * t_size; j++) {
					double t_taus = (j - t_size) * t_step * (f_start + i * f_step);
					res[i][j] = hermitian1_wavelet(t_taus) * divsqrt;
				}
			}
			break;
		case hermitian2:
			for (int i = 0; i < f_size; i++) {
				double divsqrt = 1 / sqrt(f_start + i * f_step);
				for (int j = 0; j < 2 * t_size; j++) {
					double t_taus = (j - t_size) * t_step * (f_start + i * f_step);
					res[i][j] = hermitian2_wavelet(t_taus) * divsqrt;
				}
			}
			break;
		case hermitian3:
			for (int i = 0; i < f_size; i++) {
				double divsqrt = 1 / sqrt(f_start + i * f_step);
				for (int j = 0; j < 2 * t_size; j++) {
					double t_taus = (j - t_size) * t_step * (f_start + i * f_step);
					res[i][j] = hermitian3_wavelet(t_taus);
					res[i][j] *= divsqrt;
				}
			}
			break;
		case hermitian4:
			for (int i = 0; i < f_size; i++) {
				double divsqrt = 1 / sqrt(f_start + i * f_step);
				for (int j = 0; j < 2 * t_size; j++) {
					double t_taus = (j - t_size) * t_step * (f_start + i * f_step);
					res[i][j] = hermitian4_wavelet(t_taus);
					res[i][j] *= divsqrt;
				}
			}
			break;
		case poisson2:
			for (int i = 0; i < f_size; i++) {
				double divsqrt = 1 / sqrt(f_start + i * f_step);
				for (int j = 0; j < 2 * t_size; j++) {
					double t_taus = (j - t_size) * t_step * (f_start + i * f_step);
					res[i][j] = poisson2_wavelet(t_taus) * divsqrt;
				}
			}
			break;
		case morlet:
			for (int i = 0; i < f_size; i++) {
				double divsqrt = 1 / sqrt(f_start + i * f_step);
				for (int j = 0; j < 2 * t_size; j++) {
					double t_taus = (j - t_size) * t_step * (f_start + i * f_step);
					res[i][j] = morlet_wavelet(t_taus) * divsqrt;
				}
			}
			break;
		case modified_morlet:
			for (int i = 0; i < f_size; i++) {
				double divsqrt = 1 / sqrt(f_start + i * f_step);
				for (int j = 0; j < 2 * t_size; j++) {
					double t_taus = (j - t_size) * t_step * (f_start + i * f_step);
					res[i][j] = modified_morlet_wavelet(t_taus) * divsqrt;
				}
			}
			break;
		}
	}

	// Вывод матрицы matrix типа type размерности [size1 x size2] в файл под названием outputFileName
	template <class type>
	void printMatrix(type** matrix, int size1, int size2, string outputFileName) {
		ofstream fout(outputFileName);
		for (int i = 0; i < size1; i++) {
			for (int j = 0; j < size2; j++) {
				fout << matrix[i][j] << " ";
				//fout << abs(matrix[i][j]) << " ";
				//fout << log(abs(matrix[i][j])) << " ";
			}
			fout << endl;
		}
		fout << endl;
		fout.close();
	}

	// Вывод массива arr типа type размера size в файл под названием outputFileName
	template <class type>
	void printArray(type* arr, int size, string outputFileName) {
		ofstream fout(outputFileName);
		for (int i = 0; i < size; i++)
			fout << arr[i] << endl;
		fout.close();
	}

	// Прямое вейвлет-преобразование функции, значения которой взяты из файла inputFileName, имеют временной шаг t_input_step и количество значений t_input_size;
	// В преобразовании используется каждое t_sizing значение, начиная с t_start_index и заканчивая t_end_index;
	// На выходе получается матрица размерности [f_size x t_end_index-t_start_index)/sizing], соответствующие
	// частотам с f_start по f_start+f_step*f_size с шагом f_step
	// и временам с t_input_step*t_start_index по t_input_step*t_end_index с шагом t_input_step*t_sizing.
	// Материнский вейвлет имеет значения psitaus.
	// Значения преобразования записываются в файл outputFile_name.
	// Частота f связана с масштабом s, используемом в оригинальных формулах непрерывного вейвлет-преобразования, соотношением f=1/s.
	// Обратите внимание на особенности массива, хранящего материнский вейвлет! Подробнее в описании функции построения материнского вейвлета make_waveletFunction_equalStep.
	void waveletTransform_fromIndexToIndex(string inputFIleName,
		double t_input_step, double t_input_size,
		int t_start_index, int t_end_index, int t_sizing,
		double f_start, double f_step, int f_size,
		complex<double>** psitaus,
		string outputFile_name) {
		double t_start = t_start_index * t_input_step;
		double t_step = t_input_step * t_sizing;
		int t_size = (t_end_index - t_start_index) / t_sizing;

		double** res = new double* [f_size];
		for (int i = 0; i < f_size; i++)
			res[i] = new double[t_size];

		if (t_input_step < 0) {
			cout << "Incorrect input time step!" << endl;
			return;
		}
		if (t_input_size <= 0) {
			cout << "Incorrect input time size!" << endl;
			return;
		}
		if (t_start_index < 0 || t_start_index > t_input_size) {
			cout << "Incorrect time start!" << endl;
			return;
		}
		if (t_end_index <= t_start_index || t_end_index > t_input_size) {
			cout << "Incorrect time end!" << endl;
			return;
		}
		if (t_sizing < 1 || t_sizing > t_input_size) {
			cout << "Incorrect time sizing!" << endl;
			return;
		}
		if (f_start < 0) {
			cout << "Incorrect frequency start!" << endl;
			return;
		}
		if (f_step < 0) {
			cout << "Incorrect frequency step!" << endl;
			return;
		}
		if (f_size <= 0) {
			cout << "Incorrect frequency size!" << endl;
			return;
		}

		double* X = new double[t_size];
		ifstream fin(inputFIleName);
		double gbg;
		for (int i = 0; i < t_start_index; i++)
			fin >> gbg;
		int read_size = t_end_index - t_start_index;
		double divsizing = 1. / t_sizing;
		for (int i = 0; i < read_size; i++) {
			if (i % t_sizing == 0)
				fin >> X[int(i * divsizing)];
			else fin >> gbg;
		}
		fin.close();

		int t_current_start_index = t_start / t_step; //Номер начального семпла после sizing
		int t_current_end_index = t_current_start_index + t_size; //Номер последнего семпла после sizing

		int nt = threads_number;
		omp_set_dynamic(0);
		omp_set_num_threads(nt);
#pragma omp parallel
		{
			int tn = omp_get_thread_num();
			complex<double>* f = new complex<double>[t_size];
			for (int i = tn; i < f_size; i += nt)
				for (int j = t_current_start_index; j < t_current_end_index; j++) {
					int temp = t_size - j; // Объяснения о данном сдвиге содержатся в описании функции построения материнского вейвлета make_waveletFunction_equalStep.
					for (int k = t_current_start_index; k < t_current_end_index; k++)
						f[k - t_current_start_index] = X[k - t_current_start_index] * psitaus[i][temp + k];
					res[i][j - t_current_start_index] = kotes(f, t_step, t_size).real();
				}
			delete []f;
		}

		delete []X;

		printMatrix(res, f_size, t_size, outputFile_name);
	}

	// Функция, аналогичная waveletTransform_fromIndexToIndex, но время определяется не начальным и конечным индексами, а значениями в секундах.
	void waveletTransform_fromTimeToTime(string inputFIleName,
		double t_input_step, double t_input_size,
		double t_start, double t_end, int t_sizing,
		double f_start, double f_step, int f_size,
		complex<double>** psitaus,
		string outputFile_name) {
		waveletTransform_fromIndexToIndex(inputFIleName, t_input_step, t_input_size, t_start / t_input_step, t_end / t_input_step, t_sizing, f_start, f_step, f_size, psitaus, outputFile_name);
	}

	// Обратное преобразование вейвлет-преобразования, значения которго определены в файле transformedSignal_fileName для значений
	// времени: от t_start до t_start+t_step*t_size с шагом t_step;
	// частот: от f_start до f_start+f_step*f_size с шагом f_step.
	// Материнский вейвлет имеет значения waveletFunction. Обратите внимание на особенности определения вейвлета! Подробнее в описании make_waveletFunction_equalStep.
	// Функция обратного преобразования записывается в файл originalSignal_fileName.
	// Частота f связана с масштабом s, используемом в оригинальных формулах непрерывного вейвлет-преобразования, соотношением f=1/s.
	// Внимание! Полученная функция может иметь кратно большую частоту колебаний, чем должна. Это может быть связано с принципоп неопределённости Гейзенберга.
	// В таком случае рекомендуется либо попробовать использовать другой материнский вейвлет, либо увеличить количество временных семплов в преобразовании (уменьшить t_sizing).
	void backWavelet(string transformedSignal_fileName,
		double t_start, double t_step, int t_size,
		double f_start, double f_step, int f_size,
		complex<double>** waveletFunction,
		string originalSignal_fileName) {

		double** transformedSignal = new double* [f_size];
		for (int i = 0; i < f_size; i++)
			transformedSignal[i] = new double[t_size];
		ifstream fin(transformedSignal_fileName);
		for (int i = 0; i < f_size; i++)
			for (int j = 0; j < t_size; j++)
				fin >> transformedSignal[i][j];
		fin.close();

		double* res = new double[t_size];

		double* axis_f = new double[f_size];
		for (int i = 0; i < f_size; i++)
			axis_f[i] = f_start + i * f_step;

		int nt = threads_number;
		omp_set_dynamic(0);
		omp_set_num_threads(nt);
#pragma omp parallel
		{
			complex<double>* f = new complex<double>[t_size];
			complex<double>* g = new complex<double>[f_size];
			double sqrf;
			int tn = omp_get_thread_num();
			for (int k = tn; k < t_size; k += nt) {
				for (int i = 0; i < f_size; i++) {
					sqrf = axis_f[i] * axis_f[i];
					for (int j = 0; j < t_size; j++)
						f[j] = transformedSignal[i][j] * waveletFunction[i][t_size + k - j] * sqrf;
					g[i] = kotes(f, t_step, t_size);
				}
				res[k] = kotes(g, f_step, f_size).real();
			}
			delete[] f, g;
		}

		delete[] axis_f;

		printArray(res, t_size, originalSignal_fileName);
	}
}

#endif

/*
* Using example
*/

/*
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
*/