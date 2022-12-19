//Работает только для x64, т.к. библиотеки имеют именно такую разрядность

#include <stdio.h>
#include <iostream>
#include "windows.h"
#include <conio.h>
#include "C:\Program Files\Pico Technology\SDK\inc\ps5000aApi.h"
#include <fstream>
#include <complex>
#include <chrono>
#include <thread>
#define PREF4 __stdcall

#define Sleep(a) usleep(1000*a)
#define scanf_s scanf
#define fscanf_s fscanf
#define memcpy_s(a,b,c,d) memcpy(a,c,d)

#define MAX_CHANNELS 4

#define PI 3.1415926535897932384626433832795

#pragma comment(lib, "C:\\Program Files\\Pico Technology\\SDK\\lib\\ps5000a.lib")

const float sqrt2pi = sqrt(2 * PI);

using namespace std;

typedef struct
{
	int16_t DCcoupled;
	PS5000A_RANGE range;
	int16_t enabled;
} CHANNEL_SETTINGS;

typedef enum
{
	MODEL_NONE = 0, MODEL_PS5203 = 5203, MODEL_PS5204 = 5204
} MODEL_TYPE;

typedef struct
{
	int16_t handle;
	MODEL_TYPE model;
	PS5000A_RANGE firstRange;
	PS5000A_RANGE lastRange;
	unsigned char signalGenerator;
	unsigned char external;
	int16_t ChannelCount;
	CHANNEL_SETTINGS channelSettings[MAX_CHANNELS];
	PS5000A_RANGE triggerRange;
} UNIT_MODEL;

const char* getPicoStatus(PICO_STATUS status) {
	if (status == PICO_OK)
		return "OK.";
	if (status == PICO_MAX_UNITS_OPENED)
		return "MAX UNITS OPENED.";
	if (status == PICO_MEMORY_FAIL)
		return "MEMORY FAIL.";
	if (status == PICO_NOT_FOUND)
		return "NOT FOUND.";
	if (status == PICO_NO_SAMPLES_AVAILABLE)
		return "NO SAMPLES AVAILABLE.";
	if (status == PICO_SIG_GEN_PARAM)
		return "SIG GEN PARAM.";
	if (status == PICO_SIGGEN_PK_TO_PK)
		return "SIGGEN PK TO PK.";
	if (status == PICO_INVALID_NUMBER_CHANNELS_FOR_RESOLUTION)
		return "INVALID NUMBER CHANNELS FOR RESOLUTION.";
	return "UNRECOGNIZED STATUS!";
}
PICO_STATUS optGetValues(int noOfGettings, PS5000A_RANGE range, int32_t sampleCountAfterTrigger, int16_t* maxBuffersA, int16_t* minBuffersA, float* arr, uint32_t timebase = 4) {
	UNIT_MODEL unit;
	PICO_STATUS status = ps5000aOpenUnit(
		&unit.handle,
		NULL,
		PS5000A_DEVICE_RESOLUTION::PS5000A_DR_16BIT
	);
	if (status) return status;

	int32_t nMaxSamples;
	status = ps5000aMemorySegments(
		unit.handle,
		1,
		&nMaxSamples
	);
	if (status) return status;

	status = ps5000aSetChannel(
		unit.handle,
		PS5000A_CHANNEL::PS5000A_CHANNEL_A,
		1,
		PS5000A_COUPLING::PS5000A_AC,
		range,
		0
	);
	if (status) return status;

	status = ps5000aSetSimpleTrigger(
		unit.handle,
		1,
		PS5000A_CHANNEL::PS5000A_EXTERNAL,
		20000,
		PS5000A_THRESHOLD_DIRECTION::PS5000A_RISING,
		0,
		22222
	);
	if (status) return status;

	int32_t sampleCountBeforeTrigger = sampleCountAfterTrigger / 10;
	int32_t noOfAllSamples = sampleCountAfterTrigger + sampleCountBeforeTrigger;

	ps5000aSetDataBuffers(
		unit.handle,
		PS5000A_CHANNEL::PS5000A_CHANNEL_A,
		maxBuffersA,
		minBuffersA,
		noOfAllSamples,
		0,
		PS5000A_RATIO_MODE::PS5000A_RATIO_MODE_NONE
	);
	if (status) return status;

	for (int32_t i = 0; i < noOfAllSamples; i++)
		arr[i] = 0;

	int32_t timeIndisposed;
	for (int k = 0; k < noOfGettings; k++) {
		status = ps5000aRunBlock(
			unit.handle,
			sampleCountBeforeTrigger,
			sampleCountAfterTrigger,
			timebase,
			&timeIndisposed,
			0,
			NULL,
			IntToPtr
		);
		if (status) return status;

		int16_t ready = 0;
		while (!ready)
			ps5000aIsReady(unit.handle, &ready);

		uint32_t noOfSamples = noOfAllSamples;
		int16_t overflow;
		status = ps5000aGetValues(
			unit.handle,
			0,
			&noOfSamples,
			1,
			enPS5000ARatioMode::PS5000A_RATIO_MODE_NONE,
			0,
			&overflow
		);
		if (status) return status;
	
		for (int32_t i = 0; i < noOfAllSamples; i++)
			arr[i] += maxBuffersA[i] + minBuffersA[i];
	}
	if (status) return status;

	ps5000aStop(unit.handle);

	ps5000aCloseUnit(unit.handle);

	double sum = 0;
	for (uint32_t k = 0; k < noOfAllSamples; k++)
		sum += arr[k];

	for (uint32_t k = 0; k < sampleCountBeforeTrigger + sampleCountAfterTrigger/100; k++)
		arr[k] = 0;
	
	const float avg = sum / float(noOfAllSamples);
	for (uint32_t k = sampleCountBeforeTrigger + sampleCountAfterTrigger / 100; k < noOfAllSamples; k++)
		arr[k] -= avg;

	const float aaa = 65535.0 * 2000.0;
	for (uint32_t k = sampleCountBeforeTrigger + sampleCountAfterTrigger / 100; k < noOfAllSamples; k++) {
		arr[k] = arr[k] / aaa;
	}

	return status;
}

float* getValues(int noOfGettings) {
	UNIT_MODEL unit;
	PICO_STATUS status = ps5000aOpenUnit(
		&unit.handle,
		NULL,
		PS5000A_DEVICE_RESOLUTION::PS5000A_DR_16BIT
	);
	printf("Handle: %d\n", unit.handle);

	if (status != PICO_OK)
	{
		printf("Unable to open device\n");
		printf("Error code : %d\n", status);
		while (!_kbhit())
			;
		exit(99); // exit program - nothing after this executes
	}
	printf("Device opened successfully\n\n");

	int32_t nMaxSamples;
	status = ps5000aMemorySegments(
		unit.handle,
		1,
		&nMaxSamples
	);
	cout << "Memory segments status:" << endl << getPicoStatus(status) << endl;

	status = ps5000aSetChannel(
		unit.handle,
		PS5000A_CHANNEL::PS5000A_CHANNEL_A,
		1,
		PS5000A_COUPLING::PS5000A_AC,
		PS5000A_RANGE::PS5000A_2V,
		0
	);
	cout << "Set channel status:" << endl << getPicoStatus(status) << endl;

	status = ps5000aSetSimpleTrigger(
		unit.handle,
		1,
		PS5000A_CHANNEL::PS5000A_EXTERNAL,
		20000,
		PS5000A_THRESHOLD_DIRECTION::PS5000A_RISING,
		0,
		22222
	);
	cout << "Set simple trigger status:" << endl << getPicoStatus(status) << endl;

	int32_t sampleCountAfterTrigger = 50000;
	int32_t sampleCountBeforeTrigger = sampleCountAfterTrigger / 10;
	int32_t noOfAllSamples = sampleCountAfterTrigger + sampleCountBeforeTrigger;

	int16_t *minBuffersA = new int16_t[noOfAllSamples];
	int16_t *maxBuffersA = new int16_t[noOfAllSamples];
	ps5000aSetDataBuffers(
		unit.handle,
		PS5000A_CHANNEL::PS5000A_CHANNEL_A,
		maxBuffersA,
		minBuffersA,
		noOfAllSamples,
		0,
		PS5000A_RATIO_MODE::PS5000A_RATIO_MODE_NONE
	);
	cout << "Set data buffers status:" << endl << getPicoStatus(status) << endl;

	float *maxA = new float[noOfAllSamples];
	for (int32_t i = 0; i < noOfAllSamples; i++)
		maxA[i] = 0;

	cout << "Getting values..." << endl;

	int32_t timeIndisposed;
	bool retry;
	uint32_t timebase = 4;
	auto start = clock();
	for (int k = 0; k < noOfGettings; k++) {
		do {
			retry = false;
			status = ps5000aRunBlock(
				unit.handle,
				sampleCountBeforeTrigger,
				sampleCountAfterTrigger,
				timebase,
				&timeIndisposed,
				0,
				NULL,
				IntToPtr
			);
			if (status != PICO_OK) {
				cout << "Run block status:" << endl << getPicoStatus(status) << endl;
				return NULL;
			}
			if (status == PICO_POWER_SUPPLY_CONNECTED || status == PICO_POWER_SUPPLY_NOT_CONNECTED || status == PICO_POWER_SUPPLY_UNDERVOLTAGE) {
				status = ps5000aChangePowerSource(unit.handle, status);
				retry = true;
			}
		} while (retry);

		if (status == PICO_OK) {
			int16_t ready = 0;
			while (!ready) {
				
				ps5000aIsReady(unit.handle, &ready);
			}
		}
		uint32_t noOfSamples = (uint32_t)noOfAllSamples;
		int16_t overflow;
		status = ps5000aGetValues(
			unit.handle,
			0,
			&noOfSamples,
			1,
			enPS5000ARatioMode::PS5000A_RATIO_MODE_NONE,
			0,
			&overflow
		);
		if (status != PICO_OK) {
			cout << "Get values status:" << endl << getPicoStatus(status) << endl;
			return NULL;
		}
		for (int32_t i = 0; i < noOfAllSamples; i++) {
			maxA[i] += maxBuffersA[i] + minBuffersA[i];
		}
		
	}
	auto finish = clock();
	cout << "time :" << float(finish - start)*1000.0 / float(noOfGettings) << endl;
	cout << "Get values status:" << endl << getPicoStatus(status) << endl;

	ps5000aStop(unit.handle);

	ps5000aCloseUnit(unit.handle);

	double sr = 0;
	for (uint32_t k = 0; k < noOfAllSamples; k++)
		sr += maxA[k];
	sr /= float(noOfAllSamples);
	for (uint32_t k = 0; k < noOfAllSamples; k++)
		maxA[k] -= sr;

	float aaa = 65535.0 * 2000.0;
	for (uint32_t k = 0; k < noOfAllSamples; k++) {
		maxA[k] = float(maxA[k] / aaa);
	}
	return maxA;
}

void printValues(float* maxA, uint32_t arraySize, const char* fileName) {
	ofstream fout(fileName);
	for (uint32_t i = 0; i < arraySize; i++)
		fout << maxA[i] << endl;
	fout.close();
}
void printValues(double* maxA, uint32_t arraySize, const char* fileName) {
	ofstream fout(fileName);
	for (uint32_t i = 0; i < arraySize; i++)
		fout << maxA[i] << endl;
	fout.close();
}
void printValues(complex<double>* maxA, uint32_t arraySize, const char* fileName) {
	double real, imag;
	ofstream fout(fileName);
	for (uint32_t i = 0; i < arraySize; i++) {
		real = maxA[i].real();
		imag = maxA[i].imag();
		fout << real << "\t" << imag << "\t" << sqrt(real*real+imag*imag) << endl;
	}
	fout.close();
}
void printValues(complex<float>* maxA, uint32_t arraySize, const char* fileName) {
	double real, imag;
	ofstream fout(fileName);
	for (uint32_t i = 0; i < arraySize; i++) {
		real = maxA[i].real();
		imag = maxA[i].imag();
		fout << real << "\t" << imag << "\t" << sqrt(real*real + imag * imag) << endl;
	}
	fout.close();
}

void printFrequences(float omega0, float deltaOmega, float noOmega, const char* fileName) {
	ofstream fout(fileName);
	for (uint32_t i = 0; i < noOmega; i++) {
		fout << (omega0 + i * deltaOmega) / (2 * PI) << endl;
	}
	fout.close();
}

void extSigGen() {
	UNIT_MODEL unit;
	PICO_STATUS status = ps5000aOpenUnit(
		&unit.handle,
		NULL,
		PS5000A_DEVICE_RESOLUTION::PS5000A_DR_8BIT
	);
	printf("Handle: %d\n", unit.handle);

	if (status != PICO_OK)
	{
		printf("Unable to open device\n");
		printf("Error code : %d\n", status);
		while (!_kbhit())
			;
		exit(99); // exit program - nothing after this executes
	}
	printf("Device opened successfully\n\n");

	status = ps5000aSetSigGenBuiltIn(
		unit.handle,
		0, //const
		500000, //Напряжение не выходе (выше - громче, но не больше 2В (2000000 = +-1V output))
		enPS5000AWaveType::PS5000A_SQUARE,
		3000, //Начальная частота генерируемой волны (меньше 3000Гц не слышно для синусоиды)
		10000, //Конечная частота генерируемой волны
		1000, //Инкремент частоты генерируемой волны
		10, //Продолжительность задержки?
		PS5000A_SWEEP_TYPE::PS5000A_UP,
		PS5000A_EXTRA_OPERATIONS::PS5000A_ES_OFF, //const
		1, //Количество birst?-волн
		0, //Количество свип-волн
		enPS5000ASigGenTrigType::PS5000A_SIGGEN_GATE_HIGH,
		enPS5000ASigGenTrigSource::PS5000A_SIGGEN_EXT_IN, //const
		6553 //6553=1V?
	);
	cout << getPicoStatus(status) << endl;

	ps5000aCloseUnit(unit.handle);
}

void softSigGen() {
	UNIT_MODEL unit;
	PICO_STATUS status = ps5000aOpenUnit(
		&unit.handle,
		NULL,
		PS5000A_DEVICE_RESOLUTION::PS5000A_DR_8BIT
	);
	printf("Handle: %d\n", unit.handle);

	if (status != PICO_OK)
	{
		printf("Unable to open device\n");
		printf("Error code : %d\n", status);
		while (!_kbhit())
			;
		exit(99); // exit program - nothing after this executes
	}
	printf("Device opened successfully\n\n");

	status = ps5000aSetSigGenBuiltIn(
		unit.handle,
		0, //const
		2000000, //const (+-1V output)
		enPS5000AWaveType::PS5000A_SQUARE, //const
		500, //Минимум 3000 для синусоидальных
		500,
		0,
		0,
		PS5000A_SWEEP_TYPE::PS5000A_UP,
		PS5000A_EXTRA_OPERATIONS::PS5000A_ES_OFF, //const
		1,
		0,
		enPS5000ASigGenTrigType::PS5000A_SIGGEN_GATE_HIGH,
		enPS5000ASigGenTrigSource::PS5000A_SIGGEN_SOFT_TRIG,
		0
	);
	cout << getPicoStatus(status) << endl;
	status = ps5000aSigGenSoftwareControl(unit.handle, 1);
	cout << getPicoStatus(status) << endl;
	for (int i = 0; i < 1000000000; i++)
		int a = 0;
	status = ps5000aSigGenSoftwareControl(unit.handle, 0);
	cout << getPicoStatus(status) << endl;

	ps5000aCloseUnit(unit.handle);
}
//Регистровое кэширование, развёртка циклов
complex<float>* furie3pp(float* arr, float omega0, float deltaOmega, float noOmega, float t0, float deltaT, int noT) {
	complex<float>* res = new complex<float>[noOmega];

	complex<float>* B = new complex<float>[noT];
	complex<float> temp;
	float omega = omega0;
	for (int k = 0; k < noOmega; k++, omega+=deltaOmega) {
		B[1] = exp(-1if * omega * deltaT);
		temp = arr[0] * 0.5f;
		for (int p = 1; p < noT; p++) {
			temp += B[p] * arr[p];
			B[p + 1] = B[p] * B[1];
		}
		res[k] = deltaT / sqrt2pi * temp * exp(-1if * t0 * omega);
	}

	return res;
}

float* inversionFurie3(complex<float>* arr, float omega0, float deltaOmega, int noOmega, float t0, float deltaT, int noT) {
	float* res = new float[noT];


	complex<float>* B = new complex<float>[noOmega];
	complex<float> temp;
	float t = t0;
	for (int k = 0; k < noT; k++, t += deltaT) {
		B[1] = exp(1if * t * deltaOmega);
		temp = arr[0] * 0.5f;
		for (int p = 1; p < noOmega; p++) {
			temp += B[p] * arr[p];
			B[p + 1] = B[p] * B[1];
		}
		temp *= exp(1if * omega0 * t);
		res[k] = 2 * (temp * exp(1if * omega0 * t)).real() * deltaOmega / sqrt2pi;
	}

	return res;
}

float sinc_norm(float x) {
	float xx = x * x;
	return (1 - xx) * (1 - xx / 4);
}
float sinc(float x) {
	float xx = (x * x) / (PI*PI);
	return (1 - xx) * (1 - xx / 4);
}

complex<float>* triangleFurie(float a, float omega0, float deltaOmega, float noOmega) { //a - полуширина
	complex<float>* res = new complex<float>[noOmega];
	float oa2;
	float sn;
	for (int p = 1; p < noOmega; p++) {
		oa2 = (omega0 + p * deltaOmega) * a / 2;
		sn = sin(oa2);
		if (abs(oa2)>0.1)
			res[p] = a * sn*sn / (oa2*oa2);
		else
			res[p] = a * sinc(oa2);
		res[p] /= sqrt2pi;
	}
	res[0] = res[1];

	return res;
}

float* triangle3(float a0, float a, float t0, float deltaT, int noT) {
	float* tr = new float[noT];
	int i = 0;
	const float s1 = (a0 - a - t0) / deltaT;
	const float s2 = (a0 - t0) / deltaT;
	const float s3 = (a0 + a - t0) / deltaT;
	for (; i < s1; i++)
		if (i < noT)
			tr[i] = 0;
		else
			return tr;
	for (; i < s3; i++)
		if (i < noT)
			tr[i] = 1 - abs(t0 + deltaT * i - a0) / a;
		else
			return tr;
	for (; i < noT; i++)
		if (i < noT)
			tr[i] = 0;
		else
			return tr;
	return tr;
}

float* sin(float t0, float deltaT, float noT) {
	float* res = new float[noT];
	for (int i = 0; i < noT; i++)
		res[i] = sin(t0 + deltaT * i);
	return res;
}

float transformtoabs1(float a, float b, float x) {
	return -1 + 2 * (x-a) / (b-a);
}
float transformtoab(float x, float a, float b) {
	float c = (a + b) / 2;
	return c + x * (b-a) / 2;
}

float Kotes(float* Arr, int size, float a0, float h, float a, float b) {
	int ia = (a - a0) / h; // a = a0 + ia * h
	int ib = (b - a0) / h-1; // b = a0 + ib * h
	int simple = 1;
	
	float res = Arr[ia] + Arr[ib];
	int i;
	for (i = ia+1; i < ib-1; i += 2*simple) {
		res += 4 * Arr[i];
		res += 2 * Arr[i+simple];
	}
	if (i == ib - simple) {
		res += 4 * Arr[i];
	}
	
	return res * h / 3;
}

float MHATphi(float t) {
	return (1 - t * t) * exp(-t * t / 2);
}

float underintegral(float* arr, float a, float h, float size, float tau, float s) {
	s = 1 / s;
	const float sqrts = sqrt(s);

	float* Arr = new float[size];
	for (int i = 0; i < size; i++)
		Arr[i] = arr[i] / sqrts * MHATphi((a + i * h - tau) / s); // x(t) / sqrt(s) * phi(t - tau / s)
	return Kotes(Arr, size, a, h, a, a + h * size);
}
float underintegral2(float* arr, float a, float h, float size, float tau, float s) {
	s = 1 / s;
	const float sqrts = sqrt(s);

	const float Ets0 = exp(-(pow((a - tau) / s, 2) / 2));
	const float econst = exp(-(2 * (a - h)*(h - tau) + pow(h - tau, 2)) / (s*s));
	const float e1sq = exp(-2*(h*(h - tau) / (s*s)));

	float temp = Ets0;
	
	float* Arr = new float[size];
	Arr[0] = arr[0] / sqrts * temp;
	for (int i = 1; i < size; i++) {
		temp *= econst * e1sq;
		Arr[i] = arr[i] / sqrts * temp; // x(t) / sqrt(s) * phi(t - tau / s)
	}
	return Kotes(Arr, size, a, h, a, a + h * size);
}

float* testSignal(float a, float h, float size) {
	float* res = new float[size];
	const float frq1 = 512, amp1 = 1;
	const float frq2 = 64, amp2 = 2;
	float t;

	for (int i = 0; i < size * 0.1; i++)
		res[i] = 0;
	for (int i = size * 0.1; i < size * 0.3; i++) {
		t = a + i * h;
		res[i] = amp1 * sin(2*PI * frq1 * t);
	}
	for (int i = size * 0.3; i < size * 0.6; i++)
		res[i] = 0;
	for (int i = size * 0.6; i < size * 0.9; i++) {
		t = a + i * h;
		res[i] = amp2 * sin(2*PI * frq2 * t);
	}
	for (int i = size * 0.9; i < size * 0.1; i++)
		res[i] = 0;

	return res;
}

const double sigma = 2 * PI;
const double kappa = exp(-0.5 * sigma * sigma);
const double c = pow(1 + exp(-sigma * sigma) - 2 * exp(-0.75 * sigma*sigma), -0.5) / pow(PI, -0.25);
void morlet(float a, float h, float size, float* res) {
	//const double ec = exp(-0.5*(2 * h*(a - h) + h * h));
	const double emhh = exp(-h * h);
	const complex<double> esh = exp(1i*sigma*double(h));

	double mainExp1;
	double tempExp1;
	complex<double> mainExp2;

	mainExp1 = exp(-0.5 * a * a);
	//tempExp1 = ec * emhh;
	tempExp1 = exp(-0.5*(2 * h*(a - h) + h * h) - h * h);
	mainExp2 = exp(1i*sigma*double(a));
	res[0] = ((mainExp2 - kappa) * c * mainExp1).real();
	for (int i = 1; i < size; i++) {
		mainExp1 *= tempExp1;
		tempExp1 *= emhh;
		mainExp2 *= esh;
		res[i] = ((mainExp2 - kappa) * c * mainExp1).real();
		//cout << exp(1i*sigma*double(a+i*h)) << " " << mainExp2 << endl;
	}
}
void morlet_opt(float a, float h, float size, float* res) {
	//const double ec = exp(-0.5*(2 * h*(a - h) + h * h));
	const double emhh = exp(-h * h);
	const complex<double> esh = exp(1i*sigma*double(h));

	double mainExp1;
	double tempExp1;
	complex<double> mainExp2;

	mainExp1 = exp(-0.5 * a * a) * c;
	//tempExp1 = ec * emhh;
	tempExp1 = exp(-h * (a + 0.5 * h));
	mainExp2 = exp(1i*sigma*double(a));
	res[0] = ((mainExp2 - kappa) * mainExp1).real();
	for (int i = 1; i < size; i++) {
		mainExp1 *= tempExp1;
		tempExp1 *= emhh;
		mainExp2 *= esh;
		res[i] = ((mainExp2 - kappa) * mainExp1).real();
		//cout << exp(1i*sigma*double(a+i*h)) << " " << mainExp2 << endl;
	}
}

float wavelet(float* x, float a, float h, int size, float tau, float s) {
	s = 1 / s;

	float *psitaus = new float[size];
	//psitaus = morlet((a-tau)/s, h/s, size);
	
	const float sqrts = sqrt(s);

	float* f = new float[size];
	for (int i = 0; i < size; i++) {
		f[i] = x[i] / sqrts * psitaus[i]; // x(t) / sqrt(s) * phi(t - tau / s)
	}
	return Kotes(f, size, a, h, a, a + h * size);
}

float wavelet_opt(float* x, float a, float h, int size, float tau, float s, float* psitaus, float* f ) {
	morlet_opt((a - tau) * s, h * s, size, psitaus);

	const float sqrts = sqrt(s);

	for (int i = 0; i < size; i++)
		f[i] = x[i] * sqrts * psitaus[i]; // x(t) / sqrt(s) * phi(t - tau / s)

	return Kotes(f, size, a, h, a, a + h * size);
}
float wavelet_opt2(float* x, float a, float h, int size, float tau, float s, float* psitaus, float* f) {
	//morlet start
	const float m_a = (a - tau) * s;
	const float m_h = h * s;

	const double emhh = exp(-m_h * m_h);
	const complex<double> esh = exp(1i*sigma*double(m_h));

	double mainExp1;
	double tempExp1;
	complex<double> mainExp2;

	mainExp1 = exp(-0.5 * m_a * m_a) * c;
	tempExp1 = exp(-m_h * (m_a + 0.5 * m_h));
	mainExp2 = exp(1i*sigma*double(m_a));
	psitaus[0] = ((mainExp2 - kappa) * mainExp1).real();
	for (int i = 1; i < size; i++) {
		mainExp1 *= tempExp1;
		tempExp1 *= emhh;
		mainExp2 *= esh;
		psitaus[i] = ((mainExp2 - kappa) * mainExp1).real();
	}
	//morlet end

	const float sqrts = sqrt(s);

	for (int i = 0; i < size; i++)
		f[i] = x[i] * sqrts * psitaus[i];

	//Kotes start
	float res = f[0] + f[size - 1];
	int i;
	for (i = 1; i < size - 2; i += 2) {
		res += 4 * f[i];
		res += 2 * f[i + 1];
	}
	if (i == size - 2) {
		res += 4 * f[i];
	}

	return res * h * 0.33333333f;
	//Kotes end
}
float wavelet_opt3(float* x, float a, float h, int size, float tau, float s, float* psitaus, float* f, double emhh, complex<double> esh) {
	//morlet start
	const float m_a = (a - tau) * s;
	const float m_h = h * s;

	long double mainExp1;
	long double tempExp1;
	complex<double> mainExp2;

	mainExp1 = exp(-0.5 * m_a * m_a) * c;
	tempExp1 = exp(-m_h * (m_a + 0.5 * m_h));
	mainExp2 = exp(1i*sigma*double(m_a));
	psitaus[0] = ((mainExp2 - kappa) * double(mainExp1)).real();
	for (int i = 1; i < size; i++) {
		mainExp1 *= tempExp1;
		tempExp1 *= emhh;
		mainExp2 *= esh;
		psitaus[i] = ((mainExp2 - kappa) * double(mainExp1)).real();
	}
	//morlet end

	const float sqrts = sqrt(s);

	for (int i = 0; i < size; i++)
		f[i] = x[i] * sqrts * psitaus[i];

	//Kotes start
	float res = f[0] + f[size - 1];
	int i;
	for (i = 1; i < size - 2; i += 2) {
		res += 4 * f[i];
		res += 2 * f[i + 1];
	}
	if (i == size - 2) {
		res += 4 * f[i];
	}

	return res * h * 0.33333333f;
	//Kotes end
}


//tau - time, 1/s - freq
int main() {
	/*
	float a = 0;
	float h = 0.01;
	float size = 300;
	float* s = new float[size];
	s = sin(a, h, size);
	//cout << Kotes(s, size, a, h, a + size * h/4, a + size*h/2);
	*/

	/*
	int s1 = 1000, s2 = 1000;
	float** res = new float*[s1];
	for (int i = 0; i < s1; i++)
		res[i] = new float[s2];

	const float a = 0;
	const float h = PI * 0.1;
	const float size = 1000;
	float* X = new float[size];
	X = testSignal(a, h, size);
	*/
	//for (int i=0; i<s1; i++)
	//	for (int j=0; j<s2; j++)
	//		res[i][j] = wavelet(X, a, h, size, a + i*h, j);
	/*
	ofstream fout("output1.txt");
	for (int i = 0; i < s1; i++) {
		for (int j = 0; j < s2; j++)
			fout << res[i][j] << " ";
		fout << endl;
	}
	fout.close();	
	*/
	
	/*
	float a = 0; float h = 0.000094; int t_size = 1000;
	float* X = new float[t_size];
	X = testSignal(a, h, t_size);
	*/
	float a = 0; float h = 0.000094; int t_size = 1000;
	float* X = new float[t_size];
	ifstream fin("input1.txt");
	for (int i = 0; i < t_size; i++) {
		if (!i%100)
			fin >> X[i];
	}
	fin.close();

	int s_size = 500; int tau_size = 1000;
	float s_step = 4 * 2 / PI;
	float** res = new float*[s_size];
	float* axisTau = new float[tau_size];
	float* axisS = new float[s_size];

	float *psitaus = new float[t_size];
	float* f = new float[t_size];

	for (int i = 0; i < tau_size; i++)
		axisTau[i] = a + i * h;
	for (int i = 0; i < s_size; i++)
		axisS[i] = s_step * i;

	//void* ptr = new float*[s_size*tau_size];
	for (int i = 0; i < s_size; i++)
		res[i] = new float[tau_size];
		//res[i] =(float*) ( (long long int)ptr+ i*tau_size) ;

	float s_size_div_10 = s_size * 0.1;
	float m_h; double emhh; complex<double> esh;
	for (int i = 0; i < s_size_div_10; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "10% is ready" << endl;
	for (int i = s_size_div_10; i < 2*s_size_div_10; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "20% is ready" << endl;
	for (int i = 2 * s_size_div_10; i < 3*s_size_div_10; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "30% is ready" << endl;
	for (int i = 3*s_size_div_10; i < 4*s_size_div_10; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "40% is ready" << endl;
	for (int i = 4*s_size_div_10; i < 5*s_size_div_10; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "50% is ready" << endl;
	for (int i = 5*s_size_div_10; i < 6*s_size_div_10; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "60% is ready" << endl;
	for (int i = 6*s_size_div_10; i < 7*s_size_div_10; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "70% is ready" << endl;
	for (int i = 7*s_size_div_10; i < 8*s_size_div_10; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "80% is ready" << endl;
	for (int i = 8*s_size_div_10; i < 9*s_size_div_10; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "90% is ready" << endl;
	for (int i = 9*s_size_div_10; i < s_size; i++) {
		m_h = h * s_step * i;
		emhh = exp(-m_h * m_h);
		esh = exp(1i*sigma*double(m_h));
		for (int j = 0; j < tau_size; j++)
			res[i][j] = wavelet_opt3(X, a, h, t_size, a + j * h, s_step * i, psitaus, f, emhh, esh);
	}
	cout << "100% is ready" << endl;

	delete[] psitaus;
	delete[] f;
	ofstream f1out("output1.txt");
	for (int i = 0; i < s_size; i++) {
		for (int j = 0; j < tau_size; j++)
			f1out << res[i][j] << " ";
		f1out << endl;
	}
	f1out.close();
	ofstream f2out("output2.txt");
	for (int i = 0; i < tau_size; i++) {
		f2out << i << '\t' << axisTau[i] << endl;
	}
	f2out.close();
	ofstream f3out("output3.txt");
	for (int i = 0; i < s_size; i++) {
		f3out << i << '\t' << axisS[i] << endl;
	}
	f3out.close();
	

	/*
	const float a = -10;
	const float h = 0.01;
	const int size = 2000;
	float *mor = new float[size];
	mor = morlet(a, h, size);
	printValues(mor, size, "output1.txt");
	*/

	/*
	//Переменные, которые нужно определить перед началом сбора данных.
	//Последний описанный метод - ps5000aSetDataBuffers.
	uint32_t sampleCountAfterTrigger = 50000;
	uint32_t sampleCountBeforeTrigger = sampleCountAfterTrigger / 10;
	uint32_t noOfAllSamples = sampleCountBeforeTrigger + sampleCountAfterTrigger;
	float timebase = 4;
	float* arr = new float[noOfAllSamples];
	int16_t *minBuffersA = new int16_t[noOfAllSamples];
	int16_t *maxBuffersA = new int16_t[noOfAllSamples];

	//arr = getValues(2000);
	//PICO_STATUS status = optGetValues(100, PS5000A_RANGE::PS5000A_2V, sampleCountAfterTrigger, maxBuffersA, minBuffersA, arr);
	//cout << status;

	//printValues(arr, noOfAllSamples, "output.txt");
	
	float t0 = 0;
	float deltaT = 16E-9 * (timebase - 3);
	float T = deltaT * noOfAllSamples;
	float deltaOmega = 2 * PI / T / 10;
	float omega0 = 0;
	float noOmega = 10000;
	*/

	//Отсечение помех при помощи Фурье-преобразования
	/*
	complex<float>* fArr = new complex<float>[noOmega];
	fArr = furie3pp(arr, omega0, deltaOmega, noOmega, t0, deltaT, noOfAllSamples);
	printValues(fArr, noOmega, "output1.txt");
	arr = inversionFurie3(fArr, omega0, deltaOmega, noOmega, t0, deltaT, noOfAllSamples);
	printFrequences(omega0, deltaOmega, noOmega, "omegas.txt");
	printValues(arr, noOfAllSamples, "output2.txt");
	*/
	
	//Быстрый и точный треугольный Фурье-сигнал
	/*
	float* tr = new float[noOfAllSamples];
	complex<float>* fTr = new complex<float>[noOmega];
	fTr = triangleFurie(T/2, omega0, deltaOmega, noOmega);
	//printValues(fTr, noOmega, "output.txt");
	tr = inversionFurie3(fTr, omega0, deltaOmega, noOmega, t0-T/2, deltaT, 2*noOfAllSamples);
	printValues(tr, noOfAllSamples , "output.txt");
	*/
	
	//Настраиваемый треугольный Фурье-сигнал
	/*
	float* tr = new float[noOfAllSamples];
	tr = triangle3(T/2, T/2, t0, deltaT, noOfAllSamples);
	printValues(tr, noOfAllSamples, "output1.txt");
	complex<float>* fTr = new complex<float>[noOmega];
	fTr = furie3pp(tr, omega0, deltaOmega, noOmega, t0, deltaT, noOfAllSamples);
	printValues(fTr, noOmega, "output2.txt");
	tr = inversionFurie3(fTr, omega0, deltaOmega, noOmega, t0, deltaT, noOfAllSamples);
	printValues(tr, noOfAllSamples, "output3.txt");
	*/

	//sin-сигнал
	/*
	float t0 = 0; float deltaT = 0.000094; int noOfAllSamples = 1000;
	float omega0 = 0; float deltaOmega = 10 * 2 / PI; int noOmega = 2000;
	float* s = new float[noOfAllSamples];
	s = testSignal(t0, deltaT, noOfAllSamples);
	printValues(s, noOfAllSamples, "output1.txt");
	complex<float>* fS = new complex<float>[noOmega];
	fS = furie3pp(s, omega0, deltaOmega, noOmega, t0, deltaT, noOfAllSamples);
	printValues(fS, noOmega, "output2.txt");
	s = inversionFurie3(fS, omega0, deltaOmega, noOmega, t0, deltaT, noOfAllSamples);
	printValues(s, noOfAllSamples, "output3.txt");
	*/

	//Прямоугольный сигнал
	/*
	float* r = new float[noOfAllSamples];
	for (int i = 0; i < noOfAllSamples; i++)
		r[i] = int(i > noOfAllSamples / 3 && i < noOfAllSamples * 2 / 3);
	printValues(r, noOfAllSamples, "output1.txt");
	complex<float>* fR = new complex<float>[noOmega];
	fR = furie3pp(r, omega0, deltaOmega, noOmega, t0, deltaT, noOfAllSamples);
	r = inversionFurie3(fR, omega0, deltaOmega, noOmega, t0, deltaT, noOfAllSamples);
	printValues(r, noOfAllSamples, "output2.txt");
	*/

	/*
	complex<float>* tririe = triangleFurie(50, omega0, deltaOmega/16, noOmega);
	//printValues(tririe, noOmega, "output.txt");
	float* tr = new float[4*noOfAllSamples];
	tr = inversionFurie3(tririe, omega0, deltaOmega/16, noOmega, -2*deltaT*noOfAllSamples, deltaT, 4*noOfAllSamples);
	printValues(tr, 4*noOfAllSamples, "output.txt");
	*/

	//for (int k = 0; k < noOfAllSamples; k++)
	//	test[k] = sin(2 * PI * 5 * f0 * k*deltaT);
	//printValues(arr, noOfAllSamples, "output.txt");
	//complex<float>* fArr = furie3(arr, noOfAllSamples, omega0, deltaOmega, noOmega, 0, deltaT);
	//printValues(fArr, noOmega, "output.txt");
	//printFrequences(omega0, deltaOmega, noOmega, "omegas.txt");
	//arr = inversionFurie3(tririe, omega0, deltaOmega, noOmega, t0, deltaT, noOfAllSamples);
	//printValues(arr, noOfAllSamples, "output.txt");

	//printFrequences(omega0, deltaOmega, noOmega, "omegas.txt");
	//complex<float>* F = furie3(arr, noOfAllSamples, omega0, deltaOmega, noOmega, 0, deltaT);
	//printValues(F, noOmega, "output.txt");
	/*double* Fabs = new double[500];
	for (int k = 0; k < 500; k++)
		Fabs[k] = abs(F[k]);
	printValues(Fabs, 500, "output.txt");*/
	//complex<double>* furie = DFT(fMaxA, noOfAllSamples);
	//ofstream fout("output.txt");
	//for (uint32_t i = 0; i < noOfAllSamples; i++)
	//	fout << abs(furie[i]) << endl;
	//fout.close();
	//softSigGen();
	//extSigGen();
}
