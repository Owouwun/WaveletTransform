#pragma once
#ifndef FURIE_TRANS
#define FURIE_TRANS

#include <C:\Users\user1\source\repos\owon\owon\common.h>

using namespace std;

//Регистровое кэширование, развёртка циклов
complex<float>* furie3pp(float* arr, float omega0, float deltaOmega, float noOmega, float t0, float deltaT, int noT)
{
	complex<float>* res = new complex<float>[noOmega];

	complex<float>* B = new complex<float>[noT];
	complex<float> temp;
	float omega = omega0;
	for (int k = 0; k < noOmega; k++, omega += deltaOmega) {
		B[1] = exp(-1if* omega* deltaT);
		temp = arr[0] * 0.5f;
		for (int p = 1; p < noT; p++) {
			temp += B[p] * arr[p];
			B[p + 1] = B[p] * B[1];
		}
		res[k] = deltaT / sqrt2pi * temp * exp(-1if* t0* omega);
	}
	return res;
}
float* inversionFurie3(complex<float>* arr, float omega0, float deltaOmega, int noOmega, float t0, float deltaT, int noT)
{
	float* res = new float[noT];

	complex<float>* B = new complex<float>[noOmega];
	complex<float> temp;
	float t = t0;
	for (int k = 0; k < noT; k++, t += deltaT) {
		B[1] = exp(1if* t* deltaOmega);
		temp = arr[0] * 0.5f;
		for (int p = 1; p < noOmega; p++) {
			temp += B[p] * arr[p];
			B[p + 1] = B[p] * B[1];
		}
		temp *= exp(1if* omega0* t);
		res[k] = 2 * (temp * exp(1if* omega0* t)).real()* deltaOmega / sqrt2pi;
	}
	return res;
}
float sinc(float x) {
	float xx = (x * x) / (PI * PI);
	return (1 - xx) * (1 - xx / 4);
}
complex<float>* triangleFurie(float a, float omega0, float deltaOmega, float noOmega)
{ //a - полуширина
	complex<float>* res = new complex<float>[noOmega];
	float oa2;
	float sn;
	for (int p = 1; p < noOmega; p++) {
		oa2 = (omega0 + p * deltaOmega) * a / 2;
		sn = sin(oa2);
		if (abs(oa2) > 0.1)
			res[p] = a * sn * sn / (oa2 * oa2);
		else
			res[p] = a * sinc(oa2);
		res[p] /= sqrt2pi;
	}
	res[0] = res[1];
	return res;
}
float* triangle3(float a0, float a, float t0, float deltaT, int noT)
{
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

#endif