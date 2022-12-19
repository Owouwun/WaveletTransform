#pragma once
#ifndef FURIE_TRANS
#define FURIE_TRANS

#include <C:\Users\user1\source\repos\PicoscopeCpp\PicoscopeCpp\common.h>

using namespace std;

complex<float>* furie3pp(float* arr, float omega0, float deltaOmega, float noOmega, float t0, float deltaT, int noT);
float* inversionFurie3(complex<float>* arr, float omega0, float deltaOmega, int noOmega, float t0, float deltaT, int noT);
float sinc(float x);
complex<float>* triangleFurie(float a, float omega0, float deltaOmega, float noOmega);
float* triangle3(float a0, float a, float t0, float deltaT, int noT);

#endif