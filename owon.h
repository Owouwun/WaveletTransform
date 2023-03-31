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
#include <chrono>

#include "C:\Program Files\IVI Foundation\VISA\Win64\Include\visa.h"

#pragma comment(lib, "C:/Program Files/IVI Foundation/VISA/Win64/Lib_x64/msc/visa64.lib")

using namespace std;
using namespace std::chrono;

extern "C" {
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

	#define THREADS_NUM 1

		//Режим работы канала
		enum COUPLING {
			COUPLING_AC,
			COUPLING_DC
		};

		//Вертикальное масштабирование
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

		//Режим сбора данных
		/*enum ACQUIRE_MODE {
			SAMPLE,
			PEAK
		};*/

		//Разрядность
		/*enum PRECISION {
			P8,
			P12,
			P14
		};*/

		//OFFSet (смещение по вертикальной оси):
												/*
												no. of division for vertical offset of the
												displayed signal from the specified
												channel
												*/
												//RANGe: offset 0-10M, size 1-256K (вывод данных, начиная с offset-го считанного семпла, в количестве size семплов

		void vertScale(VERTICAL_SCALE sc, float& scCoef, const char*& chSc) {
			switch (sc) {
			case VS2mV:
				chSc = "2mV";
				scCoef = 0.002;
				return;
			case VS5mV:
				chSc = "5mV";
				scCoef = 0.005;
				return;
			case VS10mV:
				chSc = "10mV";
				scCoef = 0.01;
				return;
			case VS20mV:
				chSc = "20mV";
				scCoef = 0.02;
				return;
			case VS50mV:
				chSc = "50mV";
				scCoef = 0.05;
				return;
			case VS100mV:
				chSc = "100mV";
				scCoef = 0.1;
				return;
			case VS200mV:
				chSc = "200mV";
				scCoef = 0.2;
				return;
			case VS500mV:
				chSc = "500mV";
				scCoef = 0.5;
				return;
			case VS1V:
				chSc = "1V";
				scCoef = 1;
				return;
			case VS2V:
				chSc = "2V";
				scCoef = 2;
				return;
			case VS5V:
				chSc = "5V";
				scCoef = 5;
				return;
			default:
				return;
			}
		}

		string vertScale_to_string(VERTICAL_SCALE sc) {
			switch (sc) {
			case VS2mV:
				return "2mV";
			case VS5mV:
				return "5mV";
			case VS10mV:
				return "10mV";
			case VS20mV:
				return "20mV";
			case VS50mV:
				return "50mV";
			case VS100mV:
				return "100mV";
			case VS200mV:
				return "200mV";
			case VS500mV:
				return "500mV";
			case VS1V:
				return "1V";
			case VS2V:
				return "2V";
			case VS5V:
				return "5V";
			default:
				return "";
			}
		}

		string horiScale_to_string(HORIZONTAL_SCALE hs) {
			switch (hs) {
				//case HS1ns:return"1.0ns";
				//case HS2ns:return"2.0ns";
			//case HS5ns:return"5.0ns";
			case HS10ns:return"10ns";
			case HS20ns:return"20ns";
			case HS50ns:return"50ns";
			case HS100ns:return"100ns";
			case HS200ns:return"200ns";
			case HS500ns:return"500ns";
			case HS1us:return"1.0us";
			case HS2us:return"2.0us";
			case HS5us:return"5.0us";
			case HS10us:return"10us";
			case HS20us:return"20us";
			case HS50us:return"50us";
			case HS100us:return"100us";
			case HS200us:return"200us";
			case HS500us:return"500us";
			case HS1ms:return"1.0ms";
			case HS2ms:return"2.0ms";
			case HS5ms:return"5.0ms";
			case HS10ms:return"10ms";
			case HS20ms:return"20ms";
			case HS50ms:return"50ms";
			case HS100ms:return"100ms";
			case HS200ms:return"200ms";
			case HS500ms:return"500ms";
			case HS1s:return"1.0s";
			case HS2s:return"2.0s";
			case HS5s:return"5.0s";
			case HS10s:return"10s";
			case HS20s:return"20s";
			case HS50s:return"50s";
			case HS100s:return"100s";
			default:return"";
			}
		}
		double horiScale_to_number(HORIZONTAL_SCALE hs) {
			switch (hs) {
				//case HS1ns:return"1.0ns";
				//case HS2ns:return"2.0ns";
			case HS5ns:return 5E-9;
			case HS10ns:return 10E-9;
			case HS20ns:return 20E-9;
			case HS50ns:return 50E-9;
			case HS100ns:return 100E-9;
			case HS200ns:return 200E-9;
			case HS500ns:return 500E-9;
			case HS1us:return 1E-6;
			case HS2us:return 2E-6;
			case HS5us:return 5E-6;
			case HS10us:return 10E-6;
			case HS20us:return 20E-6;
			case HS50us:return 50E-6;
			case HS100us:return 100E-6;
			case HS200us:return 200E-6;
			case HS500us:return 500E-6;
			case HS1ms:return 1E-3;
			case HS2ms:return 2E-3;
			case HS5ms:return 5E-3;
			case HS10ms:return 10E-3;
			case HS20ms:return 20E-3;
			case HS50ms:return 50E-3;
			case HS100ms:return 100E-3;
			case HS200ms:return 200E-3;
			case HS500ms:return 500E-3;
			case HS1s:return 1;
			case HS2s:return 0;
			case HS5s:return 5;
			case HS10s:return 10;
			case HS20s:return 20;
			case HS50s:return 50;
			case HS100s:return 100;
			default:return 0;
			}
		}
		HORIZONTAL_SCALE number_to_horiScale(double& hs) {
			if (hs <= 5E-9) {
				hs = 5E-9;
				return HS5ns;
			}
			if (hs <= 10E-9) {
				hs = 5E-9;
				return HS10ns;
			}
			if (hs <= 20E-9) {
				hs = 20E-9;
				return HS20ns;
			}
			if (hs <= 50E-9) {
				hs = 50E-9;
				return HS50ns;
			}
			if (hs <= 100E-9) {
				hs = 100E-9;
				return HS100ns;
			}
			if (hs <= 200E-9) {
				hs = 200E-9;
				return HS200ns;
			}
			if (hs <= 500E-9) {
				hs = 200E-9;
				return HS500ns;
			}
			if (hs <= 1E-6) {
				hs = 1E-6;
				return HS1us;
			}
			if (hs <= 2E-6) {
				hs = 2E-6;
				return HS2us;
			}
			if (hs <= 5E-6) {
				hs = 5E-6;
				return HS5us;
			}
			if (hs <= 10E-6) {
				hs = 10E-6;
				return HS10us;
			}
			if (hs <= 20E-6) {
				hs = 20E-6;
				return HS20us;
			}
			if (hs <= 50E-6) {
				hs = 50E-6;
				return HS50us;
			}
			if (hs <= 100E-6) {
				hs = 100E-6;
				return HS100us;
			}
			if (hs <= 200E-6) {
				hs = 200E-6;
				return HS200us;
			}
			if (hs <= 500E-6) {
				hs = 500E-6;
				return HS500us;
			}
			if (hs <= 1E-3) {
				hs = 1E-3;
				return HS1ms;
			}
			if (hs <= 2E-3) {
				hs = 1E-3;
				return HS2ms;
			}
			if (hs <= 5E-3) {
				hs = 5E-3;
				return HS5ms;
			}
			if (hs <= 10E-3) {
				hs = 10E-3;
				return HS10ms;
			}
			if (hs <= 20E-3) {
				hs = 20E-3;
				return HS20ms;
			}
			if (hs <= 50E-3) {
				hs = 50E-3;
				return HS50ms;
			}
			if (hs <= 100E-3) {
				hs = 100E-3;
				return HS100ms;
			}
			if (hs <= 200E-3) {
				hs = 200E-3;
				return HS200ms;
			}
			if (hs <= 500E-3) {
				hs = 500E-3;
				return HS500ms;
			}
			if (hs <= 1) {
				hs = 1;
				return HS1s;
			}
			if (hs <= 2) {
				hs = 2;
				return HS2s;
			}
			if (hs <= 5) {
				hs = 5;
				return HS5s;
			}
			if (hs <= 10) {
				hs = 10;
				return HS10s;
			}
			if (hs <= 20) {
				hs = 20;
				return HS20s;
			}
			if (hs <= 50) {
				hs = 50;
				return HS50s;
			}
			hs = 100;
			return HS100s;
		}

		const char* acqDepmem(AQUIRE_DEPMEM ad) {
			switch (ad) {
			case AD1K:
				return "1K";
			case AD10K:
				return "10K";
			case AD100K:
				return "100K";
			case AD1M:
				return "1M";
			case AD10M:
				return "10M";
			case AD100M:
				return "100M";
			default:
				return "";
			}
		}

		const char* coup(COUPLING c) {
			switch (c) {
			case COUPLING_AC:
				return "AC";
			case COUPLING_DC:
				return "DC";
			}
		}

		VERTICAL_SCALE voltageMax_to_vertScale(double& vm) {
			double vs_double = vm / 5;
			if (vs_double <= 0.002) {
				vm = 0.002 * 5;
				return VS2mV;
			}
			if (vs_double <= 0.005) {
				vm = 0.005 * 5;
				return VS5mV;
			}
			if (vs_double <= 0.01) {
				vm = 0.01 * 5;
				return VS10mV;
			}
			if (vs_double <= 0.02) {
				vm = 0.02 * 5;
				return VS20mV;
			}
			if (vs_double <= 0.05) {
				vm = 0.05 * 5;
				return VS50mV;
			}
			if (vs_double <= 0.1) {
				vm = 0.1 * 5;
				return VS100mV;
			}
			if (vs_double <= 0.2) {
				vm = 0.2 * 5;
				return VS200mV;
			}
			if (vs_double <= 0.5) {
				vm = 0.5 * 5;
				return VS500mV;
			}
			if (vs_double <= 1) {
				vm = 1. * 5;
				return VS1V;
			}
			if (vs_double <= 2) {
				vm = 2. * 5;
				return VS2V;
			}
			vm = 5. * 5;
			return VS5V;
		}

		__declspec(dllexport) int getOWONData(
			VERTICAL_SCALE verticalScale,
			HORIZONTAL_SCALE horizontalScale,
			COUPLING coupling,
			AQUIRE_DEPMEM acquireDepmem,
			float horisontalOffset,
			int waveStart, //0M-10M
			int waveSize, //1-256K
			const char* writeFileName,
			int noOfReading,
			long double* result
		) {
			int i;

			/*
			* First we must call viOpenDefaultRM to get the manager
			* handle.  We will store this handle in defaultRM.
			*/
			status1 = viOpenDefaultRM(&defaultRM);
			if (status1 < VI_SUCCESS)
			{
				printf("Could not open a session to the VISA Resource Manager!\n");
				exit(EXIT_FAILURE);
			}
			/* Find all the USB TMC VISA resources in our system and store the  */
			/* number of resources in the system in numInstrs.                  */
			ViConstString str = "USB?*INSTR";

			status2 = viFindRsrc(defaultRM, str, &findList, &numInstrs, instrResourceString);

			if (status2 < VI_SUCCESS)
			{
				printf("An error occurred while finding resources.\nHit enter to continue.");
				fflush(stdin);
				viClose(defaultRM);
				return status2;
			}

			/*
			 * Now we will open VISA sessions to all USB TMC instruments.
			 * We must use the handle from viOpenDefaultRM and we must
			 * also use a string that indicates which instrument to open.  This
			 * is called the instrument descriptor.  The format for this string
			 * can be found in the function panel by right clicking on the
			 * descriptor parameter. After opening a session to the
			 * device, we will get a handle to the instrument which we
			 * will use in later VISA functions.  The AccessMode and Timeout
			 * parameters in this function are reserved for future
			 * functionality.  These two parameters are given the value VI_NULL.
			 */

			for (i = 0; i < numInstrs; i++)
			{
				if (i > 0)
					viFindNext(findList, instrResourceString);

				status = viOpen(defaultRM, instrResourceString, VI_EXCLUSIVE_LOCK, VI_NULL, &instr);

				if (status < VI_SUCCESS)
				{
					printf("Cannot open a session to the device %d.\n", i + 1);
					continue;
				}

				/*
				 * At this point we now have a session open to the USB TMC instrument.
				 * We will now use the viWrite function to send the device the string "*IDN?\n",
				 * asking for the device's identification.
				 */
				 //strcpy(stringinput, "*IDN?\n");
				 //strcpy(stringinput, "*IDN?;*IDN?;");
				 //strcpy(stringinput, ":UPL:TXT?\r\n");
				 //strcpy(stringinput, ":FUNC:OFFS 0;\r\n");
				 //strcpy(stringinput, ":FUNC:OFFS?;\r\n");
				 // DEPMEM - количество вычисляемых семплов
				 // SCAL - шаг по времени у деления (относительно вычисленных семплов)
				 // RANG size - количество берущихся из DEPMEM семплов
				float vertScCoef;
				const char* vertScCode;
				vertScale(verticalScale, vertScCoef, vertScCode);

				//string horiOffs = to_string(horisontalOffset);
				//horiOffs = horiOffs[0] + horiOffs[1]
				string m = (
					string(":*RST;\n") +
					string(":ACQ:PREC 14;\n") +
					string(":ACQ:MODE PEAK;\n") +
					string(":MEAS:TIM 0.002;\n") +
					string(":CH1:SCAL ") + string(vertScCode) + string(";\n") +
					string(":CH1:BAND 20M;\n") +
					string(":CH1:COUP ") + coup(coupling) + string(";\n") +
					string(":CH1:OFFS 0;\n") +
					string(":CH2:COUP AC;\n") +
					string(":CH2:OFFS 0;\n") +
					string(":HORI:SCAL ") + horiScale_to_string(horizontalScale) + string(";\n") +
					string(":ACQ:DEPMEM ") + acqDepmem(acquireDepmem) + string(";\n") +
					string(":TRIG:SING:MODE EDGE;\n") +
					string(":TRIG:SING:EDGE:SOUR CH2;\n") +
					string(":TRIG:SING:EDGE:COUP AC;\n") +
					string(":TRIG:SING:EDGE:SLOP RISE;\n") +
					string(":TRIG:SING:EDGE:LEV 0;\n") +
					string(":HORI:OFFS ") + to_string(horisontalOffset) + string(";\n") +
					string(":WAV:BEG CH1;\n") +
					string(":WAV:RANG ") + to_string(waveStart) + string(",") + to_string(waveSize) + string(";\n") +
					string(":WAV:FETC?;\n") +
					string(":WAV:END;\n")
					);
				//m = string(":STOP;\n");
				//m = string("*IDN?;\n");

				const char* st = m.c_str();

				strcpy(stringinput, st);

				printf(st);

				m = m + m;

				int lastNum;
				int ii;

				//viClear(instr);

				//viSetAttribute(instr, VI_ATTR_IO_PROT, VI_PROT_USBTMC_VENDOR);

				//viSetBuf(instr, VI_IO_IN_BUF, 2 * waveSize + 16);
				//viSetBuf(instr, VI_IO_OUT_BUF, strlen(stringinput));

				for (int k = 0; k < noOfReading; k++) {
					if (k == 0)
						status = viWrite(instr, (ViBuf)stringinput, (ViUInt32)strlen(stringinput), &writeCount);
					else
						status = viWrite(instr, (ViBuf)":WAV:BEG CH1;\n:WAV:FETC?;\n:WAV:END;\n", (ViUInt32)37, &writeCount);
					if (status < VI_SUCCESS) {
						printf("Error writing to the device %d.\n", i + 1);
						status = viClose(instr);
						status = viClose(defaultRM);
						return 0;
					}
					//Sleep(1000);

					lastNum = 0;
					do {
						status = viRead(instr, (ViPBuf)buffer, 2 * waveSize + 16, &retCount);

						if (status < VI_SUCCESS) {
							printf("Error reading a response from the device %d.\n", i + 1);
							return 0;
						}
						//cout << k << endl;
						/*for (int ii = 0; ii < retCount / 2; ii++, j++)
							// >> 2: 16-PRECISION
							//0.015625=2^(PRECISION-8), PRECISION=14
							//fprintf(resultFile, "%f\n", (float(int(buffer[ii * 2])) + int(char(buffer[ii * 2 + 1]) >> 2) * 0.015625)* temp);
							result[j] += (float(int(buffer[ii * 2])) + int(char(buffer[ii * 2 + 1]) >> 2) * 0.015625) * temp;
						}*/
						//cout << k << endl;
	//#pragma omp parallel num_threads(THREADS_NUM) private(ii)
	//					{
							//for (ii = omp_get_thread_num(); ii < retCount / 2; ii += THREADS_NUM) {
						for (ii = 0; ii < retCount / 2; ii++)
							result[lastNum + ii] += (float(int(buffer[ii * 2])) + int(char(buffer[ii * 2 + 1]) >> 2) * 0.015625) * vertScCoef * 0.04;
						//						}
						//					}
											/*for (ii = omp_get_thread_num(); ii < retCount / 2; ii += THREADS_NUM) {
													result[lastNum + ii] += (float(int(buffer[ii * 2])) + int(char(buffer[ii * 2 + 1]) >> 2) * 0.015625) * temp;
												}*/

												/*int ii;
												short* ts;
												short* end = (short*)&(buffer[retCount]);
												double temp2 = temp * 0.015625;
												omp_set_dynamic(0);
												#pragma omp parallel num_threads(THREADS_NUM) private(ii, ts) shared(end)
												{
													ii = omp_get_thread_num();
													for (ts = (short*)&(buffer[2*omp_get_thread_num()]); ts < end; ts+=THREADS_NUM, ii+=THREADS_NUM)
														result[lastNum + ii] += double(char(*ts))*temp + double(*ts >> 10) * temp2;
												}*/
						lastNum += retCount / 2;
						//cout << retCount;
						//result[lastNum + ii] += (float(int(buffer[ii * 2])) + int(char(buffer[ii * 2 + 1]) >> 2) * 0.015625) * temp;
					} while (status == VI_SUCCESS_MAX_CNT);
				}

				double div = 1. / noOfReading;
				for (i = 0; i < waveSize; i++)
					result[i] *= div;

				double sumsr = 0;
				for (i = 0; i < waveSize; i++)
					sumsr += result[i];
				sumsr /= waveSize;
				for (i = 0; i < waveSize; i++)
					result[i] -= sumsr;
				status = viClose(instr);
			}
			//viClear(instr);

			//status = viClose(defaultRM);

			return 0;
		}

		AQUIRE_DEPMEM size_to_acquire(int& sample_size, int& sample_acquire) {
			if (sample_size <= 1E3) {
				sample_acquire = 1E3;
				return AD1K;
			}
			if (sample_size <= 10E3) {
				sample_acquire = 10E3;
				return AD10K;
			}
			if (sample_size <= 100E3) {
				sample_acquire = 100E3;
				return AD100K;
			}
			if (sample_size <= 1E6) {
				sample_acquire = 1E6;
				return AD1M;
			}
			sample_acquire = 10E6;
			sample_size = 10E6;
			return AD10M;
		}

		void doLog(string msg) {
			FILE* log_file = fopen("OWONLog.txt", "a");
			fprintf(log_file, "%s", msg.c_str());
			fclose(log_file);
		}

		int getOWONData2(int sample_size, double sample_step, int sample_perSec_max, double voltage_max, int trigger_level, int reading_number, double* result_arr, string result_fileName) {
			FILE* log_file = fopen("OWONLog.txt", "w");
			fclose(log_file);

			int i;

			status1 = viOpenDefaultRM(&defaultRM);
			if (status1 < VI_SUCCESS)
			{
				doLog("Could not open a session to the VISA Resource Manager!\n");
				printf("Could not open a session to the VISA Resource Manager!\n");
				return status1;
			}
			/* Find all the USB TMC VISA resources in our system and store the  */
			/* number of resources in the system in numInstrs.                  */
			ViConstString str = "USB?*INSTR";

			status2 = viFindRsrc(defaultRM, str, &findList, &numInstrs, instrResourceString);

			if (status2 < VI_SUCCESS) {
				doLog("An error occurred while finding resources.\n");
				printf("An error occurred while finding resources.\nHit enter to continue.");
				fflush(stdin);
				viClose(defaultRM);
				return status2;
			}

			/*
			 * Now we will open VISA sessions to all USB TMC instruments.
			 * We must use the handle from viOpenDefaultRM and we must
			 * also use a string that indicates which instrument to open.  This
			 * is called the instrument descriptor.  The format for this string
			 * can be found in the function panel by right clicking on the
			 * descriptor parameter. After opening a session to the
			 * device, we will get a handle to the instrument which we
			 * will use in later VISA functions.  The AccessMode and Timeout
			 * parameters in this function are reserved for future
			 * functionality.  These two parameters are given the value VI_NULL.
			 */
			for (i = 0; i < numInstrs; i++)
			{
				if (i > 0)
					viFindNext(findList, instrResourceString);

				status = viOpen(defaultRM, instrResourceString, VI_EXCLUSIVE_LOCK, VI_NULL, &instr);

				if (status < VI_SUCCESS) {
					doLog("Cannot open a session to the device " + to_string(i + 1) + "\n");
					printf("Cannot open a session to the device %d.\n", i + 1);
					return status;
				}

				int sample_acquired;
				if (sample_perSec_max == 0)
					sample_perSec_max = 125E6;
				AQUIRE_DEPMEM acquiredDepmem = size_to_acquire(sample_size, sample_acquired);
				double horizontalScale_double = sample_acquired * sample_step / 20;
				HORIZONTAL_SCALE horiScale = number_to_horiScale(horizontalScale_double);
				if (1 / (horizontalScale_double * 20 / sample_acquired) > sample_perSec_max) {
					horizontalScale_double = 1. / sample_perSec_max / 20 * sample_acquired;
					horiScale = number_to_horiScale(horizontalScale_double);
				}
				sample_step = horizontalScale_double * 20 / sample_acquired;
				VERTICAL_SCALE vertScale = voltageMax_to_vertScale(voltage_max);

				doLog("Samples per second: " + to_string(int(1 / sample_step)) + "Spsec\n");
				doLog("Max voltage: " + to_string(voltage_max) + "V\n");
				doLog("Sample size: " + to_string(sample_size) + "S\n");
				printf("Samples per second: %d Spsec\n", int(1 / sample_step));
				printf("Max voltage: %f V\n", voltage_max);
				printf("Sample size: %d S\n", sample_size);

				buffer = new char[sample_size * 2 + 16];

				/*
				 * At this point we now have a session open to the USB TMC instrument.
				 * We will now use the viWrite function to send the device the string "*IDN?\n",
				 * asking for the device's identification.
				 */
				string m = (
					string(":*RST;\n") +
					string(":ACQ:PREC 14;\n") +
					string(":ACQ:MODE SAMP;\n") +
					string(":MEAS:TIM 0.002;\n") +
					string(":CH1:SCAL ") + vertScale_to_string(vertScale) + string(";\n") +
					string(":CH1:BAND 20M;\n") +
					string(":CH1:COUP AC;\n") +
					string(":CH1:OFFS 0;\n") +
					string(":CH2:COUP AC;\n") +
					string(":CH2:OFFS 0;\n") +
					string(":HORI:SCAL ") + horiScale_to_string(horiScale) + string(";\n") +
					string(":ACQ:DEPMEM ") + acqDepmem(acquiredDepmem) + string(";\n") +
					string(":TRIG:SING:MODE EDGE;\n") +
					string(":TRIG:SING:EDGE:SOUR CH2;\n") +
					string(":TRIG:SING:EDGE:COUP AC;\n") +
					string(":TRIG:SING:EDGE:SLOP RISE;\n") +
					string(":TRIG:SING:EDGE:LEV ") + to_string(trigger_level) + string(";\n") +
					string(":HORI:OFFS 10;\n") +
					string(":WAV:BEG CH1;\n") +
					string(":WAV:RANG 0,") + to_string(sample_size) + string(";\n") +
					string(":WAV:FETC?;\n") +
					string(":WAV:END;\n")
					);

				const char* st = m.c_str();

				strcpy(stringinput, st);

				doLog("Query:\n" + m);
				printf("\nQuery:\n%s", st);

				m = m + m;

				int lastNum;
				int ii;
				for (int k = 0; k < reading_number; k++) {
					if (k == 0)
						status = viWrite(instr, (ViBuf)stringinput, (ViUInt32)strlen(stringinput), &writeCount);
					else
						status = viWrite(instr, (ViBuf)":WAV:BEG CH1;\n:WAV:FETC?;\n:WAV:END;\n", (ViUInt32)37, &writeCount);
					if (status < VI_SUCCESS) {
						doLog("Error writing to the device " + to_string(i + 1) + "\n");
						printf("Error writing to the device %d.\n", i + 1);
						status = viClose(instr);
						if (status < VI_SUCCESS) {
							doLog("Error closing the device.\n");
							printf("Error closing the device.\n");
							return status;
						}
						status = viClose(defaultRM);
						if (status < VI_SUCCESS) {
							doLog("Error closing the session.\n");
							printf("Error closing the session.\n");
							return status;
						}
						return status;
					}

					lastNum = 0;
					do {
						status = viRead(instr, (ViPBuf)buffer, 2 * sample_size + 16, &retCount);

						if (status < VI_SUCCESS) {
							doLog("Error reading a response from the device " + to_string(i + 1) + "\n");
							printf("Error reading a response from the device %d.\n", i + 1);
							return status;
						}

						for (ii = 0; ii < retCount / 2; ii++)
							result_arr[lastNum + ii] += (float(int(buffer[ii * 2])) + int(char(buffer[ii * 2 + 1]) >> 2) * 0.015625) * voltage_max * 0.008; //0.008=1/25 / 5

						lastNum += retCount / 2;
					} while (status == VI_SUCCESS_MAX_CNT);
				}

				double div = 1. / reading_number;
				for (i = 0; i < sample_size; i++)
					result_arr[i] *= div;

				double sumsr = 0;
				for (i = 0; i < sample_size; i++)
					sumsr += result_arr[i];
				sumsr /= sample_size;
				for (i = 0; i < sample_size; i++)
					result_arr[i] -= sumsr;

				status = viClose(instr);
				if (status < VI_SUCCESS) {
					doLog("Error closing the device.\n");
					printf("Error closing the device.\n");
					return status;
				}
			}

			if (result_fileName.c_str() != "") {
				FILE* result_file = fopen(result_fileName.c_str(), "w");

				doLog("Writing data...\n");
				printf("Writing data...\n");
				for (int i = 0; i < sample_size; i++)
					fprintf(result_file, "%f\n", result_arr[i]);
				
				fclose(result_file);
			}
			//status = viClose(defaultRM); //Crush
			if (status < VI_SUCCESS) {
				doLog("Error closing the session.\n");
				printf("Error closing the session.\n");
				return status;
			} else {
				OWON::doLog("Success!\n");
				printf("Success!\n");
			}

			return status;
		}
	}
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