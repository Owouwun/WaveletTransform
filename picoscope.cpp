#include <C:\Users\user1\source\repos\PicoscopeCpp\PicoscopeCpp\common.h>
#include "C:\Program Files\Pico Technology\SDK\inc\ps5000aApi.h"

#pragma comment(lib, "C:\\Program Files\\Pico Technology\\SDK\\lib\\ps5000a.lib")

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

	for (uint32_t k = 0; k < sampleCountBeforeTrigger + sampleCountAfterTrigger / 100; k++)
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
