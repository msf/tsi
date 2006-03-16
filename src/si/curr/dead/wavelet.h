#ifndef WAVELET_H_
#define WAVELET_H_

#include "configs.h"

typedef struct type_Wavelet
{
	int* points;
	//int points[WAVELET_MAX_VALUES];

	float* values;
	//double values[WAVELET_MAX_VALUES];

	int wavelet_used_values;

	// max number of wavelet values (for point > 0)
	int max_values;
} tWavelet;

extern void Wavelet(tConfigs* config, tWavelet** this);
extern void _Wavelet(tWavelet* this);

extern int Wavelet_ReadFromFile(char* file, tWavelet* this);
extern double Wavelet_PointValue(int point, tWavelet* this);
extern double Wavelet_IndexValue(int index, tWavelet* this);

#endif // WAVELET_H_
