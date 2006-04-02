#ifndef WAVELET_H_
#define WAVELET_H_

#include "configs.h"

class Wavelet
{
public:
	
	int* points;
	//int points[WAVELET_MAX_VALUES];

	float* values;
	//double values[WAVELET_MAX_VALUES];

	int wavelet_used_values;

	// max number of wavelet values (for point > 0)
	int max_values;

	Wavelet(Configs* config);
	virtual ~Wavelet(void);

	int ReadFromFile(char* file);
	double PointValue(int point);
	double IndexValue(int index);
};

#endif // WAVELET_H_


