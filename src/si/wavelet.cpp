#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "wavelet.h"

Wavelet::Wavelet(Configs* config)
{
	this->wavelet_used_values = config->GetParamInt("WAVELET", "USED_VALUES");

	this->points = (int*)malloc((this->wavelet_used_values +1) * sizeof(int));
	this->values = (float*)malloc((this->wavelet_used_values+1) * sizeof(float));
}

Wavelet::~Wavelet(void)
{
	delete[] this->points;
	delete[] this->values;
}

int Wavelet::ReadFromFile(char* file)
{
	FILE *fp;

	pdebug2("Opening file for wavelet input.");

	// open file for reading
	fp = fopen(file, "r");
	if (fp == NULL)
	{
		perror("Error reading wavelet file");
		return -1;
	}

	pdebug2("Reading wavelet values.");

	int i = 0;

	// read values from file until EOF
	while (fscanf(fp, "%d %f\n", &(this->points[i]), &(this->values[i])) != EOF)
	{
		//printf("Wavelet[%d] = Pair(%d, %f)\n", i, this->points[i], this->values[i]);
		i++;
	}
	
	// all needed wavelet values read ?
	if (i<this->wavelet_used_values)
	{
		fprintf(stderr, "EOF found in %s. Uncomplete wavelet.\nThis wavelet must have %d values.\n", file, this->wavelet_used_values);
		return -2;
	}

	// close file descriptor
	fclose(fp);

	return 0;
}

double Wavelet::PointValue(int point)
{
	//printf("Point %d is at array position %d.\n", point, (WAVELET_USED_VALUES / 2) + point);
	return this->values[(this->wavelet_used_values / 2) + point];
}

double Wavelet::IndexValue(int index)
{
	//printf("Point %d is at array position %d.\n", point, (WAVELET_USED_VALUES / 2) + point);
	return this->values[index];
}


