#ifndef CUBE_H_
#define CUBE_H_

#include "wavelet.h"

#define SEISMIC_CUBE	3
#define CORR_CUBE	2
#define DEFAULT_CUBE	1
#define SIM_CUBE	0

extern struct layers_t;

typedef struct type_Cube
{
	long x_max;
	long y_max;
	long z_max;

	unsigned long size;

	int type;

	char *name;

	float* values;
	
	double* sum_x;
	double* sum_x2;
	double* r;
	double 	GlobalSum_x;
	double  GlobalSum_x2;

	int layers_num;
	double* layersCorr;
	double 	globalCorr;
	double 	corrAvg;
} tCube;

extern void Cube(long x, long y, long z, int layers_num, int type, tCube** this);
extern void _Cube(tCube* this);

extern void Cube_freeGrid(tCube*);
extern void Cube_allocGrid(tCube*);

extern void Cube_associateFile(char*, tCube*);

extern int Cube_ReadFromFile_name(char* file, tCube* this);
extern int Cube_ReadFromFileFast(char* file, tCube* this);

extern int Cube_SaveToFile_name(char* file, int root, tCube* this);
extern int Cube_SaveToFile(tCube*);
	
extern	void Cube_Init0(tCube* this);
extern 	void Cube_Init(double, tCube* this);

extern 	void Cube_SetCell(int x, int y, int z, double value, tCube* this);
extern 	double Cube_GetCell(int x, int y, int z, tCube* this);

extern 	void Cube_DeleteValues(tCube* this);

extern 	tCube * Cube_MakeReflectionCoefsCube(tCube* this);

extern 	tCube * Cube_MakeSyntheticCube(tWavelet *, tCube* this);

extern 	void Cube_LoadNewLayers(unsigned short int, tCube* this );

extern 	void Cube_GenerateLayersCorr(struct layers_t*, tCube*, tCube* this);

extern 	void Cube_CalculateCorrData(struct layers_t*, tCube* this);

extern 	void Cube_CalculateStaticCorrSums(tCube* this);

extern 	double Cube_GetLayerCorrAverage(tCube* this);

extern 	double Cube_GetGlobalCorrAverage(tCube* this);

extern 	double* Cube_copyLayersCorr(tCube* this);

extern 	char* Cube_CorrAverageDump(tCube* this);

extern 	void Cube_DumpCRs(tCube* this);

extern 	long Cube_Coords(int x, int y, int z, tCube* this);


#endif 

