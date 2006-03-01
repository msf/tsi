#ifndef _TSI_H
#define _TSI_H


#include "configs.h"
#include "cube.h"
#include "wavelet.h"
#include "sectionsdefs.h"
#include "timer.h"
#include "log.h"
#include "dssconst.h"
#include "tsi-utils.h"

#define VERSION	"v4.2-pre"


typedef struct type_TSI
{
	// data relative to threading of TSI
	unsigned int ThreadsNumber;

	
	tConfigs *configs;

	/* data for the DSS */
	float DssParameters[DSSDLL_TOTAL_PARAMS_NUM];
	float *DssModels;

	/* HardData aka Wells */
	double *HardData;
	//float *HardData;
	int HardDataSize;

	unsigned long x_values;
	unsigned long y_values;
	unsigned long z_values;

	/* info about layers (aka sections) */
	struct layers_t	*layers;




	tWavelet *wavelet;

	tCube *seismic_cube;
	tCube *mask;
	
	// timers
	tTimer timer_global;		// total time
	
	unsigned int TotalSimulationNumber;
	unsigned int TotalIterationNumber;


	tLog* log;	

	char* filename_ai_files_dir;
	char* filename_ai_filename;
	char* filename_ai_sufix;

	char* filename_synth_files_dump_dir;
	char* filename_synth_dump_filename;
	char* filename_synth_dump_sufix;

	char* filename_corr_files_dir;
	char* filename_corr_dump_filename;
	char* filename_corr_dump_sufix;

	char* filename_bcm_fullname;
	char* filename_bai_fullname;

	char* filename_global_log;

	int dumpSynth;
	int dumpCorr;
	int dumpBCM;
	int dumpBAI;

	double global_null_value;

	int reserv_section;
	int root;
} tTSI;

	/* Methods */ 
extern void TSI(tTSI** this);
extern void _TSI(tTSI* this);
extern int TSI_Initialization(char *, tTSI* this);

extern int TSI_Dss(tCube* outputCube, tCube* DssBAI, tCube* DssBCM, int kType, int iter, int simul, tTSI* this);

extern int TSI_Si(tCube* aiData, tCube** corrCube, int iterNum, int simulNum, tTSI* this);

extern int TSI_Compare (tCube* aiData, tCube* corrCube,  struct aiGridData *, tTSI* this);

extern int TSI_Dump(tCube* SiBAI, tCube** SiBCM, struct aiGridData *,
			int iterNum, int bestGlobalSimulationNumber, tTSI* this);

extern int TSI_Shutdown(tTSI* this);

extern void TSI_ReadDSSParameters(tTSI* this);
extern int TSI_ReadHardDataFile(char* dataFilename, tTSI* this);

#endif 



