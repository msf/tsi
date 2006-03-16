#ifndef _TSI_H
#define _TSI_H


#include "configs.h"
#include "cube.h"
#include "wavelet.h"
#include "sectionsdefs.h"
#include "timer.h"
#include "log.h"
#include "dssconst.h"

#define VERSION	"v4.2-pre"

class TSI
{
public:
	// data relative to threading of TSI
	unsigned int ThreadsNumber;

	
	Configs *configs;

	/* data for the DSS */
	float DssParameters[DSSDLL_TOTAL_PARAMS_NUM];
	float *DssModels;

	double *HardData;
	//float *HardData;
	int HardDataSize;

	unsigned long x_values;
	unsigned long y_values;
	unsigned long z_values;

	/* info about layers (aka sections) */
	struct layers_t	*layers;




	Wavelet *wavelet;

	Cube *seismic_cube;
	Cube *mask;
	
	// timers
	Timer timer_global;		// total time
	
	unsigned int TotalSimulationNumber;
	unsigned int TotalIterationNumber;


	Log* log;	

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


	/* Methods */ 
	TSI();
	virtual ~TSI();
	int Initialization(char *);

	int Dss(Cube* outputCube, Cube* DssBAI, Cube* DssBCM, int kType, int iter, int simul);

	int Si(Cube* aiData, Cube** corrCube, int iterNum, int simulNum);

	int Compare(Cube* aiData, Cube* corrCube, Cube** SiBAI, Cube** SiBCM, 
				Cube** bestAiCubeGlobal, 
				int* bestGlobalSimulation,
				int iterNum, int simul);

	int Dump(Cube* SiBAI, Cube** SiBCM, Cube* bestAiCubeGlobal,
			int iterNum, int bestGlobalSimulationNumber);

	int Shutdown();

	void ReadDSSParameters();
	int ReadHardDataFile(char* dataFilename);
};

#endif 



