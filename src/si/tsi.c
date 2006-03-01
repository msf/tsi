#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dssconst.h"
#include "const.h"
#include "utils.h"
#include "tsi.h"
#include "tsi-utils.h"
#include "timer.h"

void TSI(tTSI** this)
{
	*this = (tTSI*) malloc(sizeof(tTSI));
}

void _TSI(tTSI* this)
{
	free(this);
}

/*
 * Initialization
 *
 * It will read the parameters file, and set the acording environment variables.
 */
int TSI_Initialization(char *configFile, tTSI* this)
{
	int ret_value;
		
	// read config file
	Configs(&this->configs);

	printf("Configs inited\n");

	if(configFile != NULL)
		ret_value = Configs_ReadConfigFile(configFile, this->configs);
	else
		ret_value = Configs_ReadConfigFile(CONFIGS_FILENAME, this->configs);
	if (ret_value < 0) {
		printf("Error reading config file: %s\n", CONFIGS_FILENAME);
		exit(ret_value);
	}

	printf("Configs read file\n");
	
	// GLOBAL parameters for TSI
	this->TotalIterationNumber = Configs_GetParamInt("GLOBAL", "ITERATIONS", this->configs);
	this->TotalSimulationNumber = Configs_GetParamInt("GLOBAL", "SIMULATIONS", this->configs);
	this->ThreadsNumber = (unsigned int) Configs_GetParamInt("GLOBAL","THREADS", this->configs);

	// read cube dimensions
	this->x_values = Configs_GetParamInt("GRID", "XNUMBER", this->configs);
	this->y_values = Configs_GetParamInt("GRID", "YNUMBER", this->configs);
	this->z_values = Configs_GetParamInt("GRID", "ZNUMBER", this->configs);
	
	// read sections and CORR related parameters.
	// TODO: variables that are used for "correlations" should be named: corr_name
	this->root = Configs_GetParamInt("CORR", "BCM_ROOT", this->configs);
	//
	// init layers definition
	this->layers = (struct layers_t *) malloc( sizeof(struct layers_t) );
	this->layers->minimum_number =  (unsigned short int) Configs_GetParamInt("CORR", "LAYERS_MIN", this->configs);
	this->layers->minimum_size = (unsigned short int) Configs_GetParamInt("CORR", "LAYER_SIZE_MIN", this->configs);
	this->layers->total_size = this->z_values;
	this->layers->number = 0; //sane default.

	// READ DSS RELATED PARAMETERS!
	TSI_ReadDSSParameters(this);
	printf("read dss parms\n");



	// read HardData from the hard data file
	ret_value = TSI_ReadHardDataFile(Configs_GetParam("HARDDATA", "FILENAME", this->configs), this);
	if (ret_value < 0) {
		printf("Error reading Hard Data file: %s\n", Configs_GetParam("HARDDATA", "FILENAME", this->configs));
		return ret_value;
	}
	printf("read hard data\n");

	// Read Wavelet from file to memory, including wavelet parameters
	Wavelet(this->configs, &this->wavelet);
	printf ("init wavelet\n");
	ret_value = Wavelet_ReadFromFile(Configs_GetParam("WAVELET", "FILENAME", this->configs), this->wavelet);
	if (ret_value < 0) {
		printf("Error reading Wavelet file: %s\n", Configs_GetParam("WAVELET", "FILENAME", this->configs));
		return ret_value;
	}
	printf ("read wavelet\n");
	// Read Seismic Cube from file to memory
	Cube(this->x_values, this->y_values, this->z_values, 1, CORR_CUBE, &this->seismic_cube); //the 1 is because we don't know yet how many layers we'll use
	printf ("init cube\n");	
	ret_value = Cube_ReadFromFile_name(Configs_GetParam("SEISMIC", "FILENAME", this->configs), this->seismic_cube);
	if (ret_value < 0) {
		printf("Error reading Seismic Cube file: %s\n", Configs_GetParam("SEISMIC", "FILENAME", this->configs));
		return ret_value;
	}
	printf("read seismic\n");

	// init statics values for the seismic cube (values precomputed for globalCorr computation
	Cube_CalculateStaticCorrSums(this->seismic_cube);
	printf("corr calcs\n");

	this->global_null_value = Configs_GetParamDouble("MASK", "NULL_VALUE", this->configs);
	// read mask from file
	if ((Configs_GetParamInt("MASK", "USE_MASK", this->configs)) == 1) {
		Cube(this->x_values, this->y_values, this->z_values, this->layers->number, DEFAULT_CUBE, &this->mask);
		ret_value = Cube_ReadFromFile_name(Configs_GetParam("MASK", "FILENAME", this->configs), this->mask);
		if (ret_value < 0){
			printf("Error reading Mask file: %s\n", Configs_GetParam("MASK", "FILENAME", this->configs));
			return ret_value;
		}
	} else {
		this->mask = NULL;
	}
	printf("done mask\n");


	this->filename_ai_files_dir = Configs_GetParam("FILENAME", "TEMP_DIR", this->configs);
	this->filename_ai_filename = "aiCube";
	this->filename_ai_sufix = ".out";
	

	// init all working paths (AI cubes, Synth cubes, Corr cubes, BCM, BAI, logs)
	this->filename_synth_files_dump_dir = Configs_GetParam("FILENAME", "SYNTH_DUMP_DIR", this->configs);
	this->filename_synth_dump_filename = Configs_GetParam("FILENAME", "SYNTH_DUMP", this->configs);
	this->filename_synth_dump_sufix = Configs_GetParam("FILENAME", "SYNTH_DUMP_SUFIX", this->configs);

	this->filename_corr_files_dir = Configs_GetParam("FILENAME", "CORR_DUMP_DIR", this->configs);
	this->filename_corr_dump_filename = Configs_GetParam("FILENAME", "CORR_DUMP", this->configs);
	this->filename_corr_dump_sufix = Configs_GetParam("FILENAME", "CORR_DUMP_SUFIX", this->configs);

	this->filename_bcm_fullname = Configs_GetParam("FILENAME", "BCM", this->configs);
	this->filename_bai_fullname = Configs_GetParam("FILENAME", "BAI", this->configs);

	// creates log file names and open the file	
	Log(Configs_GetParam("FILENAME", "GLOBAL_LOG_DIR", this->configs), &this->log);
	
	// init dump flags	
	this->dumpSynth = Configs_GetParamInt("DUMP", "SYNTH", this->configs);
	this->dumpCorr = Configs_GetParamInt("DUMP", "CORR", this->configs);
	this->dumpBCM = Configs_GetParamInt("DUMP", "BCM", this->configs);
	this->dumpBAI = Configs_GetParamInt("DUMP", "BAI", this->configs);

	printf("return OK\n");

	return OK;
}

/*
 * Dss
 */
int TSI_Dss(tCube* outputCube, tCube* DssBAI, tCube* DssBCM, int kType, int iter, int simul, tTSI* this)
{
	tTimer thinTimer;


	char buffer[256];

	Log_WriteResult("Dss ...starting simulation ", (unsigned int) simul, this->log);

	
	if(kType == 5) {
			Timer_Start(&thinTimer);
			dsslib(this->DssParameters, this->DssModels, this->HardData, &(this->HardDataSize), DssBCM->values, DssBAI->values, outputCube->values);
			Timer_Stop(&thinTimer);
	} else if(kType == 1) {
			Timer_Start(&thinTimer);
			dsslib(this->DssParameters, this->DssModels, this->HardData, &(this->HardDataSize), NULL, NULL, outputCube->values);
			Timer_Stop(&thinTimer);
	}
	sprintf(buffer,"Dss simulation %d cpuTime",simul);
	Log_WriteActionTime(buffer, Timer_GetElapsedTimeInSeconds(&thinTimer), this->log);
	if(Configs_GetParam("DUMP","AI", this->configs)){
		sprintf(buffer,"%s%s_%d_%d%s",this->filename_ai_files_dir,this->filename_ai_filename,iter,simul,this->filename_ai_sufix);
		Cube_SaveToFile_name(buffer, this->root, outputCube);
	}
	return OK;
}


/*
 * Si
 */
int TSI_Si(tCube* aiCube, tCube** corrCube, int iterNum, int simulNum, tTSI* this)
{
	int ret_value;
	char filename[INTERNAL_FILENAME_MAXSIZE];		// filename buffer
	//char atoi_buffer[INTERNAL_ATOI_BUFFER_MAXSIZE]; // aux atoi buffer

	tTimer timer_step;		// time of each co-simulation
	tTimer timer_thin;		// time of a micro action

	// cubes used in simulation
	tCube* crCube;		// reflection coefecients cube
	tCube* synthCube;	// synthetic cube
	
	/*
	 * Simulation Code
	 */
	Timer_Reset(&timer_step);
	Timer_Start(&timer_step);

	Log_WriteMessage("Si ...", this->log);
	
	// Create Reflection Coeficients Cube
	Timer_Reset(&timer_thin);
	Timer_Start(&timer_thin);
	crCube = Cube_MakeReflectionCoefsCube(aiCube);
	Timer_Stop(&timer_thin);
	Log_WriteActionTime("Creating reflection coeficients cube cpuTime:", Timer_GetElapsedTimeInSeconds(&timer_thin), this->log);
		
	// Apply wavelet to Reflection Coeficients Cube
	// Sintetic cube creation, by convolution process
	Timer_Reset(&timer_thin);
	Timer_Start(&timer_thin);
	synthCube = Cube_MakeSyntheticCube(this->wavelet, crCube);

	/* we don't need cr cube anymore */
	free(crCube);
	Timer_Stop(&timer_thin);
	Log_WriteActionTime("Convolution process", Timer_GetElapsedTimeInSeconds(&timer_thin), this->log);
	
	// dump synthetic cube to file (optional)
	if (this->dumpSynth) {
		Timer_Reset(&timer_thin);
		Timer_Start(&timer_thin);
		sprintf(filename, "%s%s_%d_%d%s",this->filename_synth_files_dump_dir,this->filename_synth_dump_filename,iterNum,simulNum,this->filename_synth_dump_sufix);
		ret_value = Cube_SaveToFile_name(filename, 1, synthCube);
		if (ret_value < 0) {
			fprintf(stderr,"TSI.Si(): unable to dump synthCube, aborting Si()!\n");
			return ret_value;
		}
		Timer_Stop(&timer_thin);
		Log_WriteActionTime("Saving synth cube to file", Timer_GetElapsedTimeInSeconds(&timer_thin), this->log);

		Timer_Reset(&timer_thin);
		Timer_Start(&timer_thin);
		sprintf(filename, "%s%s_CR_%d_%d%s",this->filename_synth_files_dump_dir,this->filename_synth_dump_filename,iterNum,simulNum,this->filename_synth_dump_sufix);
		ret_value = Cube_SaveToFile_name(filename, 1, crCube);
		if (ret_value < 0) {
			fprintf(stderr,"TSI.Si(): unable to dump crCube, aborting Si()!\n");
			return ret_value;
		}
		Timer_Stop(&timer_thin);
		Log_WriteActionTime("Saving CR cube to file", Timer_GetElapsedTimeInSeconds(&timer_thin), this->log);
	}


	// Create Correlations Cube by comparing Sintetic Cube <-> Seismic Cube
	Timer_Reset(&timer_thin);
	Timer_Start(&timer_thin);
	*corrCube = MakeCorrelationsCube(this->seismic_cube, synthCube, this->layers);
	Timer_Stop(&timer_thin);
	Log_WriteActionTime("Creating correlations cube", Timer_GetElapsedTimeInSeconds(&timer_thin), this->log);

	// Calculate correlations average values
	Timer_Reset(&timer_thin);
	Timer_Start(&timer_thin);
	Cube_GenerateLayersCorr(this->layers, this->mask, (*corrCube));
	(*corrCube)->globalCorr = calculateGlobalCorr(this->seismic_cube, synthCube);
	/* we don't need synth cube anymore */
	free(synthCube);

	Timer_Stop(&timer_thin);
	Log_WriteActionTime("Calculating Correlations average values", Timer_GetElapsedTimeInSeconds(&timer_thin), this->log);

	// dump correlations cube to file (optional)
	if (this->dumpCorr) {
		Timer_Reset(&timer_thin);
		Timer_Start(&timer_thin);
		sprintf(filename, "%s%s_%d_%d%s",this->filename_corr_files_dir,this->filename_corr_dump_filename,iterNum,simulNum,this->filename_corr_dump_sufix);
		ret_value = Cube_SaveToFile_name(filename, 1, (*corrCube));
		if (ret_value < 0) {
			fprintf(stderr,"TSI.Si(): unable to dump corrCube, aborting Si()!\n");
			return ret_value;
		}
		Timer_Stop(&timer_thin);
		Log_WriteActionTime("Correlations saving to file", Timer_GetElapsedTimeInSeconds(&timer_thin), this->log);
	}
	

	Timer_Stop(&timer_step);
	Log_WriteActionTime("Si  cpuTime", Timer_GetElapsedTimeInSeconds(&timer_step), this->log);
	Log_WriteMessage("\tCorrelations Average for this simulation:", this->log);
	Log_WriteDump(Cube_CorrAverageDump((*corrCube)), this->log);

	/*
	 * End of Simulation Code
	 */
	return OK;
}


/*
 * Compare is the function that selects the best Cube.
 * 
 */
int TSI_Compare(tCube* aiCube, tCube* corrCube, struct aiGridData *bestAiGrid, tTSI* this)
{
	tTimer thinTimer;		// time of a micro action
	char file[255];

	// Rebuild BCM and BAI
	//Timer_Reset(&thinTimer);
	Log_WriteMessage("Compare ...", this->log);
	Timer_Start(&thinTimer);
	//UpdateResultCubes(*SiBCM, *SiBAI, corrCube, aiCube);

	if( bestAiGrid->globalCorr < corrCube->globalCorr) {
		sprintf(file,"%sbestAiGrid-%.3f.tsi",Configs_GetParam("FILENAME", "BEST_CUBES_DIR", this->configs),corrCube->globalCorr);
		if(Cube_SaveToFile_name(file, 1, aiCube) >= 0) { 
			// if saving the new bestAI doesn't fail
			// delete the older one
			sprintf(file,"%sbestAiGrid-%.3f.tsi",Configs_GetParam("FILENAME", "BEST_CUBES_DIR", this->configs),bestAiGrid->globalCorr);
		//	removeFile(file);
		}
		bestAiGrid->globalCorr = corrCube->globalCorr;
		bestAiGrid->layersCorr = Cube_copyLayersCorr(corrCube);
		bestAiGrid->layersNum = corrCube->layers_num;

                //FIX
		//sprintf(file,"Found a new best AI grid, saved has: %sbestAiGrid-%.3f.tsi\n",Configs_GetParam("FILENAME", "BEST_CUBES_DIR", this->configs),corrCube->globalCorr);
		Log_WriteMessage(file, this->log);
		printf("%s",file);
			
	}
		
	debug("Compare: deleting aiCube & corrCube");
	_Cube(aiCube);
	_Cube(corrCube);

	Timer_Stop(&thinTimer);
	Log_WriteActionTime("Updating BCM and BAI cubes", Timer_GetElapsedTimeInSeconds(&thinTimer), this->log);

	//FIX
	//printf("Compare: Iteration: %d\t Simulation: %d\t Global correlation: %f\n",this->iterNum,this->simulNum, corrCube->globalCorr);

/*
	//find best ai cube (global)
	if(*bestAiCube == NULL) {
		aiCube->layersCorr = Cube_copyLayersCorr(corrCube);

		aiCube->corrGlobal = corrCube->corrGlobal;
		aiCube->corrAvg = corrCube->corrAvg;
		*bestAiCube = aiCube;
		*bestGlobalSimulationNumber = simulNum; 
	} else if( (*bestAiCube)->corrGlobal < corrCube->corrGlobal) {
		debug("Compare: deleting bestAiCube");
		 free(*bestAiCube);

		aiCube->layersCorr = Cube_copyLayersCorr(corrCube);
		aiCube->corrGlobal = corrCube->corrGlobal;
		aiCube->corrAvg = corrCube->corrAvg;
		*bestAiCube = aiCube;
		*bestGlobalSimulationNumber = simulNum;

		Log_WriteMessage("\tFound BEST aiCube", this->log);
	} else {
		/* aiCube is not better than current best, lets delete it and its corrCube */
/*		free(aiCube);
	}
	debug("Compare: deleting corrCube");
	free(corrCube);
*/

	return OK;
}


/* FIXME: Dump should be revised to become a checkpointing mecanism.
 * idea is to save:
 * - BAI & BCM used in the current iteration by the coDSS.
 * - BAI & BCM being created by the current iteration.
 * - number of simulations to the end of current iteration.
 *   the best probably is to checkpoint at the end of each iteration, saving only
 *   the newtBAI & nextBCM and the iteration number
 */
int TSI_Dump(tCube* SiBAI, tCube** SiBCM, struct aiGridData *bestAiCube, 
		int iterNum, int bestGlobalSimulationNumber, tTSI* this)
{
	int ret_value;
	
	char BCMFilename[1024];
	char BAIFilename[1024];

	//Timer timer_step;		// time of each co-simulation
	tTimer thinTimer;		// time of a micro action

	// Save summary to file
	FILE* summaryFp;
	char summaryFilename[1024];

	// Calculate correlations average values
	Cube_GenerateLayersCorr(this->layers, this->mask, (*SiBCM));
	

	// BCM and BAI dump to file
	if (this->dumpBCM) {
		Timer_Reset(&thinTimer);
		Timer_Start(&thinTimer);
		sprintf(BCMFilename, "%s_%d.out", this->filename_bcm_fullname, iterNum);

		ret_value = Cube_SaveToFile_name(BCMFilename, this->root, (*SiBCM));
		if (ret_value < 0)
			return ret_value;
		Timer_Stop(&thinTimer);
		Log_WriteActionTime("Dumping BCM to file", Timer_GetElapsedTimeInSeconds(&thinTimer), this->log);
	// Dump correlations average values to file
		Log_WriteMessage("BCM Cube:", this->log);
		Log_WriteDump(Cube_CorrAverageDump((*SiBCM)), this->log);

	}

	if (this->dumpBAI) {
		Timer_Reset(&thinTimer);
		Timer_Start(&thinTimer);
		sprintf(BAIFilename, "%s_%d.out", this->filename_bai_fullname, iterNum);

		ret_value = Cube_SaveToFile_name(BAIFilename, 1, SiBAI);
		if (ret_value < 0)
			return ret_value;
		Timer_Stop(&thinTimer);
		Log_WriteActionTime("Dumping BAI to file", Timer_GetElapsedTimeInSeconds(&thinTimer), this->log);
	}
	// End of BCM and BAI dump to file

	
	sprintf(summaryFilename, "%sSummary.csv", Configs_GetParam("FILENAME", "TEMP_DIR", this->configs));
	summaryFp = fopen(summaryFilename, "a");
	
	if (summaryFp == NULL) {
		printf("Resume file: %s\n", summaryFilename);
		perror("Error creating Summary file");
		return ERROR;
	}

	//	FIXME: summary is broken! review its use/purpose/utility
	fprintf(summaryFp, "%d;%d;%.6Lf;%.6f\n", 
			iterNum, bestGlobalSimulationNumber, 
			bestAiCube->globalCorr,(*SiBCM)->globalCorr);
	fclose(summaryFp);

	return OK;
}

/*
 * Shutdown
 */
int TSI_Shutdown(tTSI* this)
{
	free(this->configs);
	free(this->wavelet);
	my_free(this->layers);
	free(this->seismic_cube);
	
	if (this->mask != NULL) {
		//printf("delete this->mask\n");
		free(this->mask);
	}
	free(this);

	// terminate
	return OK;
}

/*
 * reads the parameter file, creates an Array of Models and fills the Array of Parameters
 * returns NULL if an erros occurs
 */
void TSI_ReadDSSParameters(tTSI* this)
{
	char variogram[11];
	int nparms = 0;
	int tmp = 0;
	// harddata paramaters
	this->DssParameters[NVARI]	= Configs_GetParamInt("HARDDATA", "NVARI", this->configs);
	this->DssParameters[IXL]	= Configs_GetParamInt("HARDDATA", "IXL", this->configs);
	this->DssParameters[IYL]	= Configs_GetParamInt("HARDDATA", "IYL", this->configs);
	this->DssParameters[IZL]	= Configs_GetParamInt("HARDDATA", "IZL", this->configs);
	this->DssParameters[IVRL]	= Configs_GetParamInt("HARDDATA", "IVRL", this->configs);
	this->DssParameters[IWT]	= Configs_GetParamInt("HARDDATA", "IWT", this->configs);
	this->DssParameters[ISECVR]	= Configs_GetParamInt("HARDDATA", "ISECVR", this->configs);
	this->DssParameters[NVARI]	= Configs_GetParamInt("HARDDATA", "NVARI", this->configs);
	this->DssParameters[TMIN]	= (float) Configs_GetParamDouble("HARDDATA", "TMIN", this->configs);
	this->DssParameters[TMAX]	= (float) Configs_GetParamDouble("HARDDATA", "TMAX", this->configs);

	// harddata transformations
	this->DssParameters[ITRANS]	= Configs_GetParamInt("HDTRANS", "ITRANS", this->configs);
	this->DssParameters[ISMOOTH]    = Configs_GetParamInt("HDTRANS", "ISMOOTH", this->configs);
	this->DssParameters[ISVR]	= Configs_GetParamInt("HDTRANS", "ISVR", this->configs);
	this->DssParameters[ISWT]	= Configs_GetParamInt("HDTRANS", "ISWT", this->configs);
	this->DssParameters[ZMIN]	= (float) Configs_GetParamDouble("HDTRANS", "ZMIN", this->configs);
	this->DssParameters[ZMAX]	= (float) Configs_GetParamDouble("HDTRANS", "ZMAX", this->configs);
	this->DssParameters[LTAIL]	= Configs_GetParamInt("HDTRANS", "LTAIL", this->configs);
	this->DssParameters[LTPAR]	= (float) Configs_GetParamDouble("HDTRANS", "LTPAR", this->configs);
	this->DssParameters[UTAIL]	= Configs_GetParamInt("HDTRANS", "UTAIL", this->configs);
	this->DssParameters[UTPAR]	= (float) Configs_GetParamDouble("HDTRANS", "UTPAR", this->configs);

	// simulations quality
	this->DssParameters[NTRY]	= Configs_GetParamInt("QUALITY", "NTRY", this->configs);
	this->DssParameters[ICVAR]	= Configs_GetParamInt("QUALITY", "ICVAR", this->configs);
	this->DssParameters[ICMEAN]	= Configs_GetParamInt("QUALITY", "ICMEAN", this->configs);
	
	// grid parameters
	this->DssParameters[NX]		= this->x_values;
	this->DssParameters[NY]		= this->y_values;
	this->DssParameters[NZ]		= this->z_values;
	this->DssParameters[XMN]	= Configs_GetParamInt("GRID", "XCOORD", this->configs);
	this->DssParameters[YMN]	= Configs_GetParamInt("GRID", "YCOORD", this->configs);
	this->DssParameters[ZMN]	= Configs_GetParamInt("GRID", "ZCOORD", this->configs);
	this->DssParameters[XSIZ]	= Configs_GetParamInt("GRID", "XSIZE", this->configs);
	this->DssParameters[YSIZ]	= Configs_GetParamInt("GRID", "YSIZE", this->configs);
	this->DssParameters[ZSIZ]	= Configs_GetParamInt("GRID", "ZSIZE", this->configs);

	// mask
	this->DssParameters[NOSVALUE]= (float) Configs_GetParamDouble("MASK", "NULL_VALUE", this->configs);
	this->DssParameters[IMASK]	= Configs_GetParamInt("MASK", "USE_MASK", this->configs);

	// DSS search parameters
	this->DssParameters[NDMIN]	= Configs_GetParamInt("SEARCH", "NDMIN", this->configs);
	this->DssParameters[NDMAX]	= Configs_GetParamInt("SEARCH", "NDMAX", this->configs);
	this->DssParameters[NODMAX]	= Configs_GetParamInt("SEARCH", "NODMAX", this->configs);
	this->DssParameters[SSTRAT]	= Configs_GetParamInt("SEARCH", "SSTRAT", this->configs);
	this->DssParameters[MULTS]	= Configs_GetParamInt("SEARCH", "MULTS", this->configs);
	this->DssParameters[NMULTS]	= Configs_GetParamInt("SEARCH", "NMULTS", this->configs);
	this->DssParameters[NOCT]	= Configs_GetParamInt("SEARCH", "NOCT", this->configs);
	this->DssParameters[RADIUS]	= Configs_GetParamInt("SEARCH", "RADIUS", this->configs);
	this->DssParameters[RADIUS1]    = Configs_GetParamInt("SEARCH", "RADIUS1", this->configs);
	this->DssParameters[RADIUS2]    = Configs_GetParamInt("SEARCH", "RADIUS2", this->configs);
	this->DssParameters[SANG]	= Configs_GetParamInt("SEARCH", "SANG", this->configs);
	this->DssParameters[SANG1]	= Configs_GetParamInt("SEARCH", "SANG1", this->configs);
	this->DssParameters[SANG2]	= Configs_GetParamInt("SEARCH", "SANG2", this->configs);

	// kriging parameters
	this->DssParameters[KTYPE]	= Configs_GetParamInt("KRIG", "TYPE", this->configs);
	this->DssParameters[COLOCORR]   = Configs_GetParamInt("KRIG", "COLOCORR", this->configs);
	this->DssParameters[VARRED]	= Configs_GetParamInt("KRIG", "VARRED", this->configs);
	// softdata parameters
	this->DssParameters[NVARIL]	= Configs_GetParamInt("SOFT", "NVARIL", this->configs);
	this->DssParameters[ICOLLVM]    = Configs_GetParamInt("SOFT", "ICOLLVM", this->configs);
	// variogram models
	this->DssParameters[VARNUM]	= Configs_GetParamInt("VARIOGRAM", "NUMBER", this->configs);
	this->DssParameters[NUGGET]	= Configs_GetParamInt("VARIOGRAM", "NUGGET", this->configs);

	// malloc para o nº correcto de variogramas 
	this->DssModels = (float *) malloc(((unsigned int)this->DssParameters[VARNUM]) * 8 * sizeof(float));
	// e dpois copiar os dados para as estructuras
	/*
	char variogram[11];
	int nparms = 0;
	int tmp = 0;
	*/
	// BUG: só permite até 9 variogramas diferentes!

	for(; tmp < this->DssParameters[VARNUM]; tmp++){
		// append the number
		sprintf(variogram,"VARIOGRAM%d",tmp+1);
		
		nparms = tmp*8;

		this->DssModels[0+nparms]= Configs_GetParamInt(variogram, "TYPE", this->configs);
		this->DssModels[1+nparms]= (float) Configs_GetParamDouble(variogram, "COV", this->configs);
		this->DssModels[2+nparms]= Configs_GetParamInt(variogram, "ANG1", this->configs);
		this->DssModels[3+nparms]= Configs_GetParamInt(variogram, "ANG2", this->configs);
		this->DssModels[4+nparms]= Configs_GetParamInt(variogram, "ANG3", this->configs);
		this->DssModels[5+nparms]= Configs_GetParamInt(variogram, "AA", this->configs);
		this->DssModels[6+nparms]= Configs_GetParamInt(variogram, "AA1", this->configs);
		this->DssModels[7+nparms]= Configs_GetParamInt(variogram, "AA2", this->configs);
	}
}



/*
 * reads the hard data file to the HardData array
 * if theres is no hard data file, the HardDataSize is 0
 */
int TSI_ReadHardDataFile(char* Filename, tTSI* this)
{
	FILE *fp;
	unsigned int i=0;
	unsigned short int num = (unsigned short int)this->DssParameters[NVARI]; 

	double xCoord, yCoord, zCoord, value;
	//float xCoord, yCoord, zCoord, value;

	
	if((fp = fopen(Filename, "r")) == NULL) {
		printf("ReadHardDataFile(): Can´t open Hard Data file: %s\n", Filename);
		this->HardDataSize = 0;
		return ERROR;
	} else {
		// used to pass to Fortran routines, 
		// needed for the reading of the HardDataArray in a crazy while made of goto's
		this->HardDataSize = ( num * fileLineCount(Filename));

		//TODO: CHECK 4*.. MIGHT BE MISTAKE SHOULD BE 1
		this->HardData = (double *) malloc (4 * this->HardDataSize * sizeof(double));
		
		while(fscanf(fp, "%lf %lf %lf %lf", &xCoord, &yCoord, &zCoord, &value) != EOF) {
			this->HardData[i] = xCoord;
			this->HardData[i+1] = yCoord;
			this->HardData[i+2] = zCoord;
			this->HardData[i+3] = value;
			i += num;
		}
	}		

	fclose(fp);
	return OK; 

}



