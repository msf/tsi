#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dssconst.h"
#include "const.h"
#include "utils.h"
#include "tsi.h"
#include "tsi-utils.h"

TSI::TSI()
{
	
}

TSI::~TSI(void)
{
}

/*
 * Initialization
 *
 * It will read the parameters file, and set the acording environment variables.
 */
int TSI::Initialization(char *configFile)
{
	int ret_value;
		
	// read config file
	this->configs = new Configs();


	if(configFile != NULL)
		ret_value = this->configs->ReadConfigFile(configFile);
	else
		ret_value = this->configs->ReadConfigFile(CONFIGS_FILENAME);
	if (ret_value < 0) {
		printf("Error reading config file: %s\n", CONFIGS_FILENAME);
		exit(ret_value);
	}

	// GLOBAL parameters for TSI
	this->TotalIterationNumber = this->configs->GetParamInt("GLOBAL", "ITERATIONS");
	this->TotalSimulationNumber = this->configs->GetParamInt("GLOBAL", "SIMULATIONS");
	this->ThreadsNumber = (unsigned int) configs->GetParamInt("GLOBAL","THREADS");

	// read cube dimensions
	this->x_values = this->configs->GetParamInt("GRID", "XNUMBER");
	this->y_values = this->configs->GetParamInt("GRID", "YNUMBER");
	this->z_values = this->configs->GetParamInt("GRID", "ZNUMBER");
	
	// read sections and CORR related parameters.
	// TODO: variables that are used for "correlations" should be named: corr_name
	this->root = this->configs->GetParamInt("CORR", "BCM_ROOT");
	//
	// init layers definition
	this->layers = (struct layers_t *) malloc( sizeof(struct layers_t) );
	this->layers->minimum_number =  (unsigned short int) this->configs->GetParamInt("CORR", "LAYERS_MIN");
	this->layers->minimum_size = (unsigned short int) this->configs->GetParamInt("CORR", "LAYER_SIZE_MIN");
	this->layers->total_size = this->z_values;
	this->layers->number = 0; //sane default.

	// READ DSS RELATED PARAMETERS!
	this->ReadDSSParameters();



	// read HardData from the hard data file
	ret_value = this->ReadHardDataFile(this->configs->GetParam("HARDDATA", "FILENAME"));
	if (ret_value < 0) {
		printf("Error reading Hard Data file: %s\n", this->configs->GetParam("HARDDATA", "FILENAME"));
		return ret_value;
	}

	// Read Wavelet from file to memory, including wavelet parameters
	this->wavelet = new Wavelet(this->configs);
	ret_value = this->wavelet->ReadFromFile(this->configs->GetParam("WAVELET", "FILENAME"));
	if (ret_value < 0) {
		printf("Error reading Wavelet file: %s\n", this->configs->GetParam("WAVELET", "FILENAME"));
		return ret_value;
	}

	// Read Seismic Cube from file to memory
	this->seismic_cube = new Cube(this->x_values, this->y_values, this->z_values, 1, CORR_CUBE); //the 1 is because we don't know yet how many layers we'll use
	ret_value = this->seismic_cube->ReadFromFile(this->configs->GetParam("SEISMIC", "FILENAME"));
	if (ret_value < 0) {
		printf("Error reading Seismic Cube file: %s\n", this->configs->GetParam("SEISMIC", "FILENAME"));
		return ret_value;
	}

	// init statics values for the seismic cube (values precomputed for globalCorr computation
	this->seismic_cube->CalculateStaticCorrSums();

	this->global_null_value = this->configs->GetParamDouble("MASK", "NULL_VALUE");
	// read mask from file
	if ((this->configs->GetParamInt("MASK", "USE_MASK")) == 1) {
		this->mask = new Cube(this->x_values, this->y_values, this->z_values, this->layers->number, DEFAULT_CUBE);
		ret_value = this->mask->ReadFromFile(this->configs->GetParam("MASK", "FILENAME"));
		if (ret_value < 0){
			printf("Error reading Mask file: %s\n", this->configs->GetParam("MASK", "FILENAME"));
			return ret_value;
		}
	} else {
		this->mask = NULL;
	}


	this->filename_ai_files_dir = this->configs->GetParam("FILENAME", "TEMP_DIR");
	this->filename_ai_filename = "aiCube";
	this->filename_ai_sufix = ".out";
	

	// init all working paths (AI cubes, Synth cubes, Corr cubes, BCM, BAI, logs)
	this->filename_synth_files_dump_dir = this->configs->GetParam("FILENAME", "SYNTH_DUMP_DIR");
	this->filename_synth_dump_filename = this->configs->GetParam("FILENAME", "SYNTH_DUMP");
	this->filename_synth_dump_sufix = this->configs->GetParam("FILENAME", "SYNTH_DUMP_SUFIX");

	this->filename_corr_files_dir = this->configs->GetParam("FILENAME", "CORR_DUMP_DIR");
	this->filename_corr_dump_filename = this->configs->GetParam("FILENAME", "CORR_DUMP");
	this->filename_corr_dump_sufix = this->configs->GetParam("FILENAME", "CORR_DUMP_SUFIX");

	this->filename_bcm_fullname = this->configs->GetParam("FILENAME", "BCM");
	this->filename_bai_fullname = this->configs->GetParam("FILENAME", "BAI");

	// creates log file names and open the file	
	this->log = new Log(this->configs->GetParam("FILENAME", "GLOBAL_LOG_DIR"));
	
	// init dump flags	
	this->dumpSynth = this->configs->GetParamInt("DUMP", "SYNTH");
	this->dumpCorr = this->configs->GetParamInt("DUMP", "CORR");
	this->dumpBCM = this->configs->GetParamInt("DUMP", "BCM");
	this->dumpBAI = this->configs->GetParamInt("DUMP", "BAI");


	return OK;
}

/*
 * Dss
 */
int TSI::Dss(Cube* outputCube, Cube* DssBAI, Cube* DssBCM, int kType, int iter, int simul)
{
	Timer thinTimer;


	char buffer[256];

	this->log->WriteResult("Dss ...starting simulation ", (unsigned int) simul);

	
	if(kType == 5) {
		if(this->configs->GetParamInt("MASK", "USE_MASK")) {
			thinTimer.Reset();
			thinTimer.Start();
			dssdll(this->DssParameters, this->DssModels, this->HardData, &(this->HardDataSize), DssBCM->values, DssBAI->values, this->mask->values, outputCube->values);
			thinTimer.Stop();
		} else {
			thinTimer.Reset();
			thinTimer.Start();
			dssdll(this->DssParameters, this->DssModels, this->HardData, &(this->HardDataSize), DssBCM->values, DssBAI->values, NULL, outputCube->values);
			thinTimer.Stop();
		}
	} else if(kType == 1) {
		if(this->configs->GetParamInt("MASK", "USE_MASK")) {
			thinTimer.Reset();
			thinTimer.Start();
			dssdll(this->DssParameters, this->DssModels, this->HardData, &(this->HardDataSize), NULL, NULL, this->mask->values, outputCube->values);
			thinTimer.Stop();
		} else {
			thinTimer.Reset();
			thinTimer.Start();
			dssdll(this->DssParameters, this->DssModels, this->HardData, &(this->HardDataSize), NULL, NULL, NULL, outputCube->values);
			thinTimer.Stop();
		}
	}
	sprintf(buffer,"Dss simulation %d cpuTime",simul);
	this->log->WriteActionTime(buffer, thinTimer.GetElapsedTimeInSeconds());
	if(configs->GetParam("DUMP","AI")){
		sprintf(buffer,"%s%s_%d_%d%s",filename_ai_files_dir,filename_ai_filename,iter,simul,filename_ai_sufix);
		outputCube->SaveToFile(buffer, mask, global_null_value, root);
	}
	return OK;
}


/*
 * Si
 */
int TSI::Si(Cube* aiCube, Cube** corrCube, int iterNum, int simulNum)
{
	int ret_value;
	char filename[INTERNAL_FILENAME_MAXSIZE];		// filename buffer
	//char atoi_buffer[INTERNAL_ATOI_BUFFER_MAXSIZE]; // aux atoi buffer

	Timer timer_step;		// time of each co-simulation
	Timer timer_thin;		// time of a micro action

	// cubes used in simulation
	Cube* crCube;		// reflection coefecients cube
	Cube* synthCube;	// synthetic cube
	
	/*
	 * Simulation Code
	 */
	timer_step.Reset();
	timer_step.Start();

	this->log->WriteMessage("Si ...");
	
	// Create Reflection Coeficients Cube
	timer_thin.Reset();
	timer_thin.Start();
	crCube = aiCube->MakeReflectionCoefsCube();
	timer_thin.Stop();
	this->log->WriteActionTime("Creating reflection coeficients cube cpuTime:", timer_thin.GetElapsedTimeInSeconds());
		
	// Apply wavelet to Reflection Coeficients Cube
	// Sintetic cube creation, by convolution process
	timer_thin.Reset();
	timer_thin.Start();
	synthCube = crCube->MakeSyntheticCube(this->wavelet);

	/* we don't need cr cube anymore */
	delete crCube;
	timer_thin.Stop();
	this->log->WriteActionTime("Convolution process", timer_thin.GetElapsedTimeInSeconds());
	
	// dump synthetic cube to file (optional)
	if (this->dumpSynth) {
		timer_thin.Reset();
		timer_thin.Start();
		sprintf(filename, "%s%s_%d_%d%s",this->filename_synth_files_dump_dir,this->filename_synth_dump_filename,iterNum,simulNum,this->filename_synth_dump_sufix);
		ret_value = synthCube->SaveToFile(filename, this->mask, this->global_null_value, 1);
		if (ret_value < 0) {
			fprintf(stderr,"TSI.Si(): unable to dump synthCube, aborting Si()!\n");
			return ret_value;
		}
		timer_thin.Stop();
		this->log->WriteActionTime("Saving synth cube to file", timer_thin.GetElapsedTimeInSeconds());

		timer_thin.Reset();
		timer_thin.Start();
		sprintf(filename, "%s%s_CR_%d_%d%s",this->filename_synth_files_dump_dir,this->filename_synth_dump_filename,iterNum,simulNum,this->filename_synth_dump_sufix);
		ret_value = crCube->SaveToFile(filename, this->mask, this->global_null_value, 1);
		if (ret_value < 0) {
			fprintf(stderr,"TSI.Si(): unable to dump crCube, aborting Si()!\n");
			return ret_value;
		}
		timer_thin.Stop();
		this->log->WriteActionTime("Saving CR cube to file", timer_thin.GetElapsedTimeInSeconds());
	}


	// Create Correlations Cube by comparing Sintetic Cube <-> Seismic Cube
	timer_thin.Reset();
	timer_thin.Start();
	*corrCube = MakeCorrelationsCube(this->seismic_cube, synthCube, this->layers);
	timer_thin.Stop();
	this->log->WriteActionTime("Creating correlations cube", timer_thin.GetElapsedTimeInSeconds());

	// Calculate correlations average values
	timer_thin.Reset();
	timer_thin.Start();
	(*corrCube)->GenerateLayersCorr(this->layers, this->mask);
	(*corrCube)->corrGlobal = calculateGlobalCorr(this->seismic_cube, synthCube);
	/* we don't need synth cube anymore */
	delete synthCube;

	timer_thin.Stop();
	this->log->WriteActionTime("Calculating Correlations average values", timer_thin.GetElapsedTimeInSeconds());

	// dump correlations cube to file (optional)
	if (this->dumpCorr) {
		timer_thin.Reset();
		timer_thin.Start();
		sprintf(filename, "%s%s_%d_%d%s",this->filename_corr_files_dir,this->filename_corr_dump_filename,iterNum,simulNum,this->filename_corr_dump_sufix);
		ret_value = (*corrCube)->SaveToFile(filename, this->mask, this->global_null_value, 1);
		if (ret_value < 0) {
			fprintf(stderr,"TSI.Si(): unable to dump corrCube, aborting Si()!\n");
			return ret_value;
		}
		timer_thin.Stop();
		this->log->WriteActionTime("Correlations saving to file", timer_thin.GetElapsedTimeInSeconds());
	}
	

	timer_step.Stop();
	this->log->WriteActionTime("Si  cpuTime", timer_step.GetElapsedTimeInSeconds());
	this->log->WriteMessage("\tCorrelations Average for this simulation:");
	this->log->WriteDump((*corrCube)->CorrAverageDump());

	/*
	 * End of Simulation Code
	 */
	return OK;
}


/*
 * Compare is the function that selects the best Cube.
 * 
 */
int TSI::Compare(Cube* aiCube, Cube* corrCube, Cube** SiBAI, Cube** SiBCM, Cube** bestAiCube, int* bestGlobalSimulationNumber, int iterNum, int simulNum)
{
	Timer thinTimer;		// time of a micro action

	// Rebuild BCM and BAI
	thinTimer.Reset();
	this->log->WriteMessage("Compare ...");
	thinTimer.Start();
	UpdateResultCubes( *SiBCM, *SiBAI, corrCube, aiCube);
	thinTimer.Stop();
	this->log->WriteActionTime("Updating BCM and BAI cubes", thinTimer.GetElapsedTimeInSeconds());

	printf("Compare: Iteration: %d\t Simulation: %d\t Global correlation: %f\n",iterNum,simulNum, corrCube->corrGlobal);

	//find best ai cube (global)
	if(*bestAiCube == NULL) {
		aiCube->layersCorr = corrCube->copyLayersCorr();

		aiCube->corrGlobal = corrCube->corrGlobal;
		aiCube->corrAvg = corrCube->corrAvg;
		*bestAiCube = aiCube;
		*bestGlobalSimulationNumber = simulNum; 
	} else if( (*bestAiCube)->corrGlobal < corrCube->corrGlobal) {
		debug("Compare: deleting bestAiCube");
		delete (*bestAiCube);

		aiCube->layersCorr = corrCube->copyLayersCorr();
		aiCube->corrGlobal = corrCube->corrGlobal;
		aiCube->corrAvg = corrCube->corrAvg;
		*bestAiCube = aiCube;
		*bestGlobalSimulationNumber = simulNum;

		this->log->WriteMessage("\tFound BEST aiCube");
	} else {
		/* aiCube is not better than current best, lets delete it and its corrCube */
		delete aiCube;
	}
	debug("Compare: deleting corrCube");
	delete corrCube;

	return OK;
}

int TSI::Dump(Cube* SiBAI, Cube** SiBCM, Cube* bestAiCube, 
		int iterNum, int bestGlobalSimulationNumber)
{
	int ret_value;
	
	char BCMFilename[1024];
	char BAIFilename[1024];

	//Timer timer_step;		// time of each co-simulation
	Timer thinTimer;		// time of a micro action

	// Calculate correlations average values
	(*SiBCM)->GenerateLayersCorr(this->layers, this->mask);
	

	// BCM and BAI dump to file
	if (this->dumpBCM) {
		thinTimer.Reset();
		thinTimer.Start();
		sprintf(BCMFilename, "%s_%d.out", this->filename_bcm_fullname, iterNum);

		ret_value = (*SiBCM)->SaveToFile(BCMFilename, this->mask, this->global_null_value, this->root);
		if (ret_value < 0)
			return ret_value;
		thinTimer.Stop();
		this->log->WriteActionTime("Dumping BCM to file", thinTimer.GetElapsedTimeInSeconds());
	// Dump correlations average values to file
		this->log->WriteMessage("BCM Cube:");
		this->log->WriteDump((*SiBCM)->CorrAverageDump());

	}

	if (this->dumpBAI) {
		thinTimer.Reset();
		thinTimer.Start();
		sprintf(BAIFilename, "%s_%d.out", this->filename_bai_fullname, iterNum);

		ret_value = SiBAI->SaveToFile(BAIFilename, this->mask, GLOBAL_NULL_VALUE, 1);	
		if (ret_value < 0)
			return ret_value;
		thinTimer.Stop();
		this->log->WriteActionTime("Dumping BAI to file", thinTimer.GetElapsedTimeInSeconds());
	}
	// End of BCM and BAI dump to file

	// Save summary to file
	FILE* summaryFp;
	char summaryFilename[1024];
	
	sprintf(summaryFilename, "%sSummary.csv",this->configs->GetParam("FILENAME", "TEMP_DIR"));
	summaryFp = fopen(summaryFilename, "a");
	
	if (summaryFp == NULL) {
		printf("Resume file: %s\n", summaryFilename);
		perror("Error creating Summary file");
		return ERROR;
	}

	//fprintf(summaryFp, "%03d;%03d;%.6Lf;%03d;%03d;%.6Lf;%.6Lf\n", iterNum, bestGlobalSimulationNumber, 
	//		bestAiCubeGlobal->GetGlobalCorrAverage(), bestSectionSimulationNumber, this->reserv_section,
	//		bestAiCubeSection->GetSectionCorrAverage(this->reserv_section), SiBCM->corrGlobal);
	//	FIXME: summary is broken! review its use/purpose/utility
	fprintf(summaryFp, "%d;%d;%.6Lf;%.6f\n", 
			iterNum, bestGlobalSimulationNumber, 
			bestAiCube->corrGlobal,(*SiBCM)->corrGlobal);
	fclose(summaryFp);

	return OK;
}

/*
 * Shutdown
 */
int TSI::Shutdown()
{
	delete this->configs;
	delete this->wavelet;
	my_free(this->layers);
	delete this->seismic_cube;
	
	if (this->mask != NULL) {
		//printf("delete this->mask\n");
		delete this->mask;
	}

	// terminate
	return OK;
}

/*
 * reads the parameter file, creates an Array of Models and fills the Array of Parameters
 * returns NULL if an erros occurs
 */
void TSI::ReadDSSParameters()
{
	// harddata paramaters
	this->DssParameters[NVARI]	= this->configs->GetParamInt("HARDDATA", "NVARI");
	this->DssParameters[IXL]	= this->configs->GetParamInt("HARDDATA", "IXL");
	this->DssParameters[IYL]	= this->configs->GetParamInt("HARDDATA", "IYL");
	this->DssParameters[IZL]	= this->configs->GetParamInt("HARDDATA", "IZL");
	this->DssParameters[IVRL]	= this->configs->GetParamInt("HARDDATA", "IVRL");
	this->DssParameters[IWT]	= this->configs->GetParamInt("HARDDATA", "IWT");
	this->DssParameters[ISECVR]	= this->configs->GetParamInt("HARDDATA", "ISECVR");
	this->DssParameters[NVARI]	= this->configs->GetParamInt("HARDDATA", "NVARI");
	this->DssParameters[TMIN]	= (float) this->configs->GetParamDouble("HARDDATA", "TMIN");
	this->DssParameters[TMAX]	= (float) this->configs->GetParamDouble("HARDDATA", "TMAX");

	// harddata transformations
	this->DssParameters[ITRANS]	= this->configs->GetParamInt("HDTRANS", "ITRANS");
	this->DssParameters[ISMOOTH]= this->configs->GetParamInt("HDTRANS", "ISMOOTH");
	this->DssParameters[ISVR]	= this->configs->GetParamInt("HDTRANS", "ISVR");
	this->DssParameters[ISWT]	= this->configs->GetParamInt("HDTRANS", "ISWT");
	this->DssParameters[ZMIN]	= (float) this->configs->GetParamDouble("HDTRANS", "ZMIN");
	this->DssParameters[ZMAX]	= (float) this->configs->GetParamDouble("HDTRANS", "ZMAX");
	this->DssParameters[LTAIL]	= this->configs->GetParamInt("HDTRANS", "LTAIL");
	this->DssParameters[LTPAR]	= (float) this->configs->GetParamDouble("HDTRANS", "LTPAR");
	this->DssParameters[UTAIL]	= this->configs->GetParamInt("HDTRANS", "UTAIL");
	this->DssParameters[UTPAR]	= (float) this->configs->GetParamDouble("HDTRANS", "UTPAR");

	// simulations quality
	this->DssParameters[NTRY]	= this->configs->GetParamInt("QUALITY", "NTRY");
	this->DssParameters[ICVAR]	= this->configs->GetParamInt("QUALITY", "ICVAR");
	this->DssParameters[ICMEAN]	= this->configs->GetParamInt("QUALITY", "ICMEAN");
	
	// grid parameters
	this->DssParameters[NX]		= this->x_values;
	this->DssParameters[NY]		= this->y_values;
	this->DssParameters[NZ]		= this->z_values;
	this->DssParameters[XMN]	= this->configs->GetParamInt("GRID", "XCOORD");
	this->DssParameters[YMN]	= this->configs->GetParamInt("GRID", "YCOORD");
	this->DssParameters[ZMN]	= this->configs->GetParamInt("GRID", "ZCOORD");
	this->DssParameters[XSIZ]	= this->configs->GetParamInt("GRID", "XSIZE");
	this->DssParameters[YSIZ]	= this->configs->GetParamInt("GRID", "YSIZE");
	this->DssParameters[ZSIZ]	= this->configs->GetParamInt("GRID", "ZSIZE");

	// mask
	this->DssParameters[NOSVALUE]= (float) this->configs->GetParamDouble("MASK", "NULL_VALUE");
	this->DssParameters[IMASK]	= this->configs->GetParamInt("MASK", "USE_MASK");

	// DSS search parameters
	this->DssParameters[NDMIN]	= this->configs->GetParamInt("SEARCH", "NDMIN");
	this->DssParameters[NDMAX]	= this->configs->GetParamInt("SEARCH", "NDMAX");
	this->DssParameters[NODMAX]	= this->configs->GetParamInt("SEARCH", "NODMAX");
	this->DssParameters[SSTRAT]	= this->configs->GetParamInt("SEARCH", "SSTRAT");
	this->DssParameters[MULTS]	= this->configs->GetParamInt("SEARCH", "MULTS");
	this->DssParameters[NMULTS]	= this->configs->GetParamInt("SEARCH", "NMULTS");
	this->DssParameters[NOCT]	= this->configs->GetParamInt("SEARCH", "NOCT");
	this->DssParameters[RADIUS]	= this->configs->GetParamInt("SEARCH", "RADIUS");
	this->DssParameters[RADIUS1]= this->configs->GetParamInt("SEARCH", "RADIUS1");
	this->DssParameters[RADIUS2]= this->configs->GetParamInt("SEARCH", "RADIUS2");
	this->DssParameters[SANG]	= this->configs->GetParamInt("SEARCH", "SANG");
	this->DssParameters[SANG1]	= this->configs->GetParamInt("SEARCH", "SANG1");
	this->DssParameters[SANG2]	= this->configs->GetParamInt("SEARCH", "SANG2");

	// kriging parameters
	this->DssParameters[KTYPE]	= this->configs->GetParamInt("KRIG", "TYPE");
	this->DssParameters[COLOCORR]= this->configs->GetParamInt("KRIG", "COLOCORR");
	this->DssParameters[VARRED]	= this->configs->GetParamInt("KRIG", "VARRED");
	// softdata parameters
	this->DssParameters[NVARIL]	= this->configs->GetParamInt("SOFT", "NVARIL");
	this->DssParameters[ICOLLVM]= this->configs->GetParamInt("SOFT", "ICOLLVM");
	// variogram models
	this->DssParameters[VARNUM]	= this->configs->GetParamInt("VARIOGRAM", "NUMBER");
	this->DssParameters[NUGGET]	= this->configs->GetParamInt("VARIOGRAM", "NUGGET");

	// malloc para o nº correcto de variogramas 
	this->DssModels = (float *) malloc(((unsigned int)this->DssParameters[VARNUM]) * 8 * sizeof(float));
	// e dpois copiar os dados para as estructuras
	char variogram[11];
	int nparms = 0;
	int tmp = 0;
	// BUG: só permite até 9 variogramas diferentes!

	for(; tmp < this->DssParameters[VARNUM]; tmp++){
		// append the number
		sprintf(variogram,"VARIOGRAM%d",tmp+1);
		
		nparms = tmp*8;

		this->DssModels[0+nparms]= this->configs->GetParamInt(variogram, "TYPE");
		this->DssModels[1+nparms]= (float) this->configs->GetParamDouble(variogram, "COV");
		this->DssModels[2+nparms]= this->configs->GetParamInt(variogram, "ANG1");
		this->DssModels[3+nparms]= this->configs->GetParamInt(variogram, "ANG2");
		this->DssModels[4+nparms]= this->configs->GetParamInt(variogram, "ANG3");
		this->DssModels[5+nparms]= this->configs->GetParamInt(variogram, "AA");
		this->DssModels[6+nparms]= this->configs->GetParamInt(variogram, "AA1");
		this->DssModels[7+nparms]= this->configs->GetParamInt(variogram, "AA2");
	}
}



/*
 * reads the hard data file to the HardData array
 * if theres is no hard data file, the HardDataSize is 0
 */
int TSI::ReadHardDataFile(char* Filename)
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



