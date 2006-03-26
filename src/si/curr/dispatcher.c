#include <stdio.h>
#include <stdlib.h>
#include <pthread.h> // for posix threads and mutexes
#include <semaphore.h> // for posix semaphores, used to control the number of running threads
#include "tsi.h"
#include "tsi-utils.h"
#include "timer.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))

void init_threading_environment(void); 
void * run_thread_seismic_inversion(void *);
void run_seismic_inversion(int , int);


struct thread_args {
	unsigned int simul;
	unsigned int iter;
};

tTSI *tsi;

static sem_t threads_semaphore;
static pthread_mutex_t compare_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_t	*thread;
my_sem_t save_mem_semaphore;

tCube *currBAI, *currBCM, *nextBAI, *nextBCM;

/* this structure holds 
 * data of the bestAi seen so far */
struct aiGridData bestAiData;


int bestGlobalSimulationNumber = -999;

int main(int argc, char *argv[])
{
	
	tTimer globalTimer;
	tTimer iterationTimer;

	tTimer thinTimer;
	char buf[1024];
	unsigned int iter,simul;
	int timeSeed;
	long lTime;
	char *buffer;
	int size;
	unsigned int i;

	//save bests to file
	char filenameBuffer[1024];

	printf("TSI %s - started.\n",VERSION);

	TSI(&tsi);

	printf("TSI alloced\n");

	// set a Starting point for the rand()
	lTime = time(NULL);
	timeSeed = (unsigned int) lTime/2;
	srandom(timeSeed);

	bestAiData.globalCorr = -999;
	
	//Initialization
	//
	if(argc == 2) {
		TSI_Initialization(argv[1], tsi);
	} else {
		TSI_Initialization(NULL, tsi);
	}

	Timer_Start(&globalTimer);

	printf("TSI inited\n");

	
	init_threading_environment();
	
	buffer = (char *) malloc(sizeof(char) * 256);
	sprintf(buffer,"TSI %s - Started.",VERSION);
	buffer[255] = '\0';
	Log_WriteMessage(buffer, tsi->log);
	size = tsi->x_values * tsi->y_values * tsi->z_values;
	sprintf(buffer,"Starting Simulation for grid with %u points",size);
	Log_WriteMessage(buffer, tsi->log);
	size /= (1024*1024);
	size *= sizeof(float);
	sprintf(buffer,"Memory use for simulation grid is %u Mbytes! ",size);
	Log_WriteMessage(buffer, tsi->log);
	/* 2 grids per process:
	 *  currBAI, currBCM 
	 * Plus 6 grids per Thread:
	 *  simulationCube, 
	 *  orderCube
	 *  covtable, covtable.ixnode, covtable.iynode, covtable.iznode
	 *  Plus
	 *  Wells
	 *  static overhead of config variables and other things.
	 */
	sprintf(buffer,"Predicted memory use is below %u Mbytes!",size*3+(size*6*tsi->ThreadsNumber));

	Log_WriteMessage(buffer, tsi->log);
	
	/* from here on, counts has iteration time */
	Timer_Start(&iterationTimer);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////// 

	/* lets generate new layers */
	if(generateRandomLayers(tsi->layers) != 0)
		exit(1);
	printLayers(tsi->layers);

	/* seismic cube updates number of layers, 
	 * and pre-computes data used on correlation calculus */ 
	Cube_LoadNewLayers(tsi->layers->number, tsi->seismic_cube);
	Cube_CalculateCorrData(tsi->layers, tsi->seismic_cube);

	/* nextBCM cube initialization */
	Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, CORR_CUBE, &nextBCM);
	Cube_Init(Configs_GetParamDouble("INIT", "BCM", tsi->configs), nextBCM);
	sprintf(buffer,"%snextBCM.tsi",Configs_GetParam("FILENAME", "TEMP_DIR", tsi->configs));
	Cube_associateFile(buffer, nextBCM);
	printf("nextBCM->name == %s\n",nextBCM->name);

	/* nextBAI cube initialization */
	Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, DEFAULT_CUBE, &nextBAI);
	Cube_Init(Configs_GetParamDouble("INIT", "BAI", tsi->configs), nextBAI); 
	sprintf(buffer,"%snextBAI.tsi",Configs_GetParam("FILENAME", "TEMP_DIR", tsi->configs));
	Cube_associateFile(buffer, nextBAI);
	printf("nextBAI->name == %s\n",nextBAI->name);

	// ktype = 1
	tsi->DssParameters[KTYPE] = 1;
	iter = 0;
	for(simul=1; simul <= tsi->TotalSimulationNumber; simul++) {
		run_seismic_inversion(iter, simul);
	}

	/* We must wait for the all the threads to finish */
	for(i = 0; i < tsi->TotalSimulationNumber; i++) {
		pthread_join(thread[i], NULL);
	}

	Timer_Stop(&iterationTimer);
	Log_WriteActionTime("First Iteration process time", Timer_GetElapsedTimeInSeconds(&iterationTimer), tsi->log);

	//printf("\nCalling Dump ... Iter=%d * Simul=%d\n", 0, (simul-1));
	TSI_Dump(nextBAI, &nextBCM, &bestAiData, 0, bestGlobalSimulationNumber, tsi);
	
	
	// ktype = 5
	tsi->DssParameters[KTYPE] = 5;
	for(iter=1; iter < tsi->TotalIterationNumber; iter++) {
		Log_WriteSeparator(tsi->log);
		Log_WriteMessage("Best Ai Correlations so far:", tsi->log);
		Log_WriteDump(DumpCorrData(&bestAiData), tsi->log);     //FIX
		printf("main loop (%d): Best Global  Correlation - %f\n", iter,bestAiData.globalCorr);

		Timer_Start(&iterationTimer);

		currBAI = nextBAI;
		currBCM = nextBCM;
		/* FIXME:  fazer um checkpoing da execucao aqui! */

		/* lets generate new layers */
		generateRandomLayers(tsi->layers);
		printLayers(tsi->layers);

		/* seismic cube updates number of layers, 
		 * and pre-computes data used on correlation calculus */ 
		Cube_LoadNewLayers(tsi->layers->number, tsi->seismic_cube);
		Cube_CalculateCorrData(tsi->layers, tsi->seismic_cube);
		
		// nextBCM cube initialization
		Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, CORR_CUBE, &nextBCM);
		Cube_Init(Configs_GetParamDouble("INIT", "BCM", tsi->configs), nextBCM);
		sprintf(buffer,"%snextBCM.tsi",Configs_GetParam("FILENAME", "TEMP_DIR", tsi->configs));
		Cube_associateFile(buffer, nextBCM);

		// nextBAI cube initialization
		Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, DEFAULT_CUBE, &nextBAI);
		Cube_Init(Configs_GetParamDouble("INIT", "BAI", tsi->configs), nextBAI);
		sprintf(buffer,"%snextBAI.tsi",Configs_GetParam("FILENAME", "TEMP_DIR", tsi->configs));
		Cube_associateFile(buffer, nextBAI);

		for(simul=1; simul <= tsi->TotalSimulationNumber; simul++) {
			run_seismic_inversion(iter, simul);
		}

		/* We must wait for the all the threads to finish */
		for(i = 0; i < tsi->TotalSimulationNumber; i++) {
			pthread_join(thread[i], NULL);
		}

		_Cube(currBAI);
		_Cube(currBCM);

//		printf("\nCalling Dump ... Iter=%d * Simul=%d\n", iter, (simul-1));
//		TSI_Dump(nextBAI, &nextBCM, bestAiCube, iter, bestGlobalSimulationNumber, tsi);
		
		Log_WriteResult("Iteration process time", (unsigned int) Timer_GetElapsedTimeInSeconds(&iterationTimer), tsi->log);
		Timer_Stop(&iterationTimer);
	}


	

	/* Save the best cubes has the result of this execution
	 */

	
	sprintf(buf,"TSI %s - Finished.\n\tBest AI grid is %s/bestAiGrid-%.3f.tsi\n\n",
			VERSION,Configs_GetParam("FILENAME", "BEST_CUBES_DIR", tsi->configs),bestAiData.globalCorr);
	printf(buf);
	Log_WriteMessage(buf, tsi->log);

	// TODO: substituir sprintf por snprintf
	//sprintf(filenameBuffer, "%s%s", Configs_GetParam("FILENAME", "BEST_CUBES_DIR", tsi->configs), "bestAiCube.out");
	//Cube_SaveToFile_name(filenameBuffer, tsi->root, bestAiCube);

	Log_WriteMessage("Best Ai Correlations:", tsi->log);
	Log_WriteDump(DumpCorrData(&bestAiData), tsi->log);
	printf("Best Ai Correlations:\n%s\n",DumpCorrData(&bestAiData));

	/* destroy the semaphore before exiting */
	sem_destroy(&threads_semaphore);	
	
	Timer_Stop(&globalTimer);
	Log_WriteActionTime("Global Execution Time", (unsigned long int) Timer_GetElapsedTimeInSeconds(&globalTimer), tsi->log);
	Log_Close(tsi->log);
	printf("TSI: Total Execution Time: %lu\n", (unsigned long int) Timer_GetElapsedTimeInSeconds(&globalTimer));

	TSI_Shutdown(tsi);
	free(nextBCM);
	//delete DssBAI;
	free(nextBAI);

}


void run_seismic_inversion(int i, int s)
{
	struct  thread_args * args = (struct thread_args *) malloc(sizeof(struct thread_args));
   	args->simul = s;
	args->iter = i;

	/* check number of allready running threads 
	 *  and wait for available thread
	 */
	sem_wait(&threads_semaphore);
	
	/* prepare data needed by DSS here to be "thread-safe" */
	/* lets create a thread to run DSS + SI + Compare */
	pthread_create(&thread[s-1], NULL, run_thread_seismic_inversion, (void *)args);
	
}	

void * run_thread_seismic_inversion(void * arg)
{
	struct thread_args * args;
	unsigned int iter, simul;
	tCube *corrCube;
	tCube *aiCube;
	tTimer *thinTimer;

	Timer(&thinTimer);
	args = (struct thread_args *) arg;
	iter = args->iter;
	simul = args->simul;

	Log_WriteSeparator(tsi->log);
	Log_WriteIterationNumber(iter, tsi->log);
	Log_WriteSimulationNumber(simul, tsi->log);

	/* the last thread starting a simulation should remove grids from memory */
	if(shouldSaveToFiles(&save_mem_semaphore)) {
		Log_WriteMessage("Removing Seismic, nextBCM, nextBAI from memory ...", tsi->log);
		Timer_Start(thinTimer);
		Cube_SaveToFile(nextBAI);
		Cube_SaveToFile(nextBCM);
		Timer_Stop(thinTimer);
		Log_WriteActionTime("Saving BCM and BAI grids to disk", Timer_GetElapsedTimeInSeconds(thinTimer), tsi->log);
		Timer_Start(thinTimer);
		Cube_freeGrid(tsi->seismic_cube);
		Cube_freeGrid(nextBAI);
		Cube_freeGrid(nextBCM);
		Timer_Stop(thinTimer);
		Log_WriteActionTime("Freeing Seismic, BAI, BCM grids", Timer_GetElapsedTimeInSeconds(thinTimer), tsi->log);
	}


	/*
	 * FIXME: Redu��o do consumo de mem�ria.
	 *     - na DSS o sim/aiCube e o tmp n�o s�o necess�rios ao mesmo tempo,
	 *       podemos _n�o_ alocar o sim, alocar o tmp, libertar o tmp, e depois alocar o sim.
	 *       (reduzindo efectivamente o consumo de mem�ria em 1 cubo)
	 */
	Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, DEFAULT_CUBE, &aiCube);

	if(iter == 0) {
		TSI_Dss(aiCube, NULL, NULL, 1, iter, simul, tsi); // ktype = 1
	} else {
		TSI_Dss(aiCube, currBAI, currBCM, 5, iter, simul, tsi); // ktype = 5
	}

	/* The first thread to finish dss() should read siesmic, nextBCM + nextBAI from storage */
	if(shouldRestoreFromFiles(&save_mem_semaphore)) {
		Log_WriteMessage("Reloading Seismic, nextBCM, nextBAI into memory ...", tsi->log);
		Timer_Start(thinTimer);
		Cube_allocGrid(tsi->seismic_cube);
		Cube_allocGrid(nextBAI);
		Cube_allocGrid(nextBCM);
		Timer_Stop(thinTimer);
		Log_WriteActionTime("Allocating space for Seismic, nextBCM, nextBAI grids ", Timer_GetElapsedTimeInSeconds(thinTimer), tsi->log);
		Timer_Start(thinTimer);
		Cube_ReadFromFile(tsi->seismic_cube);
		Cube_ReadFromFile(nextBCM);
		Cube_ReadFromFile(nextBAI);
		Timer_Stop(thinTimer);
		Log_WriteActionTime("Loaded Seismic, nextBCM, nextBAI grids ", Timer_GetElapsedTimeInSeconds(thinTimer), tsi->log);
		/* We must have a semaphore to "post" the arrival of the complete data to memory space. */
	}
	
	/* We must have a semaphore to "wait" for the arrival of the complete data to memory space. */

	TSI_Si(aiCube, &corrCube, iter, simul, tsi);
	printf("Compare: Iteration: %d\t Simulation: %d\t Global correlation: %f\n",iter,simul, corrCube->globalCorr);

	/* Compare must run with a mutual-exclusion lock!
	 * because Compare will update the best cubes.
	 * 	there is one alternative that is better for "big-smp" machines:
	 * 		- rw-lock per global-shared BestCube:
	 * 			- read lock for all,
	 * 		 	- aquire write lock when compare must update the global cube.
	 */
	pthread_mutex_lock(&compare_mutex);

	/* lets Compare & Update the BAI + BCM */
	Log_WriteMessage("Compare ...", tsi->log);
	Timer_Start(thinTimer);
	UpdateResultCubes(nextBCM, nextBAI, corrCube, aiCube);
	Timer_Stop(thinTimer);
	Log_WriteActionTime("Updating BCM and BAI cubes", Timer_GetElapsedTimeInSeconds(thinTimer), tsi->log);

	TSI_Compare(aiCube, corrCube, &bestAiData, tsi);

	pthread_mutex_unlock(&compare_mutex);

	/* lets post the finalization of this thread */
	sem_post(&threads_semaphore);

	/* lets free the struct with the thread arguments */
	my_free(arg);
	_Timer(thinTimer);

	return NULL;
}

void init_threading_environment()
{
	//tsi->compare_mutex = PTHREAD_MUTEX_INITIALIZER;
	unsigned int threadN;
	threadN = MIN(tsi->ThreadsNumber, tsi->TotalSimulationNumber);
	thread = (pthread_t *) malloc(sizeof(pthread_t) * tsi->TotalSimulationNumber);
	pthread_mutex_init(&save_mem_semaphore.lock, NULL);
	save_mem_semaphore.count = threadN;
	if( sem_init(&threads_semaphore, 0, threadN) )
		perror("InitThreadingEnvironment(): failed to inicialize Posix semaphore");	
}