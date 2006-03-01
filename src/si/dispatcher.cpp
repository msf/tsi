#include <stdio.h>
#include <stdlib.h>
#include <pthread.h> // for posix threads and mutexes
#include <semaphore.h> // for posix semaphores, used to control the number of running threads

#include "tsi.h"
#include "tsi-utils.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))

void threading_init_environment(void); 
void * run_thread_seismic_inversion(void *);
void run_seismic_inversion(int , int);


struct thread_args {
	unsigned int simul;
	unsigned int iter;
};

TSI *tsi;
static sem_t threads_semaphore;
static pthread_mutex_t compare_mutex =PTHREAD_MUTEX_INITIALIZER;
static pthread_t	*thread;

Cube *DssBAI, *DssBCM, *SiBAI, *SiBCM;

/* the best Ai cube for the whole simulation
 * the correlation value is in member corrGlobal */
Cube *bestAiCube = NULL;

int bestGlobalSimulationNumber = -999;

int main(int argc, char *argv[])
{
	
	printf("TSI %s - started.\n",VERSION);

	tsi = new TSI();
	Timer globalTimer;
	Timer iterationTimer;

	Timer thinTimer;

	unsigned int iter,simul;
	int timeSeed;
	long lTime;

	// set a starting point for the rand()
	lTime = time(NULL);
	timeSeed = (unsigned int) lTime/2;
	srandom(timeSeed);

	
	//Initialization
	//
	if(argc == 2) {
		tsi->Initialization(argv[1]);
	} else {
		tsi->Initialization(NULL);
	}

	globalTimer.Reset();
	globalTimer.Start();

	
	threading_init_environment();
	
	char *buffer = (char *) malloc(sizeof(char) * 256);
	sprintf(buffer,"TSI %s - started.",VERSION);
	buffer[255] = '\0';
	tsi->log->WriteMessage(buffer);
	int size = tsi->x_values * tsi->y_values * tsi->z_values;
	sprintf(buffer,"Starting Simulation for cube with %u points",size);
	tsi->log->WriteMessage(buffer);
	size /= (1024*1024);
	size *= sizeof(float);
	sprintf(buffer,"Memory use for simulation cube is %u Mbytes! ",size);
	tsi->log->WriteMessage(buffer);
	/* 4 cubes per process:
	 *  dssBAI, dssBCM 
	 *  siBAI, siBCM 
	 *  BestAIGlobal, 
	 *  SeismicCube/HardData
	 * Plus 7 cubes per Thread:
	 *  simulationCube, 
	 *  orderCube, tmpCube, (tmp é arrasável) 
	 *  covtable, covtable.ixnode, covtable.iynode, covtable.iznode
	 * FIXME: please check the number of cubes used on TSI, DSS was allready checked
	 */
	sprintf(buffer,"Predicted Peak memory use is %u Mbytes!",size*6+(size*6*tsi->ThreadsNumber));
	tsi->log->WriteMessage(buffer);
	
	/* from here on, counts has iteration time */
	iterationTimer.Reset();
	iterationTimer.Start();

	/* lets generate new layers */
	if(generateRandomLayers(tsi->layers) != 0)
		exit(1);
	printLayers(tsi->layers);

	/* seismic cube updates number of layers, 
	 * and pre-computes data used on correlation calculus */ 
	tsi->seismic_cube->LoadNewLayers(tsi->layers->number);
	tsi->seismic_cube->CalculateCorrData(tsi->layers);

	/* SiBCM cube initialization */
	SiBCM = new Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, CORR_CUBE);
	SiBCM->Init(tsi->configs->GetParamDouble("INIT", "BCM"));

	/* SiBAI cube initialization */
	SiBAI = new Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, DEFAULT_CUBE);
	SiBAI->Init(tsi->configs->GetParamDouble("INIT", "BAI")); 

	// ktype = 1
	tsi->DssParameters[KTYPE] = 1;
	iter = 0;
	for(simul=1; simul <= tsi->TotalSimulationNumber; simul++) {
		run_seismic_inversion(iter, simul);
	}

	/* We must wait for the all the threads to finish */
	for(unsigned int i = 0; i < tsi->TotalSimulationNumber; i++) {
		pthread_join(thread[i], NULL);
	}

	iterationTimer.Stop();
	tsi->log->WriteActionTime("First Iteration process time", iterationTimer.GetElapsedTimeInSeconds());

	//printf("\nCalling Dump ... Iter=%d * Simul=%d\n", 0, (simul-1));
	tsi->Dump(SiBAI, &SiBCM, bestAiCube, 0, bestGlobalSimulationNumber);
	
	
	// ktype = 5
	tsi->DssParameters[KTYPE] = 5;
	for(iter=1; iter < tsi->TotalIterationNumber; iter++) {
		tsi->log->WriteSeparator();
		tsi->log->WriteMessage("Best Ai Correlations so far:");
		tsi->log->WriteDump(bestAiCube->CorrAverageDump());
		printf("main loop (%d): Best Global  Correlation - %f\n", iter,bestAiCube->corrGlobal);

		iterationTimer.Reset();
		iterationTimer.Start();

		DssBAI = SiBAI;
		DssBCM = SiBCM;

		/* lets generate new layers */
		generateRandomLayers(tsi->layers);
		printLayers(tsi->layers);

		/* seismic cube updates number of layers, 
		 * and pre-computes data used on correlation calculus */ 
		tsi->seismic_cube->LoadNewLayers(tsi->layers->number);
		tsi->seismic_cube->CalculateCorrData(tsi->layers);
		
		// SiBCM cube initialization
		SiBCM = new Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, CORR_CUBE);
		SiBCM->Init(tsi->configs->GetParamDouble("INIT", "BCM"));
		
		// SiBAI cube initialization
		SiBAI = new Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, DEFAULT_CUBE);
		SiBAI->Init(tsi->configs->GetParamDouble("INIT", "BAI"));

		for(simul=1; simul <= tsi->TotalSimulationNumber; simul++) {
			run_seismic_inversion(iter, simul);
		}

		/* We must wait for the all the threads to finish */
		for(unsigned int i = 0; i < tsi->TotalSimulationNumber; i++) {
			pthread_join(thread[i], NULL);
		}

		delete DssBAI;
		delete DssBCM;

//		printf("\nCalling Dump ... Iter=%d * Simul=%d\n", iter, (simul-1));
		tsi->Dump(SiBAI, &SiBCM, bestAiCube, iter, bestGlobalSimulationNumber);
		
		tsi->log->WriteResult("Iteration process time", (unsigned int) iterationTimer.GetElapsedTimeInSeconds());
		iterationTimer.Stop();
	}


	//save bests to file
	char filenameBuffer[1024];
	

	/* Save the best cubes has the result of this execution
	 */

	// TODO: substituir sprintf por snprintf
	sprintf(filenameBuffer, "%s%s", tsi->configs->GetParam("FILENAME", "BEST_CUBES_DIR"), "bestAiCube.out");
	bestAiCube->SaveToFile(filenameBuffer, tsi->mask, tsi->global_null_value, tsi->root);

	tsi->log->WriteMessage("Best Ai Correlations:");
	tsi->log->WriteDump(bestAiCube->CorrAverageDump());
	printf("Best Ai Correlations:\n%s\n",bestAiCube->CorrAverageDump());

	/* destroy the semaphore before exiting */
	sem_destroy(&threads_semaphore);	
	
	globalTimer.Stop();
	tsi->log->WriteActionTime("Global Execution Time", (unsigned long int) globalTimer.GetElapsedTimeInSeconds());
	tsi->log->Close();
	printf("TSI: Total Execution Time: %lu\n", (unsigned long int) globalTimer.GetElapsedTimeInSeconds());

	tsi->Shutdown();
	//delete DssBCM;
	delete SiBCM;
	//delete DssBAI;
	delete SiBAI;

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
	struct thread_args * args = (struct thread_args *) arg;
	int iter, simul;
	iter = args->iter;
	simul = args->simul;
	Cube *corrCube;

	tsi->log->WriteSeparator();
	tsi->log->WriteIterationNumber(iter);
	tsi->log->WriteSimulationNumber(simul);

	/*
	 * FIXME: Redução do consumo de memória.
	 *     - na DSS o sim/aiCube e o tmp não são necessários ao mesmo tempo,
	 *       podemos _não_ alocar o sim, alocar o tmp, libertar o tmp, e depois alocar o sim.
	 *       (reduzindo efectivamente o consumo de memória em 1 cubo)
	 */
	Cube *simCube = new Cube(tsi->x_values, tsi->y_values, tsi->z_values, tsi->layers->number, DEFAULT_CUBE);

	if(iter == 0) {
		tsi->Dss(simCube, NULL, NULL, 1, iter, simul); // ktype = 1
	} else {
		tsi->Dss(simCube, DssBAI, DssBCM, 5, iter, simul); // ktype = 5
	}

	tsi->Si(simCube, &corrCube, iter, simul);

	/* Compare must run with a mutual-exclusion lock!
	 * because Compare will update the best cubes.
	 * 	there is one alternative that is better for "big-smp" machines:
	 * 		- rw-lock per global-shared BestCube:
	 * 			- read lock for all,
	 * 		 	- aquire write lock when compare must update the global cube.
	 */
	pthread_mutex_lock(&compare_mutex);
	tsi->Compare(simCube, corrCube, &SiBAI, &SiBCM, &bestAiCube,
			&bestGlobalSimulationNumber, iter, simul);
	pthread_mutex_unlock(&compare_mutex);


	/* lets post the finalization of this thread */
	sem_post(&threads_semaphore);

	/* lets free the struct with the thread arguments */
	my_free(arg);

	return NULL;
}

void threading_init_environment()
{
	//tsi->compare_mutex = PTHREAD_MUTEX_INITIALIZER;
	unsigned int threadN = MIN(tsi->ThreadsNumber, tsi->TotalSimulationNumber);
	thread = (pthread_t *) malloc(sizeof(pthread_t) * tsi->TotalSimulationNumber);
	sem_init(&threads_semaphore, 0, threadN);
}
