#ifndef _TSI_UTILS_H
#define _TSI_UTILS_H

#include <pthread.h> /* pthread_mutex_t */
#include "cube.h"

/* 
 * each layer has a size around [20-300]
 * and for cubes/grids of absurd sizes, 65535 layers will do
 */
struct layers_t {
	unsigned short int minimum_size;
	unsigned short int minimum_number;
	unsigned int	total_size;
	unsigned short int number;
	unsigned short int *layer;
};

struct my_sem {
	pthread_mutex_t lock;
	short int count;
};	
typedef struct my_sem my_sem_t;

/* information of ai grid used on comparisons */
struct  aiGridData {
	double globalCorr;
	double *layersCorr;
	int layersNum;
};



/* functions in tsi-utils.c */
extern int shouldRestoreFromFiles(my_sem_t *);
extern int shouldSaveToFiles(my_sem_t *);
extern double calculateGlobalCorr(tCube *, tCube *);
extern void  UpdateResultCubes(tCube *, tCube *, tCube *, tCube *);
extern tCube * MakeCorrelationsCube(tCube *, tCube *, struct layers_t *);
extern char* DumpCorrData(struct aiGridData *aiData);
extern int  generateRandomLayers(struct layers_t *);
extern void printLayers(struct layers_t *);
extern void debug(char *s);
extern void my_free(void *);

#endif 

