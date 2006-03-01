#ifndef _TSI_UTILS_H
#define _TSI_UTILS_H

#include <pthread.h> /* pthread_mutex_t */

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
int shouldRestoreFromFiles(my_sem_t *);
int shouldSaveToFiles(my_sem_t *);
double calculateGlobalCorr(Cube *, Cube *);
void  UpdateResultCubes(Cube *, Cube *, Cube *, Cube *);
Cube * MakeCorrelationsCube(Cube *, Cube *, struct layers_t *);
char* DumpCorrData(struct aiGridData *aiData);
int  generateRandomLayers(struct layers_t *);
void printLayers(struct layers_t *);
void debug(char *s);
void my_free(void *);

#endif 
