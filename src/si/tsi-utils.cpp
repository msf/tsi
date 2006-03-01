#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cube.h"
#include "tsi-utils.h"

/* shouldReadFromFiles(my_sem_t)
 *  caso o count for 0, devolver true.
 *  cc, devolver false
 *  apos o teste e antes de retornar, incrementar o count.
 */
int shouldRestoreFromFiles(my_sem_t *sem)
{

	pthread_mutex_lock(&sem->lock);
	if(sem->count++ == 0) {
		pthread_mutex_unlock(&sem->lock);
		return 1;
	} else {
		pthread_mutex_unlock(&sem->lock);
		return 0;
	}
}

/** shouldWriteToFiles(my_sem_t, int max)
 * decrementar o count.
 * caso o count for 0, devolver true.
 * cc, devolver false.
 */
int shouldSaveToFiles(my_sem_t *sem)
{

	pthread_mutex_lock(&sem->lock);
	sem->count--;
	if(sem->count == 0) {
		pthread_mutex_unlock(&sem->lock);
		return 1;
	} else {
		pthread_mutex_unlock(&sem->lock);
		return 0;
	}
}

void UpdateResultCubes(Cube *BCM, Cube *BAI, Cube *new_corr, Cube *new_ai)
{
	int x, y, z;

//	for(int i = 0; i < BCM->size; i++) {
//		if(BCM->values[i] < new_corr->values[i])
//		{
			/* new best correlation */
//			BCM->values[i] = new_corr->values[i];
//			BAI->values[i] = new_ai->values[i];
//		}
//	}
	
	for(x=0; x<BCM->x_max; x++)
		for(y=0; y<BCM->y_max; y++)
			for(z=0; z<BCM->z_max; z++) {
				if (BCM->values[BCM->Coords(x,y,z)] < new_corr->values[BCM->Coords(x,y,z)]) {
							// melhor correla��o no novo cubo
							// h� que substituir os valores nos cubos BAI e BCM
							BCM->values[BCM->Coords(x,y,z)] = new_corr->values[BCM->Coords(x,y,z)];
							BAI->values[BCM->Coords(x,y,z)] = new_ai->values[BCM->Coords(x,y,z)];
						}
			}
}

/**
 * 	
 * 	calculates correlations by layers between two cubes.
 */
Cube* MakeCorrelationsCube(Cube *seismic, Cube *syntetic, struct layers_t *layers)
{
	debug("MakeCorrelationsCube called\n");
	Cube *corr_cube = new Cube(seismic->x_max, seismic->y_max, seismic->z_max, layers->number, CORR_CUBE);

	int x, y, z;
	int i, z1, z2;
	int temp;

	double sum_xy;
	double sum_x,sum_x_2;
	double sum_y;
	double sum_y2;
	double sum_y_2;
	double denom;

	int n;

	for(x = 0; x < seismic->x_max; x++)
		for(y = 0; y < seismic->y_max; y++) {
			z1 = 0;
			for(i = 0; i < layers->number; i++) {
				z2 = z1+layers->layer[i];

				sum_xy = 0;
				sum_y = 0;
				sum_y2 = 0;
				n = layers->layer[i];

				for(z = z1; z < z2; z++) {
					temp = seismic->Coords(x,y,z);
					
					sum_xy += seismic->values[temp] * syntetic->values[temp];
					sum_y += syntetic->values[temp];
					sum_y2 += syntetic->values[temp] * syntetic->values[temp];

				}

				sum_x = seismic->sum_x[seismic->Coords(x,y,i)];
				sum_x_2 = sum_x * sum_x;
				sum_y_2 = sum_y * sum_y;

				denom = sqrt(
						((n * seismic->sum_x2[seismic->Coords(x,y,i)] ) - sum_x_2 ) * 
						((n * sum_y2 ) - sum_y_2)
						);

				if(denom == 0)
					seismic->r[seismic->Coords(x,y,i)] = 0;
				else
					seismic->r[seismic->Coords(x,y,i)] = ((n * sum_xy) - ( sum_x *sum_y ))/denom;

				if (seismic->r[seismic->Coords(x,y,i)] > 0) {
					for(z = z1; z < z2; z++) {
						corr_cube->values[seismic->Coords(x,y,z)] = seismic->r[seismic->Coords(x,y,i)];
					}
				} else {
					for(z=z1; z < z2; z++) {
						corr_cube->values[seismic->Coords(x,y,z)] = 0;
					}
				}
				
				// reposition Z pointer
				z1=z2;
			}
		}
	/* divide corrAvg for number of points to give the average. */

	return corr_cube;
}
/**
 * Calcular as corre�a��es globais entre dois cubos.
 */ 
double calculateGlobalCorr(Cube *seismic, Cube *syntetic)
{
	debug("calculateGolbalCorr called\n");
	int i;
	int size = syntetic->size;

	double sum_AB = 0;
	double sum_A =  0;
	double sum_A2 = 0;
	double sum_A_2 = 0;
	double sum_B = 0;
	double sum_B2 = 0;
	double sum_B_2 = 0;
	double denom = 0;
	double result;

/*
	sum_A = seismic->GlobalSum_x;
	sum_A2 = seismic->GlobalSum_x2;
*/

	/* FIXME: guardar o sum_X e sum_X2 global do seismic no objecto, e calcular logo na leitura do cubo
	 * para evitar calculos redundantes e acessos � memoria */

  	for(i=0; i < size; i++) {
		sum_AB += seismic->values[i] * syntetic->values[i];
		sum_B += syntetic->values[i];
		sum_B2 += syntetic->values[i] * syntetic->values[i];
		/* values are precomputed
		*/
		sum_A += seismic->values[i];
		sum_A2 += seismic->values[i] * seismic->values[i];
	}
	sum_A_2 = sum_A * sum_A;
	sum_B_2 = sum_B * sum_B;
	
	//printf("sum_SEISMIC = %f\t sum_SY = %f\nsum_SEISMIC2 = %f\t sum_SY2 = %f\nsum_AB = %f\n", sum_A,sum_B,sum_A2,sum_B2,sum_AB);
	denom = sqrt(
			((size * sum_A2) - sum_A_2) *
			((size * sum_B2) - sum_B_2));
	
	if (denom == 0) {
		fprintf(stderr,"calculateGlobalCorr: ERROR: division by zero, returning 0\n");
		return 0;
	} else {
		result =  ((size * sum_AB) - (sum_A * sum_B)) ;	
		result = result/denom;
		return result;
	}
}


char* DumpCorrData(struct aiGridData *aiData)
{
	//TODO: mudar isto para nao usar tantos arrays (2x1024)
	char buf[1024];
	char ret[1024];
	int i;

	strcpy(ret, "");
	for(i=0; i < aiData->layersNum; i++) {	
		// printf("%d of %d: %Lf\n", i, this->layers_num, this->cr_med[i]);
		strcpy(buf, "");
		sprintf(buf, "\t\t\tLayer %2d: %.5f\n", i+1, aiData->layersCorr[i]);
		strcat(ret, buf);
	}
	strcat(ret, "\t\t\t-------------------\n");
	sprintf(buf, "\t\t\t  Global: %.5f\n", aiData->globalCorr);
	strcat(ret, buf);

	return strdup(ret);
}

/**
 * Generator of randomLayers, randomizes:
 * - number of layers
 * - size of each layer
 * also assures minimum_size, minimum_number, and total size
 */
int generateRandomLayers(struct layers_t * l)
{

	unsigned short int min_size = l->minimum_size;
	unsigned short int number = l->minimum_number;
	unsigned short int max_size = l->total_size;

	/*
	printf("generateRandomLayers: min_size: %d, min_number: %d, total_size: %d\n",min_size, number, max_size);
	*/
	int delta = max_size - (min_size * number);
	unsigned short int *lays, *temp;
	unsigned int i, x, sum;
	
	if(delta < 0) { 
		/* too many layers, or layers too big */
		fprintf(stderr,"randomLayers: ERROR, too many layers, or layers too big\n");
		return -1;
	} else if(delta == 0) {
		/* all layers must be of min_size */
		lays = (unsigned short int *) malloc(sizeof(unsigned short int) * number);
		for(i = 0; i < number; i++)
			lays[i] = min_size;
		l->number = number;
		l->layer = lays;
		return 0;
	}
	/* proper random Layers above...
	 * number of layers: number < N < max_size / min_size
	 */
	
	unsigned int n = max_size / min_size;
	temp = (unsigned short int *) calloc(n,sizeof(unsigned short int));

	/* first, we warrantee the minimum number of layers */
	do {
		x = random() % (min_size * 2);
		x += min_size;
		i = (unsigned int) random() % n;

		temp[i] = x;

		x = 0;
		for(i = 0; i < n; i++) {
			if(temp[i] != 0) {
				x++;
			}
		}

	} while(  x < number); 

	/* now we assure the max_size */
	do {
		x = random() % (min_size);
		i = random() % n;
		if(temp[i] == 0)
			temp[i] += min_size;
		temp[i] += x;

		sum = 0;
		for(i = 0; i < n; i++)
			sum += temp[i];
	} while( sum < max_size );

	while(sum > max_size) {
		/* lets find the biggest layer */
		sum = 0;
		for(i = 0; i < n; i++) {
			if(temp[i] > sum) {
				sum = temp[i];
				x = i;
			}
		}
		/* now we reduce its size */
		temp[x] -= min_size;
		temp[x] /= 2;
		temp[x] += min_size;
		/* lets see if its good now */
		sum = 0;
		for(i = 0; i < n; i++)
			sum += temp[i];
	}
	if( sum != max_size) {
		/* add the diference to max_size */
		i = max_size - sum;
		temp[x] += i;
	}
	/* now lets count the layers != 0 */
	sum = 0;
	for(i = 0; i < n; i++) {
		if(temp[i] != 0)
			sum++;
	}

	/* this is the number of layers we are going to return */
	lays = (unsigned short int *) malloc(sizeof(unsigned short int) * sum);

	x = 0;
	for(i = 0; i < n; i++) {
		if(temp[i] != 0)
			lays[x++] = temp[i];
	}
	my_free(temp);

	l->number = sum;
	l->layer = lays;

	return 0;
}


void printLayers(struct layers_t *lay)
{
	for(int i = 0; i < lay->number; i++)
		printf("layer[%d] = %d\n",i,lay->layer[i]);
}

void my_free(void * ptr) 
{
	if(ptr != NULL) {
		free(ptr);
		ptr = NULL;
	}
	else 
		printf("my_free(): ERROR, tried to free null pointer");
}

void debug(char * s) 
{
#ifdef DEBUG
	printf("%s\n",s);
#endif
}
	
