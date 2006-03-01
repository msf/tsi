#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h> // needed for stat()
#include "cube.h"
#include "tsi-utils.h"
#include "tsi-io.h"
#include "utils.h"
#include "const.h"

void Cube(long x, long y, long z, int layers_n, int type, tCube** this)
{
	int i;
	tCube *aux;
	
	*this = (tCube*) malloc(sizeof(tCube));
	
	(*this)->x_max = x;
	(*this)->y_max = y;
	(*this)->z_max = z;

	(*this)->size = x * y * z;
	(*this)->type = type;
	//(*this)->corrGlobal = -100;
	(*this)->corrAvg = -100;

	(*this)->layers_num = layers_n;


	(*this)->sum_x = NULL;
	(*this)->sum_x2 = NULL;
	(*this)->r = NULL;
	(*this)->layersCorr = NULL;

	/** types of cubes:
	 * 0 - simulation cube, does not alocate space for data, since DSS will create it
	 * 1 - default cube, only alocates space for data (values)
	 * 2 - cube used for the calculus of correlations, must have adicional arrays alocated
	 * 3 - seismic cube, has the arrays for cube of type 2 and also precalculated data 
	 */

	switch ((*this)->type) {
		case SEISMIC_CUBE: 
			(*this)->sum_x = (double*) malloc(((*this)->x_max * (*this)->y_max * (*this)->layers_num) * sizeof(double));
			if ((*this)->sum_x == NULL) {
				fprintf(stderr,"Cube.constructor(): Not enough memory! Error allocating Cube\n");
			}

			(*this)->sum_x2 = (double*) malloc(((*this)->x_max * (*this)->y_max * (*this)->layers_num) * sizeof(double));
			if ((*this)->sum_x2 == NULL) {
				fprintf(stderr,"Cube.constructor(): Not enough memory! Error allocating Cube\n");
			}

		case CORR_CUBE:
			(*this)->r = (double*) malloc(((*this)->x_max * (*this)->y_max * (*this)->layers_num) * sizeof(double));
			if ((*this)->r == NULL) {
				fprintf(stderr,"Cube.constructor(): Not enough memory! Error allocating Cube\n");
			}
			(*this)->layersCorr = (double*) malloc( ((*this)->layers_num) * sizeof(double));
			if ((*this)->layersCorr == NULL) {
				fprintf(stderr,"Cube.constructor(): Not enough memory! Error allocating Cube\n");
			}

			for(i=0; i < (*this)->layers_num; i++)
				(*this)->layersCorr[i] = -100;

		case DEFAULT_CUBE:
			(*this)->values = (float *) malloc( (*this)->size * sizeof(float));
			if ((*this)->values == NULL) {
				fprintf(stderr,"Cube.constructor(): Not enough memory! Error allocating Cube\n");
			}

		case SIM_CUBE:
			break;
		
		default:
			fprintf(stderr,"Cube.constructor(): ERROR, type of cube is unknown: %d\n",(*this)->type);
	}
}

void _Cube(tCube* this)
{
	debug("Cube.Destructor() called\n");
	switch (this->type) {
		case SEISMIC_CUBE:
			debug("Cube.Destructor() seismic cube\n");
			free(this->sum_x);
			free(this->sum_x2);
		case CORR_CUBE:
			debug("Cube.Destructor() corr cube\n");
			free(this->r);
			free(this->layersCorr);
		default:
			debug("Cube.Destructor() default cube\n");
			free(this->values);
	}
}

void Cube_freeGrid(tCube* this)
{
	my_free(this->values);
}

void Cube_allocGrid(tCube* this)
{
	this->values = (float *) malloc( this->size * sizeof(float));
	if (this->values == NULL) {
		fprintf(stderr,"Cube.allocGrid(): Not enough memory! Error allocating Cube\n");
	}
}

void Cube_associateFile(char *str, tCube* this)
{
	this->name = strdup(str);
	printf("associateFile(%s)\n",this->name);
}

/**
 * LoadNewLayers frees the structures related to layers 
 * and realocates them for the new number of layers,
 * since it can vary
 * WARN: this is only used on seismic cubes!
 */
void Cube_LoadNewLayers(unsigned short int layers_n, tCube* this)
{
	int i;
	debug("Cube.LoadNewLayers() called\n");
	
	if(this->r != NULL)
		my_free(this->r);
	if(this->layersCorr != NULL)
		my_free(this->layersCorr);
	if(this->sum_x != NULL)
		my_free(this->sum_x);
	if(this->sum_x2 != NULL)
		my_free(this->sum_x2);

	this->layers_num = layers_n;
	this->r = (double*) malloc((this->x_max * this->y_max * this->layers_num) * sizeof(double));
	if (this->r == NULL) {
		printf("Not enough memory! Error allocating Cube\n");
	}

	this->layersCorr = (double*) malloc(this->layers_num * sizeof(double));
	if (this->layersCorr == NULL) {
		printf("Not enough memory! Error allocating Cube\n");
	}
	for(i=0; i < this->layers_num; i++)
		this->layersCorr[i] = -100;

	this->sum_x = (double*) malloc((this->x_max * this->y_max * this->layers_num) * sizeof(double));
	if (this->sum_x == NULL) {
		fprintf(stderr,"Cube.constructor(): Not enough memory! Error allocating Cube\n");
	}
	this->sum_x2 = (double*) malloc((this->x_max * this->y_max * this->layers_num) * sizeof(double));
	if (this->sum_x2 == NULL) {
		fprintf(stderr,"Cube.constructor(): Not enough memory! Error allocating Cube\n");
	}
}

double * Cube_copyLayersCorr(tCube* this) {
	int n;
	double *corrs;
	int i;
	n = this->layers_num;
	corrs = (double *) malloc ( sizeof(double) * n);
	for(i = 0; i < n; i++) {
		corrs[i] = this->layersCorr[i];
	}
	return corrs;
}

int Cube_ReadFromFileFast(char* file, tCube* this)
{
	FILE *fp;
	long filesize;
	int x, y, z;
	char *file_content;
	struct stat file_stat;


	pdebug2("Opening file for input.");
	printf("Opening file: %s", file);
	// open file for reading
	fp = fopen(file, "r");
	if (fp == NULL) {
		perror("Cube.ReadFromFileFast():Error reading cube file");
		return -1;
	}

	// get file size
	
	if(stat(file, &file_stat)) {
		perror("Error in getting config file size");
		return -1;
	}
	filesize = file_stat.st_size;
   	// read file to memory
	file_content = (char *)malloc(filesize + 1); // alloc space for file contents

	size_t read_bytes;
	read_bytes = fread(file_content, sizeof(char), filesize, fp);	// read full file contents
	file_content[read_bytes-1] = 0;									// terminate string


	// TODO: verificar que n�o existe nada al�m de digitos e pontos e \n

	// parse string
	char seps[] = INTERNAL_CUBE_TOKENS;
	char *token;
	token = strtok(file_content, seps);

	x = y = z = 0;
	for(z=0; z<this->z_max; z++)
		for(y=0; y<this->y_max; y++)
			for(x=0; x<this->x_max; x++) {	
				if (token != NULL) {
					//got value
					// save value to its position in cube
					this->values[Cube_Coords(x, y, z, this)] = strtod(token, NULL);
					
					token = strtok(NULL, seps);
				} else {
					// to quit the 3 cicles
					z = this->z_max + 1;
					y = this->y_max + 1;
					x = this->x_max + 1;
				}
			}

	// close file descriptor
	fclose(fp);
	my_free(file_content); 

	if ((x > this->x_max) && (y > this->y_max) && (z > this->z_max)) {
		fprintf(stderr, "\nError: Uncomplete cube in [%s]. Expecting %ld values.\n", file, this->x_max*this->y_max*this->z_max);
		return -1;
	}

	return 0;
}

int Cube_ReadFromFile_name(char* file, tCube* this)
{
	int x, y, z;
	FILE *fp;
	int ret;

	/*****************************
	double dValue;
	float fValue;

	char ch;
	char buffer[GLOBAL_LINE_SIZE];
	******************************/

	
	// open file for reading
	printf("ReadFromFile_name(%s)\n",file);
	fp = fopen(file, "r");
	if (fp == NULL) {
		perror("Cube.ReadFromFile(): Error reading cube file");
		return -1;
	}

	// TODO: make ReadFromFile() faster.
	for(z=0; z<this->z_max; z++)
		for(y=0; y<this->y_max; y++)
			for(x=0; x<this->x_max; x++) {
				// read value from file
				ret = fscanf(fp, "%f\n", &this->values[Cube_Coords(x, y, z, this)]);
				
				if (ret <= 0) {
					fprintf(stderr, "Cube.ReadFromFile(): ERROR: Read %d, %d, %d (cell %ld)", x, y, z, Cube_Coords(x, y, z,this));
					fprintf(stderr, "Cube.ReadFromFile(): Error reading file (%s) in line %ld.\n", file, Cube_Coords(x,y,z,this));
					return -3;
				}
				
				// EOF found
				if(feof(fp) != 0) {
					// in EOF case, cube is already filled
					if((x+1)*(y+1)*(z+1) < this->x_max * this->y_max * this->z_max) {
						fprintf(stderr, "Cube.ReadFromFile(): EOF found in %s. Uncomplete cube.\nThis cube must have %ld*%ld*%ld = %ld values.", file, this->x_max, this->y_max, this->z_max, this->x_max*this->y_max*this->z_max);
						return -2;
					}
				}
			}
	
	// close file descriptor
	if (fp != NULL)
		fclose(fp);

	return 0;
}

int Cube_SaveToFile_name(char* file,  int root, tCube* this)
{
	int x, y, z;
	FILE *fp;
	int ret_value = 0;
	tCube *mask = NULL;
	double null_value = 0;

	debug("Cube.SaveToFile: Opening file for output.");

	// open file for writing
	fp = fopen(file, "w");
	if (fp == NULL) {
		perror("Cube.SaveToFile(): Error writing file");
		return -1;
	}

	pdebug2("Dumping cube values to file.");

	//header of file
	fprintf(fp,"tsi\n1\ncube\n");

	// TODO: make Cube.SaveToFile() Faster!
	if (root <= 1) {
		if (mask == NULL) {
			// don't use mask
			for(z=0; z<this->z_max; z++)
				for(y=0; y<this->y_max; y++)
					for(x=0; x<this->x_max; x++) {
						// save value from file
						ret_value = fprintf(fp, INTERNAL_DUMP_DOUBLE_VALUE_FORMAT, this->values[Cube_Coords(x,y,z,this)]);
						if (ret_value < 0) {
							perror("Error writing file");
							return ret_value;
						}
					}
		} else {
			// use mask
			for(z=0; z<this->z_max; z++)
				for(y=0; y<this->y_max; y++)
					for(x=0; x<this->x_max; x++) {
						// save value from file
						if (mask->values[Cube_Coords(x,y,z,this)] == 0)
							ret_value = fprintf(fp, INTERNAL_DUMP_DOUBLE_VALUE_FORMAT, null_value);
						else
							ret_value = fprintf(fp, INTERNAL_DUMP_DOUBLE_VALUE_FORMAT, this->values[Cube_Coords(x,y,z,this)]);

						if (ret_value < 0) {
							perror("Error writing file");
							return ret_value;
						}
					}
		}
	} else {
		if (mask == NULL) {
			// don't use mask
			for(z=0; z<this->z_max; z++)
				for(y=0; y<this->y_max; y++)
					for(x=0; x<this->x_max; x++) {
						ret_value = fprintf(fp, INTERNAL_DUMP_DOUBLE_VALUE_FORMAT, nthroot(this->values[Cube_Coords(x,y,z,this)], root));

						if (ret_value < 0) {
							perror("Error writing file");
							return ret_value;
						}
					}
		} else {
			// use mask
			for(z=0; z<this->z_max; z++)
				for(y=0; y<this->y_max; y++)
					for(x=0; x<this->x_max; x++) {
						// save value from file
						if (mask->values[Cube_Coords(x,y,z,this)] == 0)
							ret_value = fprintf(fp, INTERNAL_DUMP_DOUBLE_VALUE_FORMAT, null_value);
						else
							// apply root value
							ret_value = fprintf(fp, INTERNAL_DUMP_DOUBLE_VALUE_FORMAT, nthroot(this->values[Cube_Coords(x, y, z, this)], root));

						if (ret_value < 0) {
							perror("Error writing file");
							return ret_value;
						}
					}
		}
	}

	// close file descriptor
	fclose(fp);

	return 0;
}

int Cube_ReadFromFile(tCube *this)
{
	printf("ReadFromFile(%s)\n",this->name);
	return Cube_ReadFromFile_name(this->name, this);
}

int Cube_SaveToFile(tCube *this)
{
	return Cube_SaveToFile_name(this->name, 1, this);
}

/* NEW CODE!
int Cube::ReadFromFile(char* file)
{
	int i;
	i = readbinary(file, this->values, this->size);
	if( i != 0)
		fprintf(stderr,"ReadFromFile(%s): Error, readbinary returned %d\n",file,i);

	return i;
}

/* SaveToFile.
 * the root argument saves the values rooted by "root", that is,
 *  if root is 1, the value is unmodified
 *  if root is 2, the value saved is "sqrt(value)"
 *
int Cube::SaveToFile(char* file, int root)
{
	float *tmp;
	long i;
	
//	printf("SaveToFile(%s)\n",file);
	if(root <= 1)	
		i = writebinary(file, this->values, this->size);
	else {
		fprintf(stderr,"SaveToFile(): created a new grid to save a new grid with the root values of the original grid\nPossible out of memory condition!\n"); 
		tmp = (float *) malloc( this->size * sizeof(float));
		for(i = 0; i < this->size; i++) {
			tmp[i] = nthroot(this->values[i], root);
		}
		i =	writebinary(file, tmp, this->size);
		free(tmp);
	}

	if(i != 0)
		fprintf(stderr,"SaveToFile(%s): Error, writebinary returned %ld\n",file,i);
	
	return i;
}
*/

void Cube_SetCell(int x, int y, int z, double value, tCube* this)
{
	this->values[Cube_Coords(x,y,z,this)] = value;
}

double Cube_GetCell(int x, int y, int z, tCube* this)
{
	return this->values[Cube_Coords(x,y,z,this)];
}

tCube *Cube_MakeReflectionCoefsCube(tCube* this)
{
	tCube *cr_cube;

	double value;
	int x, y, z;

	Cube(this->x_max,this->y_max,this->z_max,this->layers_num, DEFAULT_CUBE, &cr_cube);

	for(z=0; z<this->z_max-1; z++){
		for(y=0; y<this->y_max; y++){
			for(x=0; x<this->x_max; x++){
				value = (this->values[Cube_Coords(x,y,z+1,this)]-this->values[Cube_Coords(x,y,z,this)])/(this->values[Cube_Coords(x,y,z+1,this)]+this->values[Cube_Coords(x,y,z,this)]);
				cr_cube->values[Cube_Coords(x,y,z,this)] = value;
				
				//printf("MakeReflectionCoefsCube():(ai_cube)this(%d,%d,%d): %.5f\ncr_cube(%d,%d,%d): %.5f\n",x,y,z,this->values[this->Coords(x,y,z)],x,y,z,cr_cube->values[this->Coords(x,y,z)]);
			}
		}
	}

	for(x=0; x<this->x_max; x++)
		for(y=0; y<this->y_max; y++)
			cr_cube->values[Cube_Coords(x,y,this->z_max-1,this)]=0;
	return cr_cube;
}

tCube *Cube_MakeSyntheticCube(tWavelet *wave, tCube* this)
{
	tCube *synth_cube;

	int x, y, z;
	int wavelet_spots;
	int j;
	long it = 0;

	Cube(this->x_max,this->y_max,this->z_max,this->layers_num, DEFAULT_CUBE, &synth_cube);
	
	debug("Cube.MakeSyntheticCube: Calculating CR cube.");
	wavelet_spots = wave->wavelet_used_values / 2;

	for(y=0; y<this->y_max; y++)
		for(x=0; x<this->x_max; x++)
			for(z=0; z<this->z_max; z++)
				synth_cube->values[Cube_Coords(x,y,z,this)] = 0;

	for(y=0; y<this->y_max; y++) {
		for(x=0; x<this->x_max; x++) {
			for(z=0; z<this->z_max; z++) {
				it = z+1;
				for(j=-wavelet_spots+1; j<=wavelet_spots; j++) {
					if ((it+j >= 0) && (it+j < this->z_max-1)) {
						synth_cube->values[Cube_Coords(x,y,it+j,this)] += this->values[Cube_Coords(x,y,z,this)]*Wavelet_PointValue(j,wave);
					}
				}
				//printf("MakeSyntheticCube():synth_cube(%d,%d,%d): %.5f\n",x,y,z,synth_cube->values[this->Coords(x,y,z)]);
			}
		}
	}
	return synth_cube;
}

void Cube_CalculateCorrData(struct layers_t *layers, tCube *this)
{
	int x, y, z;
	int i, z1, z2;

	unsigned long int tmp;

	for(x=0; x < this->x_max; x++) {
		for(y=0; y < this->y_max; y++) {
			z1=0;
			for(i=0; i < layers->number; i++) {
				z2 = z1 + layers->layer[i];

				this->sum_x[Cube_Coords(x,y,i,this)]=0;
				this->sum_x2[Cube_Coords(x,y,i,this)]=0;

				for(z = z1; z < z2; z++) {
						tmp = Cube_Coords(x,y,z,this);
						this->sum_x[Cube_Coords(x,y,i,this)] += this->values[tmp];
						this->sum_x2[Cube_Coords(x,y,i,this)] += this->values[tmp] * this->values[tmp];
				}
				// reposition Z pointer
				z1=z2;
			}
		}
	}
}

void Cube_CalculateStaticCorrSums(tCube* this)
{
	int i;
	double  sum_x = 0, sum_x2 = 0;


	for(i=0; i < this->size; i++) {
		sum_x += this->values[i];
		sum_x2 += this->values[i] * this->values[i];
	}

	this->GlobalSum_x = sum_x;
	this->GlobalSum_x2 = sum_x2;
}

void Cube_Init0(tCube* this)
{
	Cube_Init(0, this);
}

void Cube_Init(double value, tCube* this)
{
	int x, y, z;
	for(x=0; x<this->x_max; x++)
		for(y=0; y<this->y_max; y++)
			for(z=0; z<this->z_max; z++)
				this->values[Cube_Coords(x,y,z, this)] = value;
}

void Cube_GenerateLayersCorr(struct layers_t * layers, tCube* mask, tCube* this)
{
	double acum = 0;
	int counter=0;
	int counter_global=0;
	int x, y;
	int z=0;
	int layer_limit;
	int i;

	if (mask == NULL) { // don't use mask
		for(i=0; i < layers->number; i++) {
			layer_limit = z + layers->layer[i];
			acum = 0;
			counter_global += counter;
			counter=0;
			for(; z < layer_limit; z++)
				for(x=0; x<this->x_max; x++)
					for(y=0; y<this->y_max; y++) {
							acum += this->values[Cube_Coords(x,y,z,this)];
							counter++;
					}
			//printf("Layer %3d: %Lf\n", i, acum / counter); 
			this->layersCorr[i] = acum / counter;
		}
	} else { // use mask
		for(i=0; i < layers->number; i++) {
			layer_limit = z + layers->layer[i];
			acum = 0;
			counter_global += counter;
			counter=0;
			for(; z < layer_limit; z++)
				for(x=0; x<this->x_max; x++)
					for(y=0; y<this->y_max; y++) {
						if (mask->values[Cube_Coords(x,y,z, this)] == 1) {
							acum += this->values[Cube_Coords(x,y,z,this)];
							counter++;
						} 
					}
			//printf("Layer %3d: %Lf\n", i, acum / counter); 
			this->layersCorr[i] = acum / counter;
		}
	}


	/*
	 * code to calculate correlation average from the correlations in each layer 
	 */
	this->corrAvg = 0;
	for(i=0; i<this->layers_num; i++)
		this->corrAvg += this->layersCorr[i];

	this->corrAvg /=  this->layers_num;

	/*
	// DEBUG
	for(z=0; z<1; z++)
		for(y=0; y<this->y_max; y++)
			for(x=0; x<this->x_max; x++)
				printf("value[%d][%d][%d]: %Lf\n", x, y, z, this->values[this->Coords(x,y,z)]);
	*/
}


double Cube_GetLayerCorrAverage(tCube* this)
{
	return this->corrAvg;
}

char* Cube_CorrAverageDump(tCube* this)
{
	//TODO: mudar isto para nao usar tantos arrays (2x1024)
	char buf[1024];
	char ret[1024];
	int i;

	strcpy(ret, "");
	for(i=0; i < this->layers_num; i++) {	
		// printf("%d of %d: %Lf\n", i, this->layers_num, this->cr_med[i]);
		strcpy(buf, "");
		sprintf(buf, "\t\t\tLayer %2d: %.5f\n", i+1, this->layersCorr[i]);
		strcat(ret, buf);
	}
	sprintf(buf, "\t\t\t-------------------\n");
	strcat(ret, buf);
	sprintf(buf, "\t\t\t  Global: %.5f\n", this->globalCorr);
	strcat(ret, buf);
	sprintf(buf, "\t\t\t Average: %.5f\n", this->corrAvg);
	strcat(ret, buf);

	return strdup(ret);
}

void Cube_DumpCRs(tCube* this)
{
	int i;
	for(i=0; i<this->layers_num; i++) {	
		printf("%d of %d: %f\n", i, this->layers_num, this->layersCorr[i]);
	}
}

long Cube_Coords(int x, int y, int z, tCube* this)
{
	return ((z)* this->x_max * this->y_max + (y) * this->x_max + (x+1)) - 1;
}


