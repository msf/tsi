#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


int readasciiXYZ(char *str, float *values, int xsiz, int ysiz, int zsiz) 
{
	int x,y,z;
	int zstep, ystep, pos;
	int xymax = xsiz * ysiz;
	FILE *fp = NULL;

	fp = fopen(str, "r");
	if (fp == NULL) {
		perror("Cube.ReadFromFile(): Error reading cube file");
		return -1;
	}

	for(z = 0; z < zsiz; z++) {
		zstep = xymax * z;
		for(y = 0; y < ysiz; y++) {
			ystep = xsiz * y;
			for(x = 0; x < xsiz; x++) {
				pos = zstep + ystep + (x+1) -1;
				fscanf(fp, "%f\n", &values[pos]);
			}
		}
	}

	// EOF found
	if(feof(fp) ) {
		perror("readAscii: Error reading cube file");
	}

	fclose(fp);
	return 0;
}

int readascii(char *str, float *values, unsigned int size) 
{

	unsigned int i;
	FILE *fp = NULL;

	fp = fopen(str, "r");
	if (fp == NULL) {
		perror("readAscii: Error reading cube file");
		return -1;
	}

	for(i=0; i < size; i++)
		fscanf(fp, "%f\n", &values[i]);


	fclose(fp);
	return 0;
}

int readbinary(char *str, float *values, long nelems) 
{

	int fd;
	long int err;
	long size = sizeof(float) * nelems;

	fd = open(str, O_RDONLY);
	if (fd < 0) {		
		perror("readBinary: Error reading cube file");
		return -1;
	}

	err = read(fd, values, size);
	if( err < size) {
			printf("readBinary: file ended too soon\n");
	}	

	close(fd);
	return 0;
}

int writeasciiXYZ(char *str, float *values, int xsiz, int ysiz, int zsiz) 
{
	int x,y,z;
	int zstep, ystep, pos;
	int xymax = xsiz * ysiz;
	FILE *fp;

	fp = fopen(str, "w");
	if (fp == NULL) {
		perror("writeAscii: Error writing cube file");
	}

	fprintf(fp,"tsi\n1\ngrid\n");

	for(z=0; z < zsiz; z++) {
		zstep = xymax * z;
		for(y=0; y< ysiz; y++) {
			ystep = xsiz * y;
			for(x=0; x< xsiz; x++) {
				pos = zstep + ystep + (x+1) -1;
				fprintf(fp, "%.3f\n", values[pos]);
			}
		}
	}

	fclose(fp);
	return 0;
}

int writeascii(char *str, float *values, unsigned long size)
{

	FILE *fp;
	unsigned int i;

	fp = fopen(str, "w");
	if (fp == NULL) {
		perror("writeAscii: Error opening file for write; ");
	}

	fprintf(fp,"tsi\n1\ngrid\n");

	for(i=0; i < size; i++)
		fprintf(fp, "%.3f\n", values[i]);

	fclose(fp);
	return 0;
}

int writebinary(char *str, float *values, long nelems) 
{

	int fd;
	long int err;
   	long	size = sizeof(float) * nelems;

	fd = open(str, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
	if (fd < 0) {
		perror("writeBinary: Error opening file for write; ");
		return -1;
	}

	err = write(fd, values, size );
	if( err < size ) {
			fprintf(stderr,"writeBinary: error writing file %s\n write returned: %ld\n",str,err);
	}	

	close(fd);
	return 0;
}

