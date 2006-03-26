#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "tsi_io.h"
#ifndef TSI_MPI

/***************************************************
 *
 *  Binary grid manipulation functions
 *
 ***************************************************/

TSI_FILE *open_grid_file(char *filename) {
    TSI_FILE *fp;
    fp = fopen(filename, "w+b");
    return fp;
} /* open_grid_file */



int close_grid_file(TSI_FILE *fp) {
    return fclose(fp);
} /* close_grid_file */



int write_grid_file(TSI_FILE *fp, float *address, unsigned int size) {
    if (fseek(fp, (long) 0, SEEK_SET)) return 0;
    if (fwrite(address, sizeof(float), (size_t) size, fp) < size) return 0;
    return 1;
} /* write_grid_file */



int read_grid_file(TSI_FILE *fp, float *address, unsigned int size) {
    if (fseek(fp, (long) 0, SEEK_SET)) return 0;
    if (fread(address, sizeof(float), (size_t) size, fp) < size) return 0;
    return 1;
} /* load_grid_file */



/***************************************************
 *
 *  Byte files manipulation functions
 *
 ***************************************************/

TSI_FILE *open_file(char *filename) {
    TSI_FILE *fp;
    fp = fopen(filename, "r+b");
    return fp;
} /* open_file */



TSI_FILE *create_file(char *filename) {
    return open_grid_file(filename);
} /* create_file */



int close_file(TSI_FILE *fp) {
    return close_grid_file(fp);
} /* close_file */



int read_file(TSI_FILE *fp) {
    return getc(fp);
} /* read_file */



int read_line_file(TSI_FILE *fp, char *buf) {
    int c;

    do {
        c = read_file(fp);
        if (c == EOF) break;
        *buf = (char) c;
        buf++;
    } while (c != '\n');
    *buf = 0;
    return c;
} /* read_line_file */



int read_block_file(TSI_FILE *fp, int offset, void *address, unsigned int block_size) {
    if (fseek(fp, offset, SEEK_SET)) return 0;
    if (fread(address, 1, block_size, fp) < (unsigned int) block_size) return 0;
    return offset;
} /* read_block_file */



int write_file(TSI_FILE *fp, char c) { /* TEST */
    return putc(c, fp);
} /* write_file */



int write_line_file(TSI_FILE *fp, char *buf) { /* TEST */
    int ret;

    ret = 0;
    while (*buf) {
        ret = write_file(fp, *buf);
        buf++;
        if (ret < 0) break;
    }
    return ret;
} /* write_line_file */



int write_block_file(TSI_FILE *fp, int offset, void *address, unsigned int block_size) {
    if (fseek(fp, offset, SEEK_SET)) return 0;
    if (fwrite(address, 1, block_size, fp) < (unsigned int) block_size) return 0;
    return offset;
} /* write_block_file */



/***************************************************
 *
 *  ASCII grid manipulation functions
 *
 ***************************************************/

int read_ascii_grid_file(TSI_FILE *fp, float *g, unsigned int grid_size) {
    unsigned int i;
	
	for(i=0; i < grid_size; i++)
		fscanf(fp, "%f\n", &g[i]);

	return i;
}

int write_ascii_grid_file(TSI_FILE *fp, float *g, unsigned int grid_size) { /* TEST */
    unsigned int i;

	fprintf(fp,"tsi\n1\ngrid\n"); /* ugly hack for SGEMS */

	for(i=0; i < grid_size; i++)
		fprintf(fp, "%.3f\n", g[i]);

	return i;
} /* write_ascii_grid_file */


#endif /* TSI_MPI */

/* end of tsi_io.c */
