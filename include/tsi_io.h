#ifndef _TSI_IO_H
#define _TSI_IO_H


#include <stdio.h>
#define TSI_FILE FILE


//#include "si.h"      <- fuck-up

/* open/close files */
TSI_FILE *open_file(char *filename);

TSI_FILE *create_file(char *filename);

int close_file(TSI_FILE *fp);


/* generic read/write functions */
int read_file(TSI_FILE *fp);

int write_file(TSI_FILE *fp, char c);

int read_line_file(TSI_FILE *fp, char *buf);

int write_line_file(TSI_FILE *fp, char *buf);

int read_block_file(TSI_FILE *fp, int offset, void *address, unsigned int block_size);

int write_block_file(TSI_FILE *fp, int offset, void *address, unsigned int block_size);


/* grid read/write functions */
int read_tsi_grid(TSI_FILE *fp, float *address, int x, int y, int z);

int write_tsi_grid(TSI_FILE *fp, int type, float *address, int x, int y, int z);

//int read_tsi_cmgrid(TSI_FILE *fp, cm_grid *cmg, int x, int y, int z);

//int write_tsi_cmgrid(TSI_FILE *fp, cm_grid *cmg, int x, int y);

int read_cartesian_grid(TSI_FILE *fp, float *grid, unsigned int grid_size);

int write_cartesian_grid(TSI_FILE *fp, float *grid, unsigned int grid_size);

int read_gslib_grid(TSI_FILE *fp, float *grid, unsigned int grid_size);

int write_gslib_grid(TSI_FILE *fp, float *grid, int x, int y, int z, char *desc);


/* read/write floats in binary format */
int read_float(TSI_FILE *fp, float *grid, unsigned int nelems);

int write_float(TSI_FILE *fp, float *grid, unsigned int nelems);

int dump_binary_grid(TSI_FILE *fp, float *grid, unsigned int grid_size);

int load_binary_grid(TSI_FILE *fp, float *grid, unsigned int grid_size);

#endif /* _TSI_IO_H */
