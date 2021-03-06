#ifndef _TSI_IO_H
#define _TSI_IO_H


#include <stdio.h>

#include "log.h"

#define TSI_ASCII_FILE  1 /* this is the GSLIB format */
#define TSI_BIN_FILE    2

#define TSI_FILE	FILE

/* open/close files */
FILE *open_file(char *filename);

FILE *create_file(char *filename);

int close_file(FILE *fp);


/* generic read/write functions */
int read_file(FILE *fp);

int write_file(FILE *fp, char c);

int read_line_file(FILE *fp, char *buf);

int write_line_file(FILE *fp, char *buf);

int read_block_file(FILE *fp, int offset, void *address, unsigned int block_size);

int write_block_file(FILE *fp, int offset, void *address, unsigned int block_size);


/* grid read/write functions */
int read_tsi_grid(FILE *fp, float *address, int x, int y, int z);

int write_tsi_grid(FILE *fp, float *address, int x, int y, int z);

int read_cartesian_grid(FILE *fp, float *grid, unsigned int grid_size);

int write_cartesian_grid(FILE *fp, float *grid, unsigned int grid_size);

int read_gslib_header(log_t *l, FILE *fp, int fields_per_line);

int read_gslib_grid(log_t *l, FILE *fp, float *grid, unsigned int grid_size);

int write_gslib_grid(FILE *fp, float *grid, int x, int y, int z, char *desc);


/* read/write floats in binary format */
unsigned int read_float(FILE *fp, float *grid, unsigned int nelems);

unsigned int write_float(FILE *fp, float *grid, unsigned int nelems);

int dump_binary_grid(FILE *fp, float *grid, unsigned int grid_size);

int load_binary_grid(FILE *fp, float *grid, unsigned int grid_size);

#endif /* _TSI_IO_H */
