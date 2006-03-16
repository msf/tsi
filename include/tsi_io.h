#ifndef _TSI_IO_H
#define _TSI_IO_H

#ifdef TSI_MPI
#include <mpi.h>
#define TSI_FILE MPI_File
#ifndef EOF
#define EOF (-1)
#endif /* EOF */
#else
#include <stdio.h>
#define TSI_FILE FILE
#endif /* TSI_MPI */


TSI_FILE *open_grid_file(char *filename);

int close_grid_file(TSI_FILE *fp);

int write_grid_file(TSI_FILE *fp, float *address, unsigned int size);

int read_grid_file(TSI_FILE *fp, float *address, unsigned int size);

TSI_FILE *open_file(char *filename);

TSI_FILE *create_file(char *filename);

int close_file(TSI_FILE *fp);

int read_file(TSI_FILE *fp);

int read_line_file(TSI_FILE *fp, char *buf);

int read_block_file(TSI_FILE *fp, int offset, void *address, int block_size);

int write_file(TSI_FILE *fp, char c);

int write_line_file(TSI_FILE *fp, char *buf);

int write_block_file(TSI_FILE *fp, int offset, void *address, int block_size);

int write_ascii_grid_file(TSI_FILE *fp, float *address, unsigned int size);

int read_ascii_grid_file(TSI_FILE *fp, float *address, unsigned int size);


#endif /* _TSI_IO_H */
