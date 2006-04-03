#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "memdebug.h"
#include "tsi_io.h"
#ifdef TSI_MPI

/*
 *
 *  TODO: 
 *  - improve error handling
 *  - non-blocking i/o for heap speed-up
 *    - read_grid_nowait
 *    - write_grid_nowait
 *    - test_grid_req (evals read/write completion)
 *
 */


/***************************************************
 *
 *  Binary grid manipulation functions
 *
 ***************************************************/

TSI_FILE *open_grid_file(char *filename) {
    TSI_FILE *fp;
    int amode;
    MPI_Info info;
    
    fp = tsi_malloc(sizeof(MPI_File));
    info = MPI_INFO_NULL;
    amode = MPI_MODE_RDWR | MPI_MODE_CREATE;
    if (MPI_File_open(MPI_COMM_WORLD, filename, amode, info, fp) != MPI_SUCCESS) {  /* MPI-2 */
        tsi_free(fp);
        return NULL;
    }
    return fp;
} /* open_grid_file */



int close_grid_file(TSI_FILE *fp) {
    if (MPI_File_close(fp) != MPI_SUCCESS) {  /* MPI-2 */
        tsi_free(fp);
        return 0;
    }
    tsi_free(fp);
    return 1;
} /* close_grid_file */



int write_grid_file(TSI_FILE *fp, float *address, unsigned int size) {
    MPI_Status status;
    if (MPI_File_write_at(*fp, 0, address, size, MPI_FLOAT, &status) != MPI_SUCCESS) return 0;  /* MPI-2 */
    return 1;
} /* write_grid_file */



int read_grid_file(TSI_FILE *fp, float *address, unsigned int size) {
    MPI_Status status;
    if (MPI_File_read_at(*fp, 0, address, size, MPI_FLOAT, &status) != MPI_SUCCESS) return 0;  /* MPI-2 */
    return 1;
} /* load_grid_file */



/***************************************************
 *
 *  Byte files manipulation functions
 *
 ***************************************************/

TSI_FILE *open_file(char *filename) {
    TSI_FILE *fp;
    int amode;
    MPI_Info info;
    
    fp = tsi_malloc(sizeof(MPI_File));
    info = MPI_INFO_NULL;
    amode = MPI_MODE_RDWR;
    if (MPI_File_open(MPI_COMM_WORLD, filename, amode, info, fp) != MPI_SUCCESS) {  /* MPI-2 */
        tsi_free(fp);
        return NULL;
    }
    return fp;
} /* open_file */



TSI_FILE *create_file(char *filename) {
    return open_grid_file(filename);
} /* create_file */



int close_file(TSI_FILE *fp) {
    return close_grid_file(fp);
} /* close_file */



int read_file(TSI_FILE *fp) {
    MPI_Status status;
    char c = 0;
    int n = 0;

    if (MPI_File_read(*fp, &c, 1, MPI_CHAR, &status) == MPI_SUCCESS) {  /* MPI-2 */
        if (MPI_Get_count(&status, MPI_CHAR, &n) == MPI_SUCCESS)
            if (n > 0)
                return (int) c;
            else
                return EOF;
    }
    return 0;
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
    MPI_Status status;
    if (MPI_File_read_at(*fp, offset, address, block_size, MPI_BYTE, &status) != MPI_SUCCESS) return 0;  /* MPI-2 */
    return offset;
} /* read_block_file */



int write_file(TSI_FILE *fp, char c) { /* TEST */
    MPI_Status status;
    int n = 0;

    if (MPI_File_write(*fp, &c, 1, MPI_CHAR, &status) != MPI_SUCCESS) {  /* MPI-2 */
        return EOF;
    }
    return 0;
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
    MPI_Status status;
    if (MPI_File_write_at(*fp, offset, address, block_size, MPI_BYTE, &status) != MPI_SUCCESS) return 0;  /* MPI-2 */
    return offset;
} /* write_block_file */



/***************************************************
 *
 *  ASCII grid manipulation functions
 *
 ***************************************************/

int read_ascii_grid_file(TSI_FILE *fp, float *g, unsigned int grid_size) {
    unsigned int i;
    char val[64];
    int ret;

    for (i = 0; i < grid_size; i++) {
        ret = read_line_file(fp, val);
        if (ret == EOF) return ret;
        sscanf(val, "%f\n", g+i);
    }
    return i;
}

int write_ascii_grid_file(TSI_FILE *fp, float *g, unsigned int grid_size) { /* TEST */
    unsigned int i;
    char val[64];
    int ret;

    strcpy(val, "tsi\n1\ngrid\n");  /* ugly hack for SGEMS */
    ret = write_line_file(fp, val);
    if (ret < 0) return ret;
    for (i = 0; i < grid_size; i++) {
        sprintf(val, "%.3f\n", g[i]);
        ret = write_line_file(fp, val);
        if (ret < 0) return ret;
    }
    return i;
} /* write_ascii_grid_file */


#endif /* TSI_MPI */

/* end of tsi_io.c */
