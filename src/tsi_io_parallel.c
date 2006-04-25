#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "memdebug.h"
#include "tsi_io.h"
#ifdef TSI_MPIIO

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


/**** open/close files ****/

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
} /* open_grid_file */



TSI_FILE *create_file(char *filename) {
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
} /* create_file */



int close_file(TSI_FILE *fp) {
    if (MPI_File_close(fp) != MPI_SUCCESS) {  /* MPI-2 */
        tsi_free(fp);
        return 0;
    }
    tsi_free(fp);
    return 1;
} /* close_grid_file */



/**** generic read/write functions ****/

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



int write_file(TSI_FILE *fp, char c) { /* TEST */
    MPI_Status status;
    int n = 0;

    if (MPI_File_write(*fp, &c, 1, MPI_CHAR, &status) != MPI_SUCCESS) {  /* MPI-2 */
        return EOF;
    }
    return 0;
} /* write_file */



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



int read_block_file(TSI_FILE *fp, int offset, void *address, unsigned int block_size) {
    MPI_Status status;
    if (MPI_File_read_at(*fp, offset, address, block_size, MPI_BYTE, &status) != MPI_SUCCESS) return 0;  /* MPI-2 */
    return offset;
} /* read_block_file */



int write_block_file(TSI_FILE *fp, int offset, void *address, unsigned int block_size) {
    MPI_Status status;
    if (MPI_File_write_at(*fp, offset, address, block_size, MPI_BYTE, &status) != MPI_SUCCESS) return 0;  /* MPI-2 */
    return offset;
} /* write_block_file */



/**** grid read/write functions ****/

int read_tsi_grid(TSI_FILE *fp, float *grid, int x, int y, int z) {
    char header[64];
    char tsi_h[] = "TSI";
    int err, type, x1, y1, z1, nl, i;
    unsigned int grid_size;

    /* load file header */
    if (read_line_file(fp, header) < 0) return 0;

    /* parse header */
    if (strncmp(header, tsi_h, 3)) {
        printf_dbg("\tread_tsi_grid(): unknown file format\n");
        return 0;
    }

    sscanf(header, "TSI %d %d %d %d\n", &type, &x1, &y1, &z1);
    if ((type != 1) && (type != 2)) {
        printf_dbg("\tread_tsi_grid(): incompatible format\n");
        return 0;
    }
    if ((x != x1) || (y != y1) || (z != z1)) {
        printf_dbg("\tread_tsi_grid(): incoeherent size parameters\n");
        return 0;
    }
    grid_size = (unsigned int) x * (unsigned int) y * (unsigned int) z;
    switch (type) {
        case 1:      /* ASCII data */
            read_line_file(fp, header);
            read_line_file(fp, header);
            return read_cartesian_grid(fp, grid, grid_size);

        case 2:      /* data in binary floats */
            return read_float(fp, grid, grid_size);

        default:
            printf_dbg("\tread_tsi_grid(): unknown file format\n");
            break;
    } /* switch */

    return 0;
} /* read_tsi_grid */



int write_tsi_grid(TSI_FILE *fp, int type, float *grid, int x, int y, int z) {
    char header[64];

    if (type == 1)  /* ASCII file */
        return write_gslib_grid(fp, grid, x, y, z, NULL);

    /* bin file */
    sprintf(header, "TSI 2 %d %d %d\n", x, y, z);
    if (write_line_file(fp, header) < 0) return 0;
    
    return write_float(fp, grid, (unsigned int) x * (unsigned int) y * (unsigned int) z);
} /* write_tsi_grid */



//int read_tsi_cmgrid(TSI_FILE *fp, cm_grid *cmg, int x, int y, int z) {
//    return 0;
//} /* read_tsi_cmgrid */



//int write_tsi_cmgrid(TSI_FILE *fp, cm_grid *cmg, float *address, int x, int y, int z) {
//    return 0;
//} /* read_cartesian_grid */



int read_cartesian_grid(TSI_FILE *fp, float *g, unsigned int grid_size) {
    unsigned int i;
    char val[64];
    int ret;

    for (i = 0; i < grid_size; i++) {
        ret = read_line_file(fp, val);
        if (ret == EOF) break;
        sscanf(val, "%f\n", &g[i]);
    }
    if (i < grid_size) return 0;
    return i;
}

int write_cartesian_grid(TSI_FILE *fp, float *g, unsigned int grid_size) { /* TEST */
    unsigned int i;
    char val[64];
    int ret;

    for (i = 0; i < grid_size; i++) {
        sprintf(val, "%.4f\n", g[i]);
        ret = write_line_file(fp, val);
        if (ret < 0) break;
    }
    if (i < grid_size) return 0;
    return i;
} /* write_ascii_grid_file */



int read_gslib_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    unsigned int i;
    char str[64];

    /* ignore header */
    if (read_line_file(fp, str) < 0) return 0;
    if (read_line_file(fp, str) < 0) return 0;
    if (read_line_file(fp, str) < 0) return 0;

    return read_cartesian_grid(fp, grid, grid_size);
} /* read_gslib_grid */



int write_gslib_grid(TSI_FILE *fp, float *grid, int x, int y, int z, char *desc) {
    unsigned int i, grid_size;
    char header[128];

    grid_size = (unsigned int) x * (unsigned int) y * (unsigned int) z;
    sprintf(header, "TSI 1 %d %d %d\n", x, y, z);
    if (write_line_file(fp, header) < 0) return 0;
    sprintf(header, "1\n");
    if (write_line_file(fp, header) < 0) return 0;
    if (desc)
        sprintf(header, "%s\n", desc);
    else
        sprintf(header, "grid\n");
    if (write_line_file(fp, header) < 0) return 0;

    return write_cartesian_grid(fp, grid, grid_size);
} /* write_gslib_grid */




/********************** EVAL *****************************/


int dump_binary_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    MPI_Status status;
    if (MPI_File_write_at(*fp, 0, grid, grid_size, MPI_FLOAT, &status) != MPI_SUCCESS) return 0;  /* MPI-2 */
    return 1;
} /* write_grid_file */



int load_binary_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    MPI_Status status;
    if (MPI_File_read_at(*fp, 0, grid, grid_size, MPI_FLOAT, &status) != MPI_SUCCESS) return 0;  /* MPI-2 */
    return 1;
} /* load_grid_file */



int read_float(TSI_FILE *fp, float *grid, unsigned int nelems)
{
    MPI_Status status;
    int n = 0;

    if (MPI_File_read(*fp, grid, nelems, MPI_FLOAT, &status) == MPI_SUCCESS)  /* MPI-2 */
        if (MPI_Get_count(&status, MPI_FLOAT, &n) == MPI_SUCCESS)
            if (n != MPI_UNDEFINED)
                return n;
    return 0;
} /* read_float */



int write_float(TSI_FILE *fp, float *grid, unsigned int nelems)
{
    MPI_Status status;
    int n = 0;

    if (MPI_File_write(*fp, grid, nelems, MPI_FLOAT, &status) != MPI_SUCCESS) {  /* MPI-2 */
        return 0;
    }
    return 1;
} /* write_float */



#endif /* TSI_MPIIO */

/* end of tsi_io_parallel.c */

