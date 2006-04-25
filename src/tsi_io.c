#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "tsi.h"
#include "tsi_io.h"

#ifndef TSI_MPIIO

/* I/O prototypes */
/* (none) */

/**** open/close files ****/

TSI_FILE *open_file(char *filename) {
    TSI_FILE *fp;
    fp = fopen(filename, "r+b");
    return fp;
} /* open_file */



TSI_FILE *create_file(char *filename) {
    TSI_FILE *fp;
    fp = fopen(filename, "w+b");
    return fp;
} /* create_file */



int close_file(TSI_FILE *fp) {
    return fclose(fp);
} /* close_file */


/**** grid read/write functions ****/

int read_tsi_grid(TSI_FILE *fp, float *grid, int x, int y, int z) {
    char header[64];
    char tsi_h[] = "TSI";
    int err, type, x1, y1, z1, nl, i;
    unsigned int grid_size;

    /* load file header */
	
	if ( fread(header, sizeof(char), 3, fp ) < 3)
			return 0;


    /* parse header */
    if (strncmp(header, tsi_h, 3)) {
        printf_dbg("\tread_tsi_grid(): unknown file format\n");
        return 0;
    }

    fscanf(fp, "%d %d %d %d\n", &type, &x1, &y1, &z1);
    if ((type != 1) && (type != 2)) {
        printf_dbg("\tread_tsi_grid(): incompatible format\n");
        return 0;
    }
    if ((x != x1) || (y != y1) || (z < z1)) {
        printf_dbg("\tread_tsi_grid(): incoeherent size parameters\n");
        return 0;
    }
    grid_size = (unsigned int) x * (unsigned int) y * (unsigned int) z;
    switch (type) {
        case 1:      /* ASCII data */
			fgets( header, 63, fp);
			fgets( header, 63, fp);
			fgets( header, 63, fp);
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
    if (type == 1)  /* ASCII file */
        return write_gslib_grid(fp, grid, x, y, z, NULL);

    /* bin file */
    if (!fprintf(fp,"TSI 2 %d %d %d\n", x, y, z))
        return 0;

    return write_float(fp, grid, (unsigned int) x * (unsigned int) y * (unsigned int) z);
} /* write_tsi_grid */



//int read_tsi_cmgrid(TSI_FILE *fp, cm_grid *cmg, int x, int y, int z) {
//    return 0;
//} /* read_tsi_cmgrid */



//int write_tsi_cmgrid(TSI_FILE *fp, cm_grid *cmg, float *address, int x, int y, int z) {
//    return 0;
//} /* read_cartesian_grid */



int read_cartesian_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    unsigned int i;
	
    for (i=0; i < grid_size; i++) 
		fscanf(fp, "%f\n", &grid[i]);
    return i;
} /* read_cartesian_grid */



int write_cartesian_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    unsigned int i;

    for (i = 0; i < grid_size; i++) 
		fprintf(fp, "%.3f\n", grid[i]);
    return i;
} /* write_cartesian_grid */



int read_gslib_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    unsigned int i;
    char str[64];

    /* ignore header */
    if (fgets(str, 64, fp) < 0) return 0;
    if (fgets(str, 64, fp) < 0) return 0;
    if (fgets(str, 64, fp) < 0) return 0;

    return read_cartesian_grid(fp, grid, grid_size);
} /* read_gslib_grid */



int write_gslib_grid(TSI_FILE *fp, float *grid, int x, int y, int z, char *desc) {
    unsigned int i, grid_size;

    grid_size = (unsigned int) x * (unsigned int) y * (unsigned int) z;
    if (desc) {
        if (!fprintf(fp,"TSI 1 %d %d %d\n1\n%s\n", x, y, z, desc))
            return 0; 
    } else {
        if (!fprintf(fp,"TSI 1 %d %d %d\n1\ngrid\n", x, y, z))
            return 0;
    }

    return write_cartesian_grid(fp, grid, grid_size);
} /* write_gslib_grid */



/************************* EVAL *********************************************/

int dump_binary_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    if (fseek(fp, (long) 0, SEEK_SET)) 
		return 0;
	if ( write_float(fp, grid, grid_size) < grid_size)
		return 0;

    return 1;
} /* write_grid_file */

int load_binary_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    if (fseek(fp, (long) 0, SEEK_SET)) return 0;
	if ( read_float(fp, grid, grid_size) < grid_size)
		return 0;

    return 1;
} /* load_grid_file */



int read_float(TSI_FILE *fp, float *grid, unsigned int nelems)
{
    int err;

    err = fread(grid, sizeof(float), nelems, fp);
    if (err < nelems) {
	fprintf(stderr,"read_binary: fread returned %d\n", err);
    }
    return err;
} /* read_float */



int write_float(TSI_FILE *fp, float *grid, unsigned int nelems)
{
	int err;

	err = fwrite(grid, sizeof(float), nelems, fp);
	if (err < nelems) {
		fprintf(stderr,"write_binary: fwrite returned %d\n",err);
	}

	return err;
} /* write_float */

#endif /* TSI_MPIIO */



/* end of tsi_io.c */

