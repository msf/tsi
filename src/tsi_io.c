#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "tsi.h"
#include "tsi_io.h"

/**
 * high level functions
 */

int tsi_read_grid(tsi *t, TSI_FILE *fp, float *grid, int type) {
	switch(type) {
		case CARTESIAN_FILE:
			return (read_cartesian_grid(fp, grid, t->grid_size));
		case TSI_ASCII_FILE:
		case TSI_BIN_FILE:
			return (read_tsi_grid(fp, grid, t->xsize, t->ysize, t->zsize));
		case GSLIB_FILE:
			return (read_gslib_grid(fp, grid, t->grid_size));
			/*
			   case SGEMS_FILE:
			   return read_sgems_grid(fp, grid, t->grid_size));
			 */
		default:
			fprintf(stderr, "ERROR: Unknown grid file type!\n");
			return 0;
	} /* switch */
} /* tsi_read_grid */

int tsi_write_grid(tsi *t, TSI_FILE *fp, float *grid, int type, char *desc) {
	switch(type) {
		case CARTESIAN_FILE:
			return (write_cartesian_grid(fp, grid, t->grid_size));
		case TSI_ASCII_FILE:
		case TSI_BIN_FILE:
			return (write_tsi_grid(fp, type, grid, t->xsize, t->ysize, t->zsize));
		case GSLIB_FILE:
			return (write_gslib_grid(fp, grid, t->xsize, t->ysize, t->zsize, desc));
			/*
			   case SGEMS_FILE:
			   return write_sgems_grid(fp, grid, t->grid_size));
			 */
		default:
			fprintf(stderr, "ERROR: Unknown grid file type!\n");
			return 0;
	} /* switch */ 
} /* tsi_write_grid */

/**** 
 * open/close files ***
 */

TSI_FILE *open_file(char *filename) {
    TSI_FILE *fp;
    fp = fopen(filename, "r+");
    return fp;
} /* open_file */


TSI_FILE *create_file(char *filename) {
    TSI_FILE *fp;
    fp = fopen(filename, "w+");
    return fp;
} /* create_file */


int close_file(TSI_FILE *fp) {
    return fclose(fp);
} /* close_file */


/**** 
 * grid read/write functions ***
 */

int read_tsi_grid(TSI_FILE *fp, float *grid, int x, int y, int z) {
    char header[64];
    char tsi_h[] = "TSI";
    int type, x1, y1, z1;
    unsigned int grid_size;

    /* load file header */
	
    if ( fread(header, sizeof(char), 3, fp ) < 3) {
	fprintf(stderr,"\tread_tsi_grid(): unknown file format\n");
	return 0;
    }


    /* parse header */
    if (strncmp(header, tsi_h, 3)) {
        fprintf(stderr,"\tread_tsi_grid(): unknown file format\n");
        return 0;
    }

    fscanf(fp, "%d %d %d %d\n", &type, &x1, &y1, &z1);
    if ((type != 1) && (type != 2)) {
        fprintf(stderr,"\tread_tsi_grid(): incompatible format\n");
        return 0;
    }
    if ((x != x1) || (y != y1) || (z < z1)) {
        fprintf(stderr,"\tread_tsi_grid(): incoeherent size parameters\n");
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
            fprintf(stderr,"\tread_tsi_grid(): unknown file format\n");
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


int read_cartesian_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    unsigned int i;
	
    for (i=0; i < grid_size; i++) 
	    fscanf(fp, "%f\n", &grid[i]);
    return i;
} /* read_cartesian_grid */



int write_cartesian_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    unsigned int i;

    for (i = 0; i < grid_size; i++) 
	    fprintf(fp, "%.4f\n", grid[i]);
    return i;
} /* write_cartesian_grid */



int read_gslib_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    char str[65];

    /* ignore header */
    if (fgets(str, 64, fp) == NULL) return 0;
    if (fgets(str, 64, fp) == NULL) return 0;
    if (fgets(str, 64, fp) == NULL) return 0;

    return read_cartesian_grid(fp, grid, grid_size);
} /* read_gslib_grid */



int write_gslib_grid(TSI_FILE *fp, float *grid, int x, int y, int z, char *desc) {
    unsigned int grid_size;

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
    printf_dbg2("\tdump_binary_grid(): %d\n", grid_size);
    if (fseek(fp, (long) 0, SEEK_SET)) 
		return 0;
	if ( write_float(fp, grid, grid_size) < grid_size)
		return 0;

    return 1;
} /* write_grid_file */

int load_binary_grid(TSI_FILE *fp, float *grid, unsigned int grid_size) {
    printf_dbg2("\tload_binary_grid(): %d\n", grid_size);
    if (fseek(fp, (long) 0, SEEK_SET)) return 0;
	if ( read_float(fp, grid, grid_size) < grid_size)
		return 0;

    return 1;
} /* load_grid_file */



unsigned int read_float(TSI_FILE *fp, float *grid, unsigned int nelems)
{
    unsigned int err;

    err = fread(grid, sizeof(float), nelems, fp);
    if (err < nelems) {
	fprintf(stderr,"read_binary: fread returned %u of %u\n", err, nelems);
    }
    return err;
} /* read_float */



unsigned int write_float(TSI_FILE *fp, float *grid, unsigned int nelems)
{
	unsigned int err;

	err = fwrite(grid, sizeof(float), nelems, fp);
	if (err < nelems) {
		fprintf(stderr,"write_binary: fwrite returned %u of %u\n",err,nelems);
	}

	return err;
} /* write_float */

/* end of tsi_io.c */

