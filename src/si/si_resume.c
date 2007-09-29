/* si_resume.c */

#include "tsi_resume.h"
#include "si_resume.h"


/* local prototypes */
int si_write_grid(si *eng, TSI_FILE *fp, float *grid, int type, char *desc);


int dump_synthetic_grid(si *s, float *g, int it, int sim)
{
    char filename[1024], desc[128];
    TSI_FILE *fp;

    printf_dbg("\tdump_synthetic_grid(): dumping SY...\n");
    sprintf(desc, "SY_%d_%d", it, (sim * s->n_procs + s->proc_id));        
    sprintf(filename, "%s%s.tsi", s->dump_path, desc);
    fp = create_file(filename);
    if (!si_write_grid(s, fp, g, s->dump_file, desc)) {
        printf_dbg("\tdump_synthetic_grid(): failed to dump SY\n");
        return 0;
    }
    close_file(fp);
    return 1;
} /* dump_synthetic_grid */



int dump_reflections_grid(si *s, float *g, int it, int sim)
{
    char filename[1024], desc[128];
    TSI_FILE *fp;

    printf_dbg("\tdump_reflections_grid(): dumping RG...\n");
    sprintf(desc, "SY_CR_%d_%d", it, (sim * s->n_procs + s->proc_id));        
    sprintf(filename, "%s%s.tsi", s->dump_path, desc);
    fp = create_file(filename);
    if (!si_write_grid(s, fp, g, s->dump_file, desc)) {
        printf_dbg("\tdump_synthetic_grid(): failed to dump RG\n");
        return 0;
    }
    close_file(fp);
    return 1;
} /* dump_reflections_grid */



int si_write_grid(si *eng, TSI_FILE *fp, float *grid, int type, char *desc) {
	switch(type) {
//		case CARTESIAN_FILE:
//			return write_cartesian_grid(fp, grid, eng->grid_size);
		case TSI_ASCII_FILE:
			return write_gslib_grid(fp, grid, eng->xsize, eng->ysize, eng->zsize, desc);
		case TSI_BIN_FILE:
			return write_tsi_grid(fp, grid, eng->xsize, eng->ysize, eng->zsize);
			/*
			   case SGEMS_FILE:
			   return write_sgems_grid(fp, grid, t->grid_size));
			 */
		default:
			fprintf(stderr, "ERROR: Unknown grid file type!\n");
			return 0;
	} /* switch */ 
} /* si_write_grid */

/* end of file si_resume.c */
