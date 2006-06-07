/* tsi_resume.c */

#include <tsi_resume.h>


int tsi_backup_simulation(tsi *t, int i, int s)
{
    if (t->dump_ai) {
        printf_dbg("tsi_backup_simulation(): dumping AI...\n");
        /* dump AI */
    }
    if (t->dump_cm) {
        printf_dbg("tsi_backup_simulation(): dumping CC...\n");
        /* expand and dump CC */
    }
    if (t->dump_bcm) {
        printf_dbg("tsi_backup_simulation(): dumping BCM...\n");
        if (s) {
            /* expand and dump nextBCM */
        } else {
            /* expand and dump CC as first BCM */
        }
    }
    if (t->dump_bai) {
        printf_dbg("tsi_backup_simulation(): dumping BAI...\n");
        if (s) {
            /* dump nextBAI */
        } else {
            /* dump AI as first BAI */
        }
    }
    if (t->resume) {
        printf_dbg("tsi_backup_simulation(): dumping best results...\n");
        if (s) {
            /* dump nextBCM_c and nextBAI in binary */
        } else {
            /* dump CC_c as first BCM_c */
            /* dump AI as first BAI in binary */
        }
        if (t->global_best.value > t->last_corr.value) {
            printf_dbg("tsi_backup_simulation(): dumping BestAI...\n");
            t->last_corr.value = t->global_best.value; 
            t->last_corr.proc_id = t->global_best.proc_id; 
            /* dump bestAI (sim numbered for sequence eval on resume) */
            /* include correlation on 3rd line for fast resume (what about cart-grid format?) */
        }
    }
    return 1;
} /* tsi_backup_simulation */



int tsi_restore_simulation(tsi *t, int i, int s)
{
    if (!t->resume) return 0;
    /* recurse each node sim sequence for best results -> fuck up... develop new model for sim sequence... */
    return 1;
} /* tsi_restore_simulation */



int tsi_backup_iteration(tsi *t, int i)
{
    if (t->dump_bcm) {
        printf_dbg("tsi_backup_iteration(): dumping BCM...\n");
        /* expand and dump nextBCM */
    }
    if (t->dump_bai) {
        printf_dbg("tsi_backup_iteration(): dumping BAI...\n");
        /* dump nextBAI */
    }
    return 1;
} /* tsi_backup_iteration */



int tsi_restore_iteration(tsi *t, int i)
{
    if (!t->resume) return 0;
    return 1;
} /* tsi_restore_iteration */



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
}



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
}

/* end of file tsi_resume.c */
