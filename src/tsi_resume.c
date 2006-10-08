/* tsi_resume.c */

#include <string.h>
#include "tsi_resume.h"


int tsi_backup_simulation(tsi *t, int i, int s)
{
    float *g;
    cm_grid *cmg;
    char filename[1024], desc[128];
    int g_idx;
    TSI_FILE *fp;

    if (t->dump_ai) {
        printf_dbg("\ttsi_backup_simulation(): dumping AI...\n");
        /* dump AI */
        g = load_grid(t->heap, t->ai_idx);
        sprintf(desc, "AI_%d_%d", i, (s * t->n_procs + t->proc_id));        
        sprintf(filename, "%s%s.tsi", t->dump_path, desc);
        fp = create_file(filename);
        if (!tsi_write_grid(t, fp, g, t->dump_file, desc)) {
            printf_dbg("\ttsi_backup_simulation(): failed to dump AI\n");
            return 0;
        }
        close_file(fp);
        clear_grid(t->heap, t->ai_idx);
    }

    if (t->dump_bai) {
        printf_dbg("\ttsi_backup_simulation(): dumping BAI...\n");
        sprintf(desc, "BAI_%d_%d", i, (s * t->n_procs + t->proc_id));        
        sprintf(filename, "%s%s.tsi", t->dump_path, desc);
        if (s) {
            /* dump nextBAI */
            g_idx = t->nextBAI_idx;
        } else {
            /* dump AI as first BAI */
            g_idx = t->ai_idx;
        }
        g = load_grid(t->heap, g_idx);
        fp = create_file(filename);
        if (!tsi_write_grid(t, fp, g, t->dump_file, desc)) {
            printf_dbg("\ttsi_backup_simulation(): failed to dump BAI\n");
            return 0;
        }
        close_file(fp);
        clear_grid(t->heap, g_idx);
    }

    if (t->dump_cm) {
        printf_dbg("\ttsi_backup_simulation(): dumping CC...\n");
        /* expand and dump CC */
        cmg = get_cmgrid(t->si_eng);
        load_cmgrid(cmg);
        g_idx = new_grid(t->heap);
        g = load_grid(t->heap, g_idx);
        expand_correlations_grid(cmg, g);
        sprintf(desc, "CC_%d_%d", i, (s * t->n_procs + t->proc_id));        
        sprintf(filename, "%s%s.tsi", t->dump_path, desc);
        fp = create_file(filename);
        if (!tsi_write_grid(t, fp, g, t->dump_file, desc)) {
            printf_dbg("\ttsi_backup_simulation(): failed to dump CC\n");
            return 0;
        }
        close_file(fp);
        delete_grid(t->heap, g_idx);
        clear_cmgrid(cmg);
    }

    if (t->dump_bcm) {
        printf_dbg("\ttsi_backup_simulation(): dumping BCM...\n");
        sprintf(desc, "BCM_%d_%d", i, (s * t->n_procs + t->proc_id));        
        sprintf(filename, "%s%s.tsi", t->dump_path, desc);
        if (s) {
            /* expand and dump nextBCM */
            cmg = t->nextBCM_c;
        } else {
            /* expand and dump CC as first BCM */
            cmg = get_cmgrid(t->si_eng);
        }
        load_cmgrid(cmg);
        g_idx = new_grid(t->heap);
        g = load_grid(t->heap, g_idx);
        expand_correlations_grid(cmg, g);
        fp = create_file(filename);
        if (!tsi_write_grid(t, fp, g, t->dump_file, desc)) {
            printf_dbg("\ttsi_backup_simulation(): failed to dump CC\n");
            return 0;
        }
        close_file(fp);
        delete_grid(t->heap, g_idx);
        clear_cmgrid(cmg);
    }

    if (t->resume) {
        printf_dbg("\ttsi_backup_simulation(): dumping best results...\n");
        if (s) {
            /* dump nextBCM_c and nextBAI in binary */
        } else {
            /* dump CC_c as first BCM_c */
            /* dump AI as first BAI in binary */
        }
        if (t->global_best.value > t->last_corr.value) {
            printf_dbg("\ttsi_backup_simulation(): dumping BestAI...\n");
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
    float *g;
    cm_grid *cmg;
    char filename[1024], desc[128];
    int g_idx;
    TSI_FILE *fp;

    if (t->dump_bai) {
        printf_dbg("\ttsi_backup_iteration(): dumping BAI...\n");
        sprintf(desc, "BAI_%d", i);        
        sprintf(filename, "%s%s.tsi", t->dump_path, desc);
        g_idx = t->nextBAI_idx;
        g = load_grid(t->heap, g_idx);
        fp = create_file(filename);
        if (!tsi_write_grid(t, fp, g, t->dump_file, desc)) {
            printf_dbg("\ttsi_backup_iteration(): failed to dump BAI\n");
            return 0;
        }
        close_file(fp);
        clear_grid(t->heap, g_idx);
    }

    if (t->dump_bcm) {
        printf_dbg("\ttsi_backup_simulation(): dumping BCM...\n");
        sprintf(desc, "BCM_%d", i);        
        sprintf(filename, "%s%s.tsi", t->dump_path, desc);
        g_idx = t->nextBCM_idx;
        g = load_grid(t->heap, g_idx);
        fp = create_file(filename);
        if (!tsi_write_grid(t, fp, g, t->dump_file, desc)) {
            printf_dbg("\ttsi_backup_simulation(): failed to dump CC\n");
            return 0;
        }
        close_file(fp);
        clear_grid(t->heap, g_idx);
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