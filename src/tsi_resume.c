/* tsi_resume.c */

#include <string.h>
#include "tsi.h"
#include "tsi_resume.h"
#include "tsi_io.h"


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

    if (t->dump_cm) {
        printf_dbg("\ttsi_backup_simulation(): dumping CC... type: %d\n", t->dump_file);
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
	char buf[2048];
	char *fname;

    if (!t->resume) return 0;
	/* resume allways starts at sim = 0, iteration = 1 */
	if( s != 0 && s!= 0) return 0;

	/* load currBCM, currBAI
	 * at the begginning of the iteration, nextBAI & nextBCM will be 
	 * adopted has the currBCM and currBAI, so what we load here, is the
	 * next*
	 */

	t->nextBCM_idx = new_grid(t->heap);
	t->nextBCM = load_grid(t->nextBCM_idx);


	if ((k = get_key(reg, "RESUME", "BCM")) == NULL) {
		printf_dbg("new_tsi(%d): failed to allocate seismic grid\n", t->proc_id);
		delete_tsi(t);
		return NULL;
	}
	sprintf(buf, "%s%s", t->seismic_path, get_string(k));
	if ((fp = fopen(buf, "r")) == NULL) {
		printf("ERROR: Failed to open the BCM grid file: %s\n", buf);
		delete_tsi(t);
		return NULL;
	}
	if (!tsi_read_grid(t, fp, t->nextBCM, t->seismic_file)) {
		printf("ERROR: Failed to load BCM file: %s\n", buf);
		delete_tsi(t);
		return NULL;
	}
	dirty_grid(t->heap, t->nextBCM_idx);

	t->nextBAI_idx = new_grid(t->heap);
	t->nextBAI = load_grid(t->nextBAI_idx);

	if ((k = get_key(reg, "RESUME", "BAI")) == NULL) {
		printf_dbg("new_tsi(%d): failed to allocate seismic grid\n", t->proc_id);
		delete_tsi(t);
		return NULL;
	}
	sprintf(buf, "%s%s", t->seismic_path, get_string(k));
	if ((fp = fopen(buf, "r")) == NULL) {
		printf("ERROR: Failed to open the BAI grid file: %s\n", buf);
		delete_tsi(t);
		return NULL;
	}
	if (!tsi_read_grid(t, fp, t->nextBAI, t->seismic_file)) {
		printf("ERROR: Failed to load BAI file: %s\n", buf);
		delete_tsi(t);
		return NULL;
	}
	dirty_grid(t->heap, t->nextBAI_idx);

    return 1;
} /* tsi_restore_iteration */




/* end of file tsi_resume.c */
