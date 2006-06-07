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
        printf_dbg("tsi_backup_simulation(): dumping BestAI...\n");
        if (s) {
            /* dump BCM_c and CC_c */
        } else {
            /* dump CC_c twice (first BCM_c) */
        }
    }
    return 1;
} /* tsi_backup_simulation */



int tsi_restore_simulation(tsi *t, int i, int s)
{
    if (!t->resume) return 0;
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

/* end of file tsi_resume.c */
