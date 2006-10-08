#ifndef _TSI_RESUME_H
#define _TSI_RESUME_H

#include <tsi.h>


int tsi_backup_simulation(tsi *t, int i, int s);

int tsi_restore_simulation(tsi *t, int i, int s);

int tsi_backup_iteration(tsi *t, int i);

int tsi_restore_iteration(tsi *t, int i);

#endif
