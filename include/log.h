#ifndef __LOG_H
#define __LOG_H

#include "tsi.h"

typedef struct log_type {
    char buf[192];
    TSI_FILE *fp;
} log_t;

log_t *new_log(char *filename);

char *logbuf(log_t *l);

char *reset_buffer(log_t *l);

int commit_log(log_t *l);

void delete_log(log_t *l);
 
#endif
