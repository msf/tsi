#ifndef LOG_H
#define LOG_H 

#include <stdio.h>
#include "tsi_io.h"


typedef struct log_type {
    TSI_FILE *log_file;
} log_t;


log_t *new_log(filename);

void close_log(log_t *log);

int log(log_t *l, char *msg);

#endif

