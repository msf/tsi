#ifndef __LOG_H
#define __LOG_H

#include "registry.h"

#define ERROR(_log_,_func_x_,_str_x_) {\
        fprintf( (_log_)->logFile, "ERROR: %s failed in %s, file %s, line %d\n",\
         (_func_x_), (_str_x_), __FILE__, __LINE__); \
	}


typedef struct log_type {
    FILE *logFile;
    char *logBuf;
    int verbose;
    int procID;
    int iterNum;
    int simulNum;
} log_t;

log_t *new_log(registry *, int proc_id);

void delete_log(log_t *l);

void log_msg(log_t *l, char *fmt, ...);

void log_iteration_number(log_t *l, int iterNum);
void log_simulation_number(log_t *l, int simulNum);
void log_message(log_t *l, int level, char* message);
void log_action_time(log_t *l, int level, char* action, float time);
void log_string(log_t *l, char* dump);
void log_result(log_t *l, int level, char* name, double value);
void log_separator(log_t *l);
void log_close(log_t *l);

#endif
