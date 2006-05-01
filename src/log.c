/* log.c */

#include <stdio.h>
#include <string.h>

#include "debug.h"
#include "registry.h"
#include "log.h"

#define BUFFER_SIZE	1024

log_t *new_log(registry *reg, int proc_id) {
    char *log_path, empty_path;
    log_t *new_log;
    reg_key *k;
    
    if ((new_log = (log_t *) tsi_malloc(sizeof(log_t))) == NULL) {
        printf_dbg("new_log: failed to allocate space for new log\n");
        return NULL;
    }
    
    if ((new_log->logBuf = (char *) tsi_malloc(BUFFER_SIZE * sizeof(char))) == NULL) {
        printf_dbg("new_log: failed to allocate space for new log\n");
        delete_log(new_log);
        return NULL;
    }

    new_log->procID = proc_id;
	
    /* use log or output path */
    empty_path = 0;
    if ((k = get_key(reg, "GLOBAL", "LOG_PATH")) == NULL) {
       if ((k = get_key(reg, "GLOBAL", "OUTPUT_PATH")) == NULL) {
           log_path = &empty_path;
       } else {
           log_path = get_string(k);
       }
    } else {
        log_path = get_string(k);
    }

    sprintf(new_log->logBuf,"%stsi-proc-%d.log", log_path, proc_id);

    if ((k = get_key(reg, "GLOBAL", "VERBOSE")) == NULL)
        new_log->verbose = 0;
    else
        new_log->verbose = get_int(k);

    new_log->logFile = fopen(new_log->logBuf, "w");
    if (new_log->logFile == NULL) {
        if (new_log->verbose) {
            fprintf(stderr,"new_log(): ERROR, could not open %s for writting\n using only stdout..\n", new_log->logBuf);
        } else {
            fprintf(stderr,"new_log(): ERROR, could not open %s for writting\n enabling verbosive mode..\n", new_log->logBuf);
            new_log->logFile = stdout;
            new_log->verbose = 0;
        }   
    }

    new_log->simulNum = 0;
    new_log->iterNum = 0;
    return new_log;
}


void delete_log(log_t *l) {
    if (l) {
        /* close file */
        if (l->logFile) log_close(l);
        if (l->logBuf) tsi_free(l->logBuf);
        tsi_free(l);
    } else {
        printf_dbg("delete_log: received NULL as parameter!\n");
    }
}
 
void log_indent(char *buf, int level)
{
    int i;

    if (level > 8) level = 9; /* limit identation to 9 tabs */
    for (i = 0; i < level; i++) {
        buf[i] = '\t';
    }
    buf[i] = 0;
}

void log_iteration_number(log_t *l, int iterNum)
{
	l->simulNum = 0;
	l->iterNum = iterNum;
	fprintf(l->logFile, "(%d,%d,-) - starting Iteration: %2d\n", l->procID,l->iterNum, iterNum);
	if(l->verbose){
		printf("(%d,%d,-) - starting Iteration: %2d\n", l->procID, l->iterNum, iterNum);
		fflush(stdout);
	}
	fflush(l->logFile);
}

void log_simulation_number(log_t *l, int simulNum)
{
	l->simulNum = simulNum;
	fprintf(l->logFile, "(%d,%d,%d) - starting simulation: %2d\n", l->procID, l->iterNum, l->simulNum, l->simulNum);
	if(l->verbose) {
		printf("(%d,%d,%d) - starting simulation: %2d\n", l->procID, l->iterNum, l->simulNum, l->simulNum);
		fflush(stdout);
	}
	fflush(l->logFile);
}

void log_message(log_t *l, int level, char* message)
{
	char buf[20];
			
	log_indent(buf, level);
	fprintf(l->logFile, "(%d,%d,%d) %s+ - %s\n", l->procID, l->iterNum, l->simulNum, buf, message);
	if(l->verbose) {
		printf("(%d,%d,%d) %s+ - %s\n", l->procID, l->iterNum, l->simulNum, buf, message);
		fflush(stdout);
	}
	fflush(l->logFile);
}

void log_action_time(log_t *l, int nivel, char* action, float time)
{
	char buf[20];
			
	log_indent(buf, nivel);
	fprintf(l->logFile, "(%d,%d,%d) %s+ - %s ..... [ %.3fs ]\n", l->procID, l->iterNum, l->simulNum, buf, action, time);
	if(l->verbose) {
		printf("(%d,%d,%d) %s+ - %s ..... [ %.3fs ]\n", l->procID, l->iterNum, l->simulNum, buf, action, time);
		fflush(stdout);
	}
	fflush(l->logFile);
}

void log_result(log_t *l, int nivel, char* name, double value)
{
	char buf[20];
			
	log_indent(buf, nivel);
	fprintf(l->logFile, "(%d,%d,%d) %s+ - %s = %.5f\n", l->procID, l->iterNum, l->simulNum, buf, name, value);
	if(l->verbose) {
		printf("(%d,%d,%d) %s+ - %s = %.5f\n", l->procID, l->iterNum, l->simulNum, buf, name, value);
		fflush(stdout);
	}
	fflush(l->logFile);
}

void log_string(log_t *l, char* dump)
{
	fwrite(dump, sizeof(char), strlen(dump), l->logFile);
	if(l->verbose) {
		fwrite(dump, sizeof(char), strlen(dump), stdout);
		fflush(stdout);
	}
	fflush(l->logFile);
}

void log_separator(log_t *l)
{
	char buf[100];
	int i;

	for(i = 0; i < 80; i++)
		buf[i] = '=';
	buf[i] = 0;

	fprintf(l->logFile, "%s\n",buf);
	if(l->verbose) {
		printf("%s\n",buf);
		fflush(stdout);
	}
	fflush(l->logFile);
}

void log_close(log_t *l)
{
	fflush(l->logFile);
	fflush(stdout);
	fclose(l->logFile);
}


/* end of file log.c */
