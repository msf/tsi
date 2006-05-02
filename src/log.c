/* log.c */

#include <stdio.h>
#include <string.h>

#include "debug.h"
#include "registry.h"
#include "log.h"

#define MAXNAME	1024

log_t *new_log(registry *reg, int proc_id) {
	char logName[MAXNAME];
    log_t *new_log;
	reg_key *k;
    
    new_log = (log_t *) tsi_malloc(sizeof(log_t));
    if (new_log == NULL) {
        printf_dbg("new_log: failed to allocate space for new log\n");
        return NULL;
    }

	new_log->procID = proc_id;
	
			/* use log or output path */
	if ( ((k = get_key(reg, "GLOBAL", "LOG_PATH")) == NULL) ) 
		if( ((k = get_key(reg, "GLOBAL", "OUTPUT_PATH")) == NULL) ) {
		fprintf(stderr, "MISSING LOG_PATH AND OUTPUT PATH on config file, aborting\n");	
		delete_log(new_log);
		return NULL;
	}

	sprintf(logName,"%s/tsi-proc-%d.log", get_string(k),proc_id);

	new_log->verbose = 0;
	k = get_key(reg, "GLOBAL", "VERBOSE");
	if(k != NULL) {
		new_log->verbose = get_int(k);
	}

	new_log->logFile = fopen(logName, "w");	
	if(new_log->logFile == NULL){
		if(new_log->verbose){
			fprintf(stderr,"new_log(): ERROR, could not open %s for writting\n using only stdout..\n",logName);
		}
		else
			fprintf(stderr,"new_log(): ERROR, could not open %s for writting\n enabling verbosive mode..\n",logName);
		new_log->logFile = stdout;
		new_log->verbose = 0;
	}

	new_log->simulNum = 0;
	new_log->iterNum = 0;
    return new_log;
}


void delete_log(log_t *l) {
    if (l) {
        /* close file */
		log_close(l);
        tsi_free(l);
    } else {
        printf_dbg("delete_log: received NULL as parameter!\n");
    }
}
 
void log_indent(char *buf, int nivel) 
{
	int i;

	for(i = 0; i < nivel && i < 9; i++) {
		buf[i] = '\t';
	}
	buf[i] = '\0';
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

void log_message(log_t *l, int nivel, char* message)
{
	char buf[10];
			
	log_indent(buf, nivel);
	fprintf(l->logFile, "(%d,%d,%d) %s+ - %s\n", l->procID, l->iterNum, l->simulNum, buf, message);
	if(l->verbose) {
		printf("(%d,%d,%d) %s+ - %s\n", l->procID, l->iterNum, l->simulNum, buf, message);
		fflush(stdout);
	}
	fflush(l->logFile);
}

void log_action_time(log_t *l, int nivel, char* action, float time)
{
	char buf[10];
			
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
	char buf[10];
			
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
	char buf[80];
	int i;
	for(i = 0; i < 80; i++)
		buf[i] = '=';
	buf[i] = '\0';

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
