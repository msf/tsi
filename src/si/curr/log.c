#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "log.h"

#define MAXNAME	1024
void Log(char* folderName, tLog** this)
{
	char logNameBuffer[MAXNAME];

	*this = (tLog *) malloc(sizeof(tLog));
	
	strcpy(logNameBuffer, folderName);
	strcat(logNameBuffer, "tsi_log");
	strcat(logNameBuffer, ".log");

	(*this)->logFile = fopen(logNameBuffer, "w");	
	if((*this)->logFile == NULL){
		fprintf(stderr,"TSI.Log could not open %s for writting\n Exiting..\n",logNameBuffer);
		exit(1);
	}
}

void _Log(tLog* this)
{
}

void Log_WriteIterationNumber(int iterNum, tLog* this)
{
	fprintf(this->logFile, "[%d] - Iteration: %2d\n", iterNum, iterNum);
	fflush(this->logFile);
}

void Log_WriteSimulationNumber(int simulNum, tLog* this)
{
	fprintf(this->logFile, "\t[%d] - Simulation: %2d\n", simulNum, simulNum);
	fflush(this->logFile);
}

void Log_WriteMessage(char* message, tLog* this)
{
	fprintf(this->logFile, "+ - %s\n", message);
	fflush(this->logFile);
}

void Log_WriteActionTime(char* action, float time, tLog* this)
{
	fprintf(this->logFile, "\t+ - %s ..... [ %.3fs ]\n", action, time);
	fflush(this->logFile);
}

void Log_WriteDump(char* dump, tLog* this)
{
	fwrite(dump, sizeof(char), strlen(dump), this->logFile);
	fflush(this->logFile);
}

void Log_WriteResult(char* name, double value, tLog* this)
{
	fprintf(this->logFile, "+ - %s = %.5f\n", name, value);
	fflush(this->logFile);
}

void Log_WriteSeparator(tLog* this)
{
	fprintf(this->logFile, "=======================================================================\n");
	fflush(this->logFile);
}

void Log_Close(tLog* this)
{
	fflush(this->logFile);
	fclose(this->logFile);
}


