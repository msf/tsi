#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "log.h"

#define MAXNAME	1024
Log::Log(char* folderName)
{
	char logNameBuffer[MAXNAME];

	strcpy(logNameBuffer, folderName);
	strcat(logNameBuffer, "tsi_log");
	strcat(logNameBuffer, ".log");

	this->logFile = fopen(logNameBuffer, "w");	
	if(logFile == NULL){
		fprintf(stderr,"TSI.Log could not open %s for writting\n Exiting..\n",logNameBuffer);
		exit(1);
	}
}

Log::~Log(void)
{
}

void Log::WriteIterationNumber(int iterNum)
{
	fprintf(this->logFile, "[%d] - Iteration: %2d\n", iterNum, iterNum);
	fflush(this->logFile);
}

void Log::WriteSimulationNumber(int simulNum)
{
	fprintf(this->logFile, "\t[%d] - Simulation: %2d\n", simulNum, simulNum);
	fflush(this->logFile);
}

void Log::WriteMessage(char* message)
{
	fprintf(this->logFile, "+ - %s\n", message);
	fflush(this->logFile);
}

void Log::WriteActionTime(char* action, float time)
{
	fprintf(this->logFile, "\t+ - %s ..... [ %.3fs ]\n", action, time);
	fflush(this->logFile);
}

void Log::WriteDump(char* dump)
{
	fwrite(dump, sizeof(char), strlen(dump), this->logFile);
	fflush(this->logFile);
}

void Log::WriteResult(char* name, double value)
{
	fprintf(this->logFile, "+ - %s = %.5f\n", name, value);
	fflush(this->logFile);
}

void Log::WriteSeparator()
{
	fprintf(this->logFile, "=======================================================================\n");
	fflush(this->logFile);
}

void Log::Close()
{
	fflush(this->logFile);
	fclose(this->logFile);
}


