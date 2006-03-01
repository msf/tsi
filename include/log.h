
#ifndef _LOG_H
#define _LOG_H 


typedef struct type_Log
{
	FILE* logFile;
} tLog;


	/* Methods */
extern	void Log(char* folderName, tLog** this);
extern	void _Log(tLog* this);
extern	void Log_WriteIterationNumber(int iterNum, tLog* this);
extern	void Log_WriteSimulationNumber(int simulNum, tLog* this);
extern	void Log_WriteMessage(char* message, tLog* this);
extern	void Log_WriteActionTime(char* action, float time, tLog* this);
extern	void Log_WriteDump(char* dump, tLog* this);
extern	void Log_WriteResult(char* name, double value, tLog* this);
extern	void Log_WriteSeparator(tLog* this);
extern	void Log_Close(tLog* this);


#endif

