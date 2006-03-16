
#ifndef _LOG_H
#define _LOG_H 


class Log
{
public:
	FILE* logFile;

	/* Methods */
	Log(char* folderName);
	virtual ~Log();
	void WriteIterationNumber(int iterNum);
	void WriteSimulationNumber(int simulNum);
	void WriteMessage(char* message);
	void WriteActionTime(char* action, float time);
	void WriteDump(char* dump);
	void WriteResult(char* name, double value);
	void WriteSeparator();
	void Close();
};

#endif

