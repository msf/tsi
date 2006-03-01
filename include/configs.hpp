// Filename.h
#ifndef CONFIGS_H_
#define CONFIGS_H_
//

#include "const.h"

class Configs
{
	char *group[CONFIGS_MAX_PARAMS];
	char *name[CONFIGS_MAX_PARAMS];
	char *value[CONFIGS_MAX_PARAMS];
	int elements;

	//char *file_content;

public:
	Configs(void);
	virtual ~Configs(void);

	int ReadConfigFile(char *file);

	char* GetParam(char* group, char* param_name);
	float GetParamDouble(char* group, char* param_name);
	int GetParamInt(char* group, char* param_name);
};


#endif // CONFIGS_H_
