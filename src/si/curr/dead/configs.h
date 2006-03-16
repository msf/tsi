// Filename.h
#ifndef CONFIGS_H_
#define CONFIGS_H_
//

#include "const.h"

typedef struct type_Configs
{
	char *group[CONFIGS_MAX_PARAMS];
	char *name[CONFIGS_MAX_PARAMS];
	char *value[CONFIGS_MAX_PARAMS];
	int elements;
	//char *file_content;
} tConfigs;

extern	void Configs(tConfigs** this);
extern  void _Configs(tConfigs* this);

extern	int Configs_ReadConfigFile(char *file, tConfigs* this);

extern	char* Configs_GetParam(char* group, char* param_name, tConfigs* this);
extern	float Configs_GetParamDouble(char* group, char* param_name, tConfigs* this);
extern	int Configs_GetParamInt(char* group, char* param_name, tConfigs* this);



#endif // CONFIGS_H_
