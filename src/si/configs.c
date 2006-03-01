#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "configs.h"


void Configs(tConfigs** this)
{
	*this = (tConfigs*) malloc(sizeof(tConfigs));
}

void _Configs(tConfigs* this)
{
	free(this);
}


int Configs_ReadConfigFile(char *file, tConfigs * this)
{
    FILE *fp;
    char file_line[GLOBAL_LINE_SIZE];
    

    /* open file for reading */
    fp = fopen(file, "r");
    if (fp == NULL) {
        perror("Error reading config file");
        return -1;
    }

    char *token;
    unsigned int line = 0;

    this->elements = 0;

    /* parse string */
    while(fgets(file_line, GLOBAL_LINE_SIZE, fp) != NULL) {
        line++;

        if(file_line[0] != '[') {
            continue;
        }
        

        /* get group name */
        token = (char *) strtok(&file_line[1],"]");
        this->group[this->elements] = strdup(token);

        /* get parameter name */
        token = (char *) strtok(NULL,"=");
        if (token == NULL){
            fprintf(stderr,"Config file with bad syntax line(%d), aborting.\n",line);
            return -1;
        }
        this->name[this->elements] = strdup(token);
        
        //get value name
        token = (char *) strtok(NULL,"\n\r");
        if (token == NULL){
            fprintf(stderr,"Config file with bad syntax line(%d), aborting.\n",line);
            return -1;
        }
        this->value[this->elements] = strdup(token);

		//printf("Config file: %s, %s, %s\n",this->group[this->elements],this->name[this->elements],this->value[this->elements]);
        
        this->elements++;
    }


    // close file descriptor
    fclose(fp);

    return 0;
}



char* Configs_GetParam(char* group, char* param_name, tConfigs* this)
{
	int i;
	for(i=0; i<this->elements; i++){
		
		if ((strcmp(this->group[i], group)==0) && (strcmp(this->name[i], param_name)==0)) 
			return strdup(this->value[i]);
	}
	return NULL;
}


float Configs_GetParamDouble(char* group, char* param_name, tConfigs* this)
{
	int i;
	for(i=0; i<this->elements; i++)
		if ((strcmp(this->group[i], group)==0) && (strcmp(this->name[i], param_name)==0))
			return atof(this->value[i]);

	return ERROR;
}

int Configs_GetParamInt(char* group, char* param_name, tConfigs* this)
{
	int i;
	for( i=0; i<this->elements; i++)
		if ((strcmp(this->group[i], group)==0) && (strcmp(this->name[i], param_name)==0))
			return atoi(this->value[i]);

	return ERROR;
}


