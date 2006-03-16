#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sectionsdefs.h"

void SectionsDefs(int sections_max, int w_c, tSectionsDefs* this)
{
	this->max_sections = sections_max;

	this->sections = (int*)malloc((this->max_sections + 1) * sizeof(int));

	this->w_c = w_c;

}

void _SectionsDefs(tSectionsDefs* this)
{
	free(this->sections);
}

int SectionsDefs_String2Sections(char* sections, tSectionsDefs* this)
{
	char* string;
	int ret;

	string = sections;

	int i;
	i = 1;

	sscanf(string, "%d,", &this->sections[i]);
	this->acum = this->sections[i];

	// printf("Sections[%d] = %d\n", i, this->sections[i]);

	for(i++; i<=this->max_sections; i++)
	{
		string = strpbrk(string, ",");

		if (string == NULL)
		{
			printf("ERROR: Problem in SECTIONS definition. Check that the number of sections is %d.\n", this->max_sections);
			return -1;
		}

		string++;

		ret = sscanf(string, "%d,", &this->sections[i]);
		this->acum = this->acum + this->sections[i];

		if (ret <= 0)
		{
			printf("ERROR: Problem in SECTIONS definition. Using [%s].\n", sections);
			return -1;
		}

		// printf("Sections[%d] = %d\n", i, this->sections[i]);
	}

	return 0;
}

int SectionsDefs_GetAcum(tSectionsDefs* this)
{
	return (this->acum);
}


