#ifndef SECTIONSDEFS_H_
#define SECTIONSDEFS_H_

typedef struct type_SectionsDefs
{
	int* sections;
	//int sections[Z_VALUES];
	int max_sections;

	int w_c;

	int acum;
} tSectionsDefs;

extern void SectionsDefs(int sections_max, int w_c, tSectionsDefs* this);
extern void _SectionsDefs(tSectionsDefs* this);

extern int SectionDefs_String2Sections(char* sections, tSectionsDefs* this);
extern int SectionDefs_GetAcum(tSectionsDefs* this);

#endif // SECTIONSDEFS_H_

