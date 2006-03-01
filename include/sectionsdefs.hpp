#ifndef SECTIONSDEFS_H_
#define SECTIONSDEFS_H_

class SectionsDefs
{
public:
	int* sections;
	//int sections[Z_VALUES];
	int max_sections;

	int w_c;

	int acum;

	SectionsDefs(int sections_max, int w_c);
	~SectionsDefs(void);

	int String2Sections(char* sections);
	int GetAcum();
};

#endif // SECTIONSDEFS_H_

