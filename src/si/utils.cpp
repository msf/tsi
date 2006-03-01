#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "const.h"


/*void WaitForKeyPress()
{
	fflush(stdin);
	while(1)
	{
		if(kbhit())
		return;
	}
}*/

void pdebug1(char *s)
{
	if (DEBUG_LEVEL >= 1)
		fprintf(stdout, "%s\n", s);
}

void pdebug2(char *s)
{
	if (DEBUG_LEVEL >= 2)
		fprintf(stdout, "%s\n", s);
}

char *Trim (char *s) {
  int n;
  while (*s && *s <= 32) ++ s;
  for (n = (int)strlen(s) - 1; n >= 0 && s[n] <= 32; -- n) s[n] = 0;
  return s;
}

/******************************************
 *
 *  nthroot(x,n) - Cal the nth root of x, x**(1/n).
 *
 *
 *******************************************/
double nthroot(double xx, int nn)
{
	int vi, i, doitagain;
	double vdelta,diff,testval,midpnt,lower,upper;
	double xu,epsilon,sign;

	epsilon = 0.00000001;
	xu = xx;
	sign = 1.0;
	if( xu < 0.0){
		sign = -1.0;
		xu = -xu;
	}
	lower = 0.0;
	if(xu >= 1.0) {
		upper = 1.5*sqrt(xu); // use tagged sqrt().
	} else {
		upper = 1.0;
	}
	doitagain = 1;

	i = 0;                           /* binary search loop */
	while( i < 1000 && doitagain == 1){
		midpnt = lower + ((upper - lower)/2.0);
		testval = 1.0;
		for(vi = 0; vi < nn; vi ++)
			testval = midpnt*testval;
		diff = testval - xu;
		vdelta = diff;
		if( vdelta < 0.0)
			vdelta =  0.0 - vdelta;
		if( vdelta > epsilon)
			doitagain = 1;
		if( vdelta <= epsilon)
			doitagain = 0;
		if( diff < 0.0 )
			lower = midpnt;
		if( diff > 0.0 )
			upper = midpnt;
		if( diff == 0.0 )
			doitagain = 0;
		i = i + 1;
	}
	return(sign*midpnt);
}

int fileLineCount(char* Filename)
{
	FILE* fp;
	char line[GLOBAL_LINE_SIZE];
	int n = 0;

	if ((fp=fopen(Filename, "r")) == NULL) 
	{
		printf("Canno open file: %s\n", Filename);
		return ERROR;
	}
	else
	{
		while(fgets(line, GLOBAL_LINE_SIZE, fp) != NULL)
			n++;
	}
	fclose(fp);
	return n;
}


