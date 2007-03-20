#include "timer.h"

double getElapsedTime(struct my_time *start, struct my_time *stop)
{
	double d = stop->tv_sec - start->tv_sec;
	double c = stop->tv_usec - start->tv_usec;
	if(c < 0) {
		d--;
		c = 1000000 - start->tv_usec + stop->tv_usec;
	}
	c /= 1000000;
	return  d + c;
}
	
int getCurrTime(struct my_time *tv)
{
    return gettimeofday(tv,'\0');
}
