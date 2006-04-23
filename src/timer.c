#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "debug.h"
#include "timer.h"

double getElapsedTime(struct timeval *start, struct timeval *stop)
{
	double d = stop->tv_sec - start->tv_sec;
	double c = stop->tv_usec - start->tv_usec;
	if(c < 0) {
		printf_dbg2("microsecs has wrapped around\n");
		d--;
		c =  start->tv_usec - stop->tv_usec;
	}
	c /= 1000000;
	return  d + c;
}
	
int getCurrTime(struct timeval *tv)
{
    return gettimeofday(tv,NULL);
}

