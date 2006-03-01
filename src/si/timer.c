#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "timer.h"

void Timer(tTimer **this)
{
	*this = (tTimer *) malloc(sizeof(tTimer));
}

void _Timer(tTimer *this)
{
	free(this);
}


double Timer_getElapsedTimeC(struct timeval *start, struct timeval *stop)
{
	double d = stop->tv_sec - start->tv_sec;
	double c = stop->tv_usec - start->tv_usec;
	if(c < 0) {
		d--;
		c =  start->tv_usec - stop->tv_usec;
	}
	c /= 1000000; // number of usecs in a second
	return  d + c;
}

void Timer_Start(tTimer* this)
{
	gettimeofday(&this->tv1,NULL);
}

void Timer_Stop(tTimer* this)
{
	gettimeofday(&this->tv2,NULL);
}

void Timer_Reset(tTimer* this)
{
	this->tv1.tv_sec = 0;
	this->tv1.tv_usec = 0;
	this->tv2.tv_sec = 0;
	this->tv2.tv_usec = 0;
}

double Timer_GetElapsedTime(tTimer* this)
{
	return Timer_getElapsedTimeC(&this->tv1,&this->tv2);
}

double Timer_GetElapsedTimeInSeconds(tTimer* this)
{
	return Timer_getElapsedTimeC(&this->tv1,&this->tv2);
}

