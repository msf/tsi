#include "timer.h"

int getTimeMillis(struct my_time start, struct my_time stop)
{
	return (int) (stop.tv_usec - start.tv_usec);
}

double getTimeDiff(struct my_time start, struct my_time stop)
{
	double d = stop.tv_sec - start.tv_sec;
	double c = stop.tv_usec - start.tv_usec;
	if(c < 0) {
		d--;
		c = 1000000 - start.tv_usec + stop.tv_usec;
	}
	c /= 1000000;
	return  d + c;
}

double getElapsedTime(struct my_time *start, struct my_time *stop)
{
	return getTimeDiff(*start, *stop);
}
	
struct my_time getTime(void)
{
	struct my_time tv;
	gettimeofday(&tv, '\0');
	return tv;
}

int getCurrTime(struct my_time *tv)
{
    return gettimeofday(tv,'\0');
}


