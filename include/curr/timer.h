#ifndef TIMER_H_
#define TIMER_H_

#include <sys/time.h>
#include <time.h>


typedef struct type_Timer
{
	struct timeval tv1, tv2;
} tTimer;


extern void Timer(tTimer **this);
extern void _Timer(tTimer *this);

extern void Timer_Start(tTimer* this);

extern void Timer_Stop(tTimer* this);

extern void Timer_Reset(tTimer* this);

extern double Timer_GetElapsedTime(tTimer* this);

extern double Timer_GetElapsedTimeInSeconds(tTimer* this);


#endif 



