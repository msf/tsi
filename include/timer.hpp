#ifndef TIMER_H_
#define TIMER_H_

#include <time.h>


class Timer
{
private:
	clock_t tv1, tv2;
	double time;

public:

	Timer(void)
	{
	}

	~Timer(void)
	{
	}

	void Start()
	{
		tv1 = clock();
	}

	void Stop()
	{
		tv2 = clock();
	}

	void Reset()
	{
		tv1 = 0;
		tv2 = 0;
	}

	double GetElapsedTime()
	{
		return ((tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0));
	}

	double GetElapsedTimeInSeconds()
	{
		return ((tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0)) / 1000;
	}

	double GetInstantTime()
	{
		return ((clock() - tv1)/(CLOCKS_PER_SEC / (double) 1000.0));
	}

	double GetInstantTimeInSeconds()
	{
		return ((clock() - tv1)/(CLOCKS_PER_SEC / (double) 1000.0)) / 1000;
	}
};

#endif 



