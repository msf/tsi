/* timer.h 
   simple cpu timer api.
   */
#ifdef WIN32
#include <sys/types.h> 
#include <sys/timeb.h> 
#include <winsock2.h> //timeval
#else
#include <sys/time.h>
#include <time.h>
#endif

#define my_time			timeval
double getTimeDiff(struct my_time *start, struct my_time *stop);
int getCurrTime(struct my_time *tv);


#ifdef WIN32
/* define gettimeofday() */
static inline int gettimeofday(struct timeval *tv, void *ignore) 
{
	struct _timeb buf;
	_ftime_s(&buf);
	tv->tv_sec = buf.time;
	tv->tv_usec = buf.millitm * 1000;
	return 0;
}
#endif

