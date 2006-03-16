/* timer.h 
   simple cpu timer api.
   */


double getElapsedTime(struct timeval *start, struct timeval *stop);
int getCurrTime(struct timeval *tv);
