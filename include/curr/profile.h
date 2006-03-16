/*
 * Data structures used to profile the algorithm 
 *
 *
 *
 */
 
#ifndef _PROFILE_H
#define _PROFILE_H
 
 
struct {
   unsigned long acorni, newAcorni,
                      backtr, powint, locate, gcum,
                      cova3,
                      ctable,
                      dssdll,
                      gauinv,
                      getindx,
                      krige,
                      ksol,
                      picksup,
                      readdata,
                      readparm,
                      sdsim,
                      setrot,
                      setsupr,
                      sortem,
                      sqdist,
                      srchnd,
                      srchsupr;
    int ident;
} profile;
 
 
void newProfile (void);

void showResults (void);

void profBegin(char *);

void profEnd(char *);

void profPrint(char *);
 
#endif /* _PROFILE_H */

