/*
 * Data structures used to profile the algorithm 
 *
 *
 *
 */
 
#include <stdio.h>
#include "profile.h" 
 

void newProfile ()
{
    profile.acorni =
    profile.newAcorni =
    profile.backtr =
    profile.cova3 =
    profile.ctable =
    profile.dssdll =
    profile.gauinv =
    profile.gcum =
    profile.getindx =
    profile.krige =
    profile.ksol =
    profile.locate =
    profile.picksup =
    profile.powint =
    profile.readdata =
    profile.readparm =
    profile.sdsim =
    profile.setrot =
    profile.setsupr =
    profile.sortem =
    profile.sqdist =
    profile.srchnd =
    profile.srchsupr = 0;

    profile.ident = 0;
} /* newProfile */
 
 
void showResults()
{
    printf("\nProfiling results:\n");
    printf("acorni ____ %lu\n", profile.acorni);
    printf("newAcorni _ %lu\n", profile.newAcorni);
    printf("backtr ____ %lu\n", profile.backtr);
    printf("cova3 _____ %lu\n", profile.cova3);
    printf("ctable ____ %lu\n", profile.ctable);
    printf("dssdll ____ %lu\n", profile.dssdll);
    printf("gauinv ____ %lu\n", profile.gauinv);
    printf("gcum ______ %lu\n", profile.gcum);
    printf("getindx ___ %lu\n", profile.getindx);
    printf("krige _____ %lu\n", profile.krige);
    printf("ksol ______ %lu\n", profile.ksol);
    printf("locate ____ %lu\n", profile.locate);
    printf("picksup ___ %lu\n", profile.picksup);
    printf("powint ____ %lu\n", profile.powint);
    printf("readdata __ %lu\n", profile.readdata);
    printf("readparm __ %lu\n", profile.readparm);
    printf("sdsim _____ %lu\n", profile.sdsim);
    printf("setrot ____ %lu\n", profile.setrot);
    printf("setsupr ___ %lu\n", profile.setsupr);
    printf("sortem ____ %lu\n", profile.sortem);
    printf("sqdist ____ %lu\n", profile.sqdist);
    printf("srchnd ____ %lu\n", profile.srchnd);
    printf("srchsupr __ %lu\n", profile.srchsupr);
    printf("\n\n");

} /* showResults */



void profBegin (char *x)
{
   int i;
   
   for (i = 0; i < profile.ident; i++) printf("    ");
   printf("--> %s\n", x);
   profile.ident++;
} /* profBegin */



void profEnd (char *x)
{
   int i;
   
   profile.ident--;
   for (i = 0; i < profile.ident; i++) printf("    ");
   printf("<-- %s\n", x);
} /* profEnd */


void profPrint (char *x)
{
   int i;
   
   for (i = 0; i < profile.ident; i++) printf("    ");
   printf("--- %s\n", x);
} /* profPrint */

