#include <stdio.h>
#include "debug.h"
#include "registry.h"
#include "tsi.h"

#ifdef WIN32
void debug_check() { return ; }
#endif

int main(int argc, char *argv[])
{
    /* main TSI objects */
    registry  *r;
    tsi *t;

    int i, res;
    char default_registry_file[] = "tsi.conf";
    char *reg_file;

    if (argc > 1)  /* evaluate if a registry file (or more) was passed as parameter */
        reg_file = argv[1];
    else     /* switch to default file */
    {
        printf_dbg("No registry file passed as argument! Trying default file...\n");
        reg_file = default_registry_file;
    }
    
    r = new_registry(reg_file);   /* attempts to load registry file */

    if (r)   /* evaluate if the first registry file was loaded successfully */
    {
        i = 2;
        while (i < argc)   /* try to load any other parameter files */
        {
            res = 0;
            reg_file = argv[i];
            res = merge_registry(&r, reg_file);
            if (!res) return 1;
        }
    }
    else
        return 1;

    printf_dbg("Registry loaded...\n");

    /* starting new program from registry */
    if (!(t = new_tsi(r))) {
        printf_dbg("Failed to load TSI!\n");
        return 2;
    }

    /* start execution */
    if (!run_tsi(t)) {
        printf_dbg("TSI Failed to execute!\n");
        delete_tsi(t);
        return 3;
    }

    delete_tsi(t);
    return 0;
} /* main */

