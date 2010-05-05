#include <stdio.h>
#include <string.h>

#include "mask.h"

int mask_test1(mask_t *m, int len)
{
    int i;
    for(i = 0; i < len; i++) {
        if(mask_isset(m, i)) {
            printf("ERROR 1 in: %d\n", i);
            return 1;
        }
    }
    return 0;
}

int mask_test2(mask_t *m, int len)
{
    int i;
    mask_setall(m);
    for(i = 0; i < len; i++) {
        if(!mask_isset(m, i)) {
            printf("ERROR 2 in: %d\n", i);
            return 1;
        }
    }
    return 0;
}

int mask_test3(mask_t *m, int len)
{
    int i;
    mask_unsetall(m);
    for(i = 0; i < len; i++) {
        if(mask_isset(m, i)) {
            printf("ERROR 3 in: %d\n", i);
            return 1;
        }
    }
    return 0;
}

int mask_test4(mask_t *m, int len)
{
    int i = len / 2;
    mask_unsetall(m);
    mask_set(m, i);
    if( mask_isset(m, i-1) || mask_isset(m, i+1) || !mask_isset(m, i) ) {
        printf("ERROR 4 in: %d\n", i);
        return 1;
    }
    return 0;
}

int mask_test(int len)
{
    int i;
    mask_t *m;

    m = mask_new(len);

    if( mask_test1(m, len) ||
        mask_test2(m, len) ||
        mask_test3(m, len) ||
        mask_test4(m, len) ) {
        mask_free(m);
        return 1;
    }
    return 0;
}


int main(int argc, char *argv[])
{
    int len, i;

    len = atoi(argv[1]);
    for(i = 1 << 30; i > len; i >>= 1)
        printf("%d -> %d\n",i, mask_test(i));

}
