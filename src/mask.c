#include <stdlib.h>

#include "mask.h"

mask_t * mask_new(int bit_count)
{
    int map_size;
    mask_t *msk = malloc(sizeof(mask_t));
    if(!msk)
        return NULL;

    msk->bit_count = bit_count;

    map_size = (bit_count / 8) + 1;
    msk->map = (char *) calloc(map_size, sizeof(char));
    if(!msk->map) {
        free(msk);
        return NULL;
    }
    msk->map_size = map_size;

    return msk;
}

void mask_free(mask_t *mask)
{
    if( !mask  || !mask->map)
        return ;

    free(mask->map);
    free(mask);
}

int mask_set(mask_t *mask, int pos)
{
    if( !mask  || !mask->map)
        return -1;

    int n = pos / 8;
    int off = pos % 8;
    mask->map[n] |= 1 << off;
    return 0;
}

int mask_unset(mask_t *mask, int pos)
{
    if( !mask  || !mask->map)
        return -1;

    int n = pos / 8;
    int off = pos % 8;
    mask->map[n] &= ~(1 << off);
    return 0;
}

int mask_isset(mask_t *mask, int pos)
{
    if( !mask  || !mask->map)
        return -1;

    int n = pos / 8;
    int off = pos % 8;
    return ( mask->map[n] >> off ) & 1 ;
}






