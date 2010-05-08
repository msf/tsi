#include <assert.h>

#include "math_random.h"

typedef unsigned long ulong_t;

ulong_t get_ulong(long seed)
{
    genrand_t *g = genrand_new(seed);
    ulong_t l = genrand_ulong(g);
    genrand_free(g);
    return l;
}

double get_real(long seed)
{
    genrand_t *g = genrand_new(seed);
    double l = genrand_real(g);
    genrand_free(g);
    return l;
}


void test_mt_long(long seed)
{
    ulong_t l1, l2;

    l1 = get_ulong(seed);
    l2 = get_ulong(seed);

    assert(l1 == l2);
}

void test_mt_real(long seed)
{
    double l1, l2;

    l1 = get_real(seed);
    l2 = get_real(seed);

    assert(l1 == l2);
}

void test_static_long(long seed)
{
    ulong_t l1, l2;

    tsi_seed_random(seed);

    l1 = get_ulong( seed );
    l2 = tsi_random();

    assert(l1 == l2);
}

void test_static_real(long seed)
{
    double l1, l2;

    tsi_seed_random(seed);

    l1 = get_real(seed);
    l2 = tsi_random_real();

    assert(l1 == l2);
}

void test_many_long(long seed, int count)
{
    ulong_t l1, l2;
    int i;
    genrand_t * g = genrand_new(seed);

    tsi_seed_random(seed);

    for(i = 0; i < count; i++) {
        l1 = genrand_ulong(g);
        l2 = tsi_random();
        assert(l1 == l2);
    }
    genrand_free(g);
}

void test_many_real(long seed, int count)
{
    double l1, l2;
    int i;
    genrand_t * g = genrand_new(seed);

    tsi_seed_random(seed);

    for(i = 0; i < count; i++) {
        l1 = genrand_real(g);
        l2 = tsi_random_real();
        assert(l1 == l2);
    }
    genrand_free(g);
}


void test_noneq_long(long seed)
{
    ulong_t l1, l2;

    l1 = get_ulong(seed);
    l2 = get_ulong(seed-3);

    assert(l1 != l2);
}

void test_noneq_real(long seed)
{
    double l1, l2;

    l1 = get_real(seed);
    l2 = get_real(seed-3);

    assert(l1 != l2);
}


int main(int argc, char *argv[])
{

    long a[] = { 0, -1, 9, 32766, 5794, 3883, 30914, 3159 };
    int i;

    for(i =0; i < sizeof(a); i++) {
        test_static_long( a[i] );
        test_static_real( a[i] );
        test_mt_long( a[i] );
        test_mt_real( a[i] );
        test_many_long( a[i], 100 );
        test_many_long( a[i], 100 );
        test_noneq_long( a[i] );
        test_noneq_real( a[i] );
    }
    return 0;
}


