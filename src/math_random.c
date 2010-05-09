/*
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's real version.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <stdlib.h>

#include "math_random.h"

/* Period parameters */
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

static mtrand_t mt_state;

/* initializes state[N] with a seed */
static void mtrand_init(mtrand_t *mt, unsigned long s)
{
    unsigned long *state = mt->state;
    int j;
    state[0]= s & 0xffffffffUL;
    for (j=1; j<N; j++) {
        state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array state[].                     */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        state[j] &= 0xffffffffUL;  /* for >32 bit machines */
    }
    mt->left = 1;
    mt->initf = 1;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
static void init_by_array(mtrand_t *mt, unsigned long init_key[], int key_length)
{
    unsigned long *state = mt->state;
    int i, j, k;

    mtrand_init(mt, 19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) {
		   	state[0] = state[N-1];
		   	i=1;
	   	}
        if (j>=key_length)
			j=0;
    }
    for (k=N-1; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) {
			state[0] = state[N-1];
		   	i=1;
	   	}
    }

    state[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
    mt->left = 1;
    mt->initf = 1;
}

static void mtrand_next_state(mtrand_t *mt)
{
    unsigned long *p=mt->state;
    int j;

    /* if init_genrand() has not been called, */
    /* a default initial seed is used         */
    if (mt->initf==0)
		mtrand_init(mt, 5489UL);

    mt->left = N;
    mt->next = mt->state;

    for (j=N-M+1; --j; p++)
        *p = p[M] ^ TWIST(p[0], p[1]);

    for (j=M; --j; p++)
        *p = p[M-N] ^ TWIST(p[0], p[1]);

    *p = p[M-N] ^ TWIST(p[0], mt->state[0]);
}

/* generates a random number on [0,0xffffffff]-interval */
static unsigned long mtrand_int32(mtrand_t *mt)
{
    unsigned long y;

    if (--mt->left == 0)
		mtrand_next_state( mt );
    y = (*mt->next)++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
static long mtrand_int31(mtrand_t *mt)
{
    return (long)(mtrand_int32(mt)>>1);
}

/* generates a random number on [0,1)-real-interval */
static double mtrand_real2(mtrand_t *mt)
{
    return ((double)mtrand_int32(mt)) * (1.0 / 4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
static double mtrand_real3(mtrand_t *mt)
{
    return ((double)mtrand_int32(mt) + 0.5) * (1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
static double mtrand_res53(mtrand_t *mt)
{
    unsigned long a=mtrand_int32(mt)>>5, b=mtrand_int32(mt)>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */


/*  CODE ADDED BY MIG.  #####################################  */
#ifdef WIN32
#undef OLD_RAND
#endif

#ifdef OLD_RAND
#include <stdlib.h>
#endif

void tsi_seed_random(long seed)
{
#ifdef OLD_RAND
    srandom(seed);
#else
    mt_state.initf = 0;
    mt_state.left = 1;
	mtrand_init( &mt_state, seed);
#endif
}

unsigned long tsi_random(void)
{
#ifdef OLD_RAND
    return (unsigned long)random();
#else
	//TODO: try int28()
	return mtrand_int31(&mt_state);
#endif
}

double tsi_random_real(void)
{
#ifdef OLD_RAND
    return ((double)random()) / RAND_MAX;
#else
	return mtrand_real2(&mt_state);
#endif
}



mtrand_t *genrand_new(long seed)
{
    genrand_t *state = malloc(sizeof(genrand_t));
    if( ! state )
        return NULL;

    state->initf = 0;
    state->left = 1;
    mtrand_init(state, seed);
    return state;
}

void genrand_free(genrand_t *state)
{
    if(state)
        free(state);
    state = NULL;
}


unsigned long genrand_ulong(genrand_t *state)
{
    return mtrand_int31(state);
}

double genrand_real(genrand_t *state)
{
    return mtrand_real2(state);
}

