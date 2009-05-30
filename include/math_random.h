/*
 * header file for the use of the internal random number generator.
 * it uses the mersenne_twister.
 */

#define N 624

/* state structure */
typedef struct {
    unsigned long state[N]; /* the array for the state vector  */
    int left;
    int initf;
    unsigned long *next;
} mtrand_t;

/* give abstract name to the random generator */
typedef mtrand_t genrand_t;

genrand_t *genrand_new(long seed);
void genrand_free(genrand_t *state);


unsigned long genrand_ulong(genrand_t *state);

double genrand_real(genrand_t *state);


void tsi_seed_random(long seed);

unsigned long tsi_random(void);

double	tsi_random_real(void);


