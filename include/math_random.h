/*
 * header file for the use of the internal random number generator.
 * it uses the mersenne_twister.
 */


void tsi_seed_random(long seed);

unsigned long tsi_random(void);

double	tsi_random_real(void);


