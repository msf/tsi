#ifndef _REGISTRY_H_
#define _REGISTRY_H_

#include <stdio.h>

/* struct that defines a 64 bit resolution registry key */
typedef struct reg_key_type {
    char *name;
    char *value;
    int type;   /* 0-string; 1-integer; 2-floating point */
    union {
        long long int lli;
        double dbl;
    } kval;
    struct reg_key_type *next;
} reg_key;

/* struct for registry lists */
typedef struct registry_type {
    char *name;
    reg_key *klist;
    struct registry_type *next;
} registry;

/* builds a new registry */
registry *new_registry(char *filename);

/* deletes a registry from memory */
void delete_registry(registry *reg);

/* merges a file to given registry */
int merge_registry(registry **reg, char *filename);

/* searches for an given entry and returns the key */
reg_key *get_key(registry *reg, char *section, char *parm);

/* parses the key as string and returns a pointer */
char *get_string(reg_key *k);

/* parses the key as an integer and returns the value */
int get_int(reg_key *k);

/* parses the key as a long integer and returns the value */
long get_long(reg_key *k);

/* parses the key as a long long integer and returns the value */
long long int get_llint(reg_key *k);

/* parses the key as a float and returns the value */
float get_float(reg_key *k);

/* parses the key as a double and returns the value */
double get_double(reg_key *k);

#endif
