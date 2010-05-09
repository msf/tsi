
typedef struct {
    char *map;
    int bit_count;
    int map_size;
} mask_t;


/* creates a new (all bits unset) mask */
mask_t * mask_new(int bit_count);

void mask_free(mask_t *mask);
int mask_set(mask_t *mask, int pos);
int mask_unset(mask_t *mask, int pos);
int mask_isset(mask_t *mask, int pos);

void mask_unsetall(mask_t *mask);
void mask_setall(mask_t *mask);
