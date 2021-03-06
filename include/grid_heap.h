#ifndef _GRID_HEAP_H_
#define _GRID_HEAP_H_

#include "tsi_io.h"

#ifdef WIN32
#define strdup _strdup
#endif

typedef struct grid_type {
    float *grid,       /* alligned grid */
          *pointer;    /* original grid pointer */

    int swappable,     /* swappable flag. set to null when grid is in use */
        dirty,         /* dirty grid flag. set to null if grid contents are disposable */
        next_grid,     /* next free grid on stack */
        valid,
        first_load;    /* flag to zero the grid */
    unsigned int size; /* new size if different from default */

    char *filename;    /* filename to use when creating file */
    TSI_FILE *fp;      /* file where the grid is stored */

} grid;


typedef struct grid_heap_type {
    grid *g;                 /* grid records array */
    unsigned int grid_size;  /* absolute grid size */

    int nodes,               /* number of threads */
        rank;                /* thread rank */

    int heap_size,           /* number of grids in heap */
        alloc_grids,         /* number of allocated grids in the heap */
        threshold,           /* swap threshold */
        use_fs;              /* swap to file flag */

    int curr_grids,          /* number of grids in memory */
        reads,               /* number of reads in files */
        writes,              /* number of writes to files */
        max_grids;           /* max number of grids allocated */

    int next_grid;           /* next free grid in stack */
} grid_heap;

/* starts a new grid heap */
grid_heap *new_heap(int nodes, int rank, int heap_size, int swap_thr, int use_fs, char *path, unsigned int grid_size);

/* allocates a new grid */
int new_grid(grid_heap *h);

/* return a pointer to the grid */
float *load_grid(grid_heap *h, int g);

/* grid not in use and has not been modified */
int clear_grid(grid_heap *h, int g);

/* grid not in use and has changes */
int dirty_grid(grid_heap *h, int g);

/* delete the grid */
void delete_grid(grid_heap *h, int g);

/* delete the heap */
void delete_heap(grid_heap *h);

/* print heap data */
void print_heap_data(grid_heap *h);

/* set a grid size different from default */
void set_grid_size(grid_heap *h, int idx, unsigned int size);

#endif /* _GRID_HEAP_H */
