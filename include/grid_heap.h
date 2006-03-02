#ifndef _GRID_HEAP_H_
#define _GRID_HEAP_H_

#include <stdio.h>

typedef struct grid_type {
    FILE *fp;       /* file where the grid is stored */

    float *grid,    /* alligned grid */
          *pointer; /* original grid pointer */

    int swappable,  /* swappable flag. set to null when grid is in use */
        dirty,      /* dirty grid flag. set to null if grid contents are disposable */
        next_grid;  /* next free grid on stack */
} grid;


typedef struct grid_heap_type {
    grid *g;                 /* grid records array */
    unsigned int grid_size;  /* absolute grid size */
    
    int nodes,               /* number of threads */
        rank;                /* thread rank */

    int heap_size;           /* number of grids in heap */
    int alloc_grids;         /* number of allocated grids in the heap */
    int curr_grids;          /* number of grids in memory */

    int next_grid;           /* next free grid in stack */
    int use_fs;              /* swap to file flag */
} grid_heap;

/* starts a new grid heap */
grid_heap *new_heap(int ngrids, int nodes, int rank, int use_fs,
                    unsigned int xsize, unsigned int ysize, unsigned int zsize);

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

#endif
