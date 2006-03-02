#include <stdlib.h>
#include <stdio.h>
#include "grid_heap.h"
#include "debug.h"

#define ALLIGNMENT 16          /* byte allignment */
#define WHAT ???               /* what define is missing here?!?!?! */


grid_heap *new_heap(int heap_size, int nodes, int rank, int use_fs, int xsize, int ysize, int zsize) {

    grid_heap *h;
    int i, j;
    char filename[32];
    
    /* initializes the heap record */
    h = (grid_heap *) malloc(sizeof(grid_heap));
    if (!h) {
        printf_dbg2("new_heap(): failed to allocate heap!\n");
        return NULL;
    }

    /* machine data */
    h->nodes = nodes;
    h->rank = rank;

    /* heap constants */
    h->heap_size = heap_size;   /* number of grids in heap */
    h->grid_size = (unsigned int)xsize * (unsigned int)ysize * (unsigned int)zsize; // + ALLIGNMENT/sizeof(float);

    /* free grids stack */
    h->next_grid = 0;      /* free grids stack */
    
    /* grid counters*/
    h->alloc_grids = 0;    /* allocated grids counter */
    h->curr_grids = 0;     /* grids in use counter */
    
    /* initializes the grid records array */
    h->g = (grid *) calloc(heap_size, sizeof(grid));
    if (!h->g) {
        printf_dbg2("new_heap(): failed to allocate grid data array!\n");
        delete_heap(h);
        return NULL;
    } /* if */

    /* create the swap files for each grid */
    if (use_fs) {
        for(i = 0; i < heap_size; i++) {
            sprintf(filename, "/tmp/tsi_grid.node%d.%d", rank, i);
            printf_dbg2("new_heap(): filename string >%s<\n", filename);
            h->g[i].fp = fopen(filename, "w+b");
            if (!h->g[i].fp) {
                printf_dbg2("new_heap(): failed to open grid file!\n");
                delete_heap(h);
                return NULL;
            } /* if */
            h->g[i].next_grid = i+1;  /* set free grid stack */
            h->g[i].swappable = 1;
        } /* for */
    } /* if */

    return h;
} /* new_heap */



/* allocates a new grid */
int new_grid(grid_heap *h) {

    grid *g;
    float *f;
    int idx;

    idx = 0;
    
    if (h->alloc_grids < h->heap_size) {       /* check if there are any grids available */

        h->alloc_grids++;        /* increase the number of allocated grids */
        idx = h->next_grid;      /* get index from free grids stack */

        g = h->g + idx;
        h->next_grid = g->next_grid;   /* update free grids stack */
        g->dirty = 1;    /* mark new grid to prevent being cleared and overwritten */

        printf_dbg2("new_grid(): grid %d, next %d\n", idx, h->next_grid);
        printf_dbg2("new_grid(): curr_grids=%d alloc_grids=%d\n", h->curr_grids, h->alloc_grids);
    } else {
        printf_dbg2("new_grid(): no more grids available!\n");
        return -1;
    }

    return idx;
} /* new_grid */



/* return a pointer to the grid */
float *load_grid(grid_heap *h, int idx) {
    int i, j;
    
    if (idx < h->heap_size) {
        printf_dbg2("load_grid(): curr_grids=%d alloc_grids=%d\n", h->curr_grids, h->alloc_grids);

        if (h->g[idx].grid) {    /* check if this grid still has its space allocated */
           h->g[idx].swappable = 0;   /* mark grid in use flag */
           return h->g[idx].grid;     /* return pointer to grid space */
        }

        if (h->curr_grids < h->alloc_grids) {
            /* find a cleared grid */
            j = h->heap_size;
            for (i = 0; i < h->heap_size; i++) {
                if (h->g[i].swappable && h->g[i].grid)
                    j = i;
                if (h->g[i].swappable && !h->g[i].dirty && h->g[i].grid)
                    break;
            }
            if (i != j) { /* found dirty grid */
                /* dump j grid */
                if (fseek(h->g[j].fp, 0, SEEK_SET)) {
                    printf_dbg2("load_grid(): failed to reset file position for dirty grid %d\n", j);
                    return NULL;
                }
                printf_dbg2("load_grid(): swapping grid %d to disk\n", j);
                if (fwrite(h->g[j].grid, sizeof(float), h->grid_size, h->g[j].fp) < h->grid_size) {
                    printf_dbg2("load_grid(): failed to dump grid %d\n", j);
                    return NULL;
                }
                h->g[j].dirty = 0;
            }

            if (i == h->heap_size) {
                /* no grid available, allocate more ram */
                h->curr_grids++;
                h->g[idx].grid = (float *) malloc(sizeof(float));
                h->g[idx].pointer = h->g[idx].grid;
                h->g[idx].swappable = 0;
                if (!h->g[idx].dirty) {
                    /* load grid */
                    if (fseek(h->g[idx].fp, 0, SEEK_SET)) {
                        printf_dbg2("load_grid(): failed to reset file position for grid %d\n", idx);
                        return NULL;
                    }
                    printf_dbg2("load_grid(): loading grid %d to memory\n", idx);
                    if (fread(h->g[idx].grid, sizeof(float), h->grid_size, h->g[idx].fp) < h->grid_size) {
                        printf_dbg2("load_grid(): failed to dump grid %d\n", idx);
                        return NULL;
                    }
                }
            } else {
                h->g[idx].grid = h->g[j].grid;
                h->g[idx].pointer = h->g[j].pointer;
                h->g[idx].swappable = 0;
                h->g[j].grid = h->g[j].pointer = NULL;
                /* load grid */
                if (fseek(h->g[idx].fp, 0, SEEK_SET)) {
                    printf_dbg2("load_grid(): failed to reset file position for grid %d\n", idx);
                    return NULL;
                }
                printf_dbg2("load_grid(): loading grid %d to memory\n", idx);
                if (fread(h->g[idx].grid, sizeof(float), h->grid_size, h->g[idx].fp) < h->grid_size) {
                    printf_dbg2("load_grid(): failed to dump grid %d\n", idx);
                    return NULL;
                }
            }
            return h->g[idx].grid;
        }
    } else {
        printf_dbg2("load_grid(): requested grid %d is off range\n", idx);
    }   
    return NULL;
} /* load_grid */



int clear_grid(grid_heap *h, int idx) {
    if ((idx < h->heap_size) && (h->g[idx].grid != NULL)) {
        h->g[idx].swappable = 1;   /* mark grid as not needed */
        printf_dbg2("clear_grid(): curr_grids=%d alloc_grids=%d\n", h->curr_grids, h->alloc_grids);
        return 1;
    } else {
        printf_dbg2("clear_grid(): requested grid %d is off range or empty\n", idx);
    }
    return 0;
} /* clear_grid */



int dirty_grid(grid_heap *h, int idx) {
    if ((idx < h->heap_size) && (h->g[idx].grid != NULL)) {
        h->g[idx].swappable = 1;   /* mark grid as not needed */
        h->g[idx].dirty = 1;       /* mark grid to be saved before being reused */
        printf_dbg2("dirty_grid(): curr_grids=%d alloc_grids=%d\n", h->curr_grids, h->alloc_grids);
    } else {
        printf_dbg2("dirty_grid(): requested grid %d is off range or empty\n", idx);
    }
    return 0;
} /* dirty_grid */



void delete_grid(grid_heap *h, int idx) {
    if (idx < h->heap_size) {
        h->alloc_grids--;
        h->g[idx].next_grid = h->next_grid;
        h->next_grid = idx;
        h->g[idx].swappable = 1;
        h->g[idx].dirty = 0;
        printf_dbg2("delete_grid(): curr_grids=%d alloc_grids=%d\n", h->curr_grids, h->alloc_grids);
    } else {
        printf_dbg2("delete_grid(): requested grid %d is off range\n", idx);
    }
} /* delete_grid */



void delete_heap(grid_heap *h) {
    /* TODO */
} /* delete_heap */

/* end of grid_heap.c */
