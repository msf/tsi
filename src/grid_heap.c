#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "grid_heap.h"
#include "debug.h"
#include "memdebug.h"
#include "tsi_io.h"

#define ALLIGNMENT 16          /* byte allignment */

grid_heap *new_heap(int nodes, int rank, int heap_size, int swap_thr, int use_fs, int xsize, int ysize, int zsize) {
    grid_heap *h;
    int i;
    char filename[32];
    
    /* initializes the heap record */
    h = (grid_heap *) my_malloc(sizeof(grid_heap));
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
    h->threshold = swap_thr;
    h->use_fs = use_fs;

    /* free grids stack */
    h->next_grid = 0;      /* free grids stack */
    
    /* grid counters*/
    h->alloc_grids = 0;    /* allocated grids counter */
    h->curr_grids = 0;     /* grids in use counter */
    h->reads = 0;
    h->writes = 0;
    h->max_grids = 0;
    
    /* initializes the grid records array */
    h->g = (grid *) calloc((size_t)heap_size, sizeof(grid));
    if (!h->g) {
        printf_dbg2("new_heap(): failed to allocate grid data array!\n");
        delete_heap(h);
        return NULL;
    } /* if */

    /* create the swap files for each grid */
    for(i = 0; i < heap_size; i++) {
        if (h->use_fs) {      
            sprintf(filename, "/tmp/tsi_grid.node%d.%d", rank, i);
            printf_dbg2("new_heap(): filename string >%s<\n", filename);
            h->g[i].fp = open_grid_file(filename);
            if (!h->g[i].fp) {
                printf_dbg2("new_heap(): failed to open grid file!\n");
                delete_heap(h);
                return NULL;
            } /* if */
        }
        h->g[i].next_grid = i+1;  /* set free grid stack */
        h->g[i].swappable = 1;
    } /* for */

    printf_dbg2("new_heap(%d): heap started\n", rank);
    return h;
} /* new_heap */



/* allocates a new grid */
int new_grid(grid_heap *h) {
    grid *g;
    int idx;

    idx = 0;
    
    if (h->alloc_grids < h->heap_size) {       /* check if there are any grids available */

        h->alloc_grids++;        /* increase the number of allocated grids */
        if (h->alloc_grids > h->max_grids) h->max_grids = h->alloc_grids;
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
    
    printf_dbg2("load_grid(): grid req=%d curr_grids=%d alloc_grids=%d\n", idx, h->curr_grids, h->alloc_grids);
    if (idx < h->heap_size) {

        if (h->g[idx].grid) {   /* check if this grid still has its space allocated */
            h->g[idx].swappable = 0;   /* mark grid in use flag */
            return h->g[idx].grid;     /* return pointer to grid space */
        } else if (h->curr_grids < h->threshold) {  /* if #grids is bellow threshold */
            h->curr_grids++;
            printf("Number of grids raised to %d\n", h->curr_grids);
            h->g[idx].grid = (float *) my_malloc(h->grid_size*sizeof(float));
            h->g[idx].pointer = h->g[idx].grid;
            h->g[idx].swappable = 0;
            return h->g[idx].grid;
        }
        
        /* seek for a cleared grid */
        j = h->heap_size;
        for (i = 0; i < h->heap_size; i++) {
            if (h->g[i].swappable && h->g[i].grid)
                j = i;
            if (h->g[i].swappable && !h->g[i].dirty && h->g[i].grid)
                break;
        }
        printf_dbg2("load_grid(): j = %d, i = %d\n", j, i);

        if ((h->use_fs) && (i != j)) { /* found dirty grid */
            /* dump j grid */
            h->writes++;
            printf_dbg2("load_grid(): swapping dirty grid %d to disk\n", j);
            if (!write_grid_file(h->g[j].fp, h->g[j].grid, h->grid_size)) {
                printf_dbg2("load_grid(): failed to dump grid %d\n", j);
                return NULL;
            }
            h->g[j].dirty = 0;
            i = j;
        }

        if (i == h->heap_size) {
            /* no grid available, allocate more ram */
            h->curr_grids++;
            printf("Number of grids raised to %d\n", h->curr_grids);
            h->g[idx].grid = (float *) my_malloc(h->grid_size*sizeof(float));
            h->g[idx].pointer = h->g[idx].grid;
            h->g[idx].swappable = 0;
        } else {
            /* use available grid */
            h->g[idx].grid = h->g[j].grid;
            h->g[idx].pointer = h->g[j].pointer;
            h->g[idx].swappable = 0;
            h->g[j].grid = h->g[j].pointer = NULL;
        }
        
        if ((h->use_fs) && (!h->g[idx].dirty)) {
            /* load grid */
            h->reads++;
            printf_dbg2("load_grid(): loading grid %d to memory\n", idx);
            if (!read_grid_file(h->g[idx].fp, h->g[idx].grid, h->grid_size)) {
                 printf_dbg2("load_grid(): failed to load grid %d\n", idx);
                 return NULL;
            }
        }
        return h->g[idx].grid;
    } else {
        printf_dbg2("load_grid(): grid %d off range\n", idx);
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
    int i;
    if (h) {
        if (h->g) {
            for (i = 0; i < h->heap_size; i++) {
                if (h->g[i].fp) close_file(h->g[i].fp);
                if (h->g[i].grid) my_free(h->g[i].grid);
            }
            my_free(h->g);
        }
        my_free(h);
    }
} /* delete_heap */


/* end of grid_heap.c */
