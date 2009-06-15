#ifndef _SI_UTILS_H_
#define _SI_UTILS_H_

/* give absolute index offset in grid based on x, y, z coordinates */
/* coordinate start at 0,0,0 */
unsigned int get_grid_point(int x, int y, int z, int nx, int ny);


#endif
