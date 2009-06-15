
unsigned int get_grid_point(int x, int y, int z, int nx, int ny)
{
    int nxy = nx * ny;
    return (z*nxy)+(nx*y)+x;
}

