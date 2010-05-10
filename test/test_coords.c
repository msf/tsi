#include <assert.h>

#include "dssutils.h"
#include "si_utils.h"



#define NX 3
#define NY 4
#define NZ 5


int test_origin_coords(int nxy, int nx)
{
    int i, j;

    // dss origin coord is 1x1x1, tsi origin coord is 0x0x0
    // So, we should allways have a +1 offset in each coordinate
    // in dss, grid starts at 1, up to MAX_LEN
    // in tsi, grid starts at 0, up to MAX_LEN-1
    i = get_grid_point( 0, 0, 0, nx, nxy);
    assert(i == 0);

    j = getPos( 1, 1, 1, nx, nxy);
    assert(j == 1);

    return 0;
}

int test_max_coords(int nx, int ny, int nz)
{
    int nxy = nx * ny;
    int max = nxy * nz;
    int i, j;

    j = get_grid_point(nx-1, ny-1, nz-1, nx, nxy);
    assert(j == max -1);

    j = getPos(nx, ny, nz, nx, nxy);
    assert(j ==  max);

    return 0;
}


int main(int arvc, char *argv[])
{
    int nxy = NX * NY;
    int x, y, z, i, j, t =0;


    int max = NX * NY * NZ;
    int min = 0;

    test_origin_coords(nxy, NX);
    test_max_coords(NX, NY, NZ);
    test_ilegal_coords(NX, NY, NZ);

    for(x = 0; x < NX; x++)
        for(y = 0; y < NY; y++)
            for(z = 0; z < NZ; z++) {
                j = get_grid_point(x, y, z, NX, nxy);
                assert(j < max);
                assert(j >= min);
                // diff origin point
                t = getPos(x+1, y+1, z+1, NX, nxy);
                assert(t <= max);
                assert(t > min);

                // diff is +1 array offset
                assert(j+1 == t);
            }


    return 0;
}

