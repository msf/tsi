#include <math.h>
#include "dss.h"

static 
float my_roundf(const float val)
{
    return ((val - floorf(val)) < 0.5) ? floorf(val) : ceilf(val);
}


int getPos(int x, int y, int z, int xlen, int xylen)
{	
	return (x + (xlen * (y -1)) + (xylen * (z-1)));
}


void get3Dcoords(int ind, int xlen, int xylen, int *x, int *y, int *z)
{
	int iy, iz;
	int t;
	
	/* the same way that we subtract y & z in getPos, we increment here! */
	ind++;

	iz = ind / xylen;
	*z = iz +1;

	t = (ind % xylen);
	iy = t / xlen;
	*y = iy +1;

	*x = ind % xlen;
}


/*     Gets the coordinate index location of a point within a grid 
 min     origin at the center of the first cell 
 siz     size of the cells 
 loc     location of the point being considered 

 coordinate index should be checked to see if its inside given grid by the caller.
 */

int getIndex(float min, float siz, float loc)
{
	return (int) my_roundf( 1 + (loc - min) / siz );
}

/* getting absolute coords from relative coords */
int getAbsolutePos(float base, float siz, int index)
{
	return base + (float) (index -1) * siz;
}

/* !Calculo do equivalente valor gaussiano        (SDSIM) */
float compute_gaussian_equiv(float cmean, unsigned size, harddata_point_t *point)
{
    unsigned low, i, j;

    float vmy;
    if (cmean <= point[0].val) {
        return	point[0].gauss_cprob;
    }
    if (cmean >= point[size-1].val) {
        return point[size-1].gauss_cprob;
    }

    low = 0;
    i = size/2;
    // binary search for value closer to cmean in global histogram 
    do {
        if(point[i].val < cmean) {
            j = (i-low)/2;
            low = i;
            i += j;
        } else if(point[i].val > cmean) {
            i -= (i-low)/2;

        } else
            break;
    } while( low + 2 < i);

    vmy = point[i-1].gauss_cprob + 
		(cmean - point[i-1].val) * 
        (point[i].gauss_cprob - point[i-1].gauss_cprob) / 
		(point[i].val - point[i-1].val);

    return vmy;
}

int cmpfloat(const void *a, const void *b)
{
	float r = *(float *)a - *(float *)b ; 

	return (int) my_roundf(r);
}

int cmpharddata_point_val(const void *a, const void *b)
{
	harddata_point_t *p1 = (harddata_point_t *) a;
	harddata_point_t *p2 = (harddata_point_t *) b;

	return (int) my_roundf( p1->val - p2->val );
}

int cmpharddata_point_gauss_cprob(const void *a, const void *b)
{
	harddata_point_t *p1 = (harddata_point_t *) a;
	harddata_point_t *p2 = (harddata_point_t *) b;
	float r = p1->gauss_cprob - p2->gauss_cprob;
	return (int) my_roundf(r);
}

int cmpvalue_index(const void *a, const void *b)
{
	value_index_t *v1 = (value_index_t *) a;
	value_index_t *v2 = (value_index_t *) b;

	return v1->index - v2->index;
}

