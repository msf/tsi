#include <math.h>
#include "dss.h"

static
int my_roundf(const float val)
{
    return (int) ((val - floorf(val)) < 0.5) ? floorf(val) : ceilf(val);
}

int getPos(int x, int y, int z, int xlen, int xylen)
{
	return (x + (xlen * (y -1)) + (xylen * (z-1)));
}


void get3Dcoords(int ind, int xlen, int xylen, int *x, int *y, int *z)
{
	int iy, iz;
	int t;

	ind--;

	iz = ind / xylen;
	*z = iz +1;

	t = (ind % xylen);
	iy = t / xlen;
	*y = iy +1;

	*x = (ind % xlen) + 1;
}


/*     Gets the coordinate index location of a point within a grid
 min     origin at the center of the first cell
 siz     size of the cells
 loc     location of the point being considered

 coordinate index should be checked to see if its inside given grid by the caller.
 */

int getIndex(float min, float siz, float loc)
{
	return my_roundf( 1 + (loc - min) / siz );
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

    /* div by 0 bug fix */
    while(point[i-1].val == point[i].val)
        i++;

    vmy = point[i-1].gauss_cprob + (cmean - point[i-1].val) *
        (point[i].gauss_cprob - point[i-1].gauss_cprob) /
		(point[i].val - point[i-1].val);

    return vmy;
}

int cmpfloat(const void *a, const void *b)
{
	float r = *(float *)a - *(float *)b ;
	if( r < 0)
		return -1;
	else if(r > 0)
		return 1;
	else
		return 0;
}

int cmpharddata_point_val(const void *a, const void *b)
{
	harddata_point_t *p1 = (harddata_point_t *) a;
	harddata_point_t *p2 = (harddata_point_t *) b;
	float x = p1->val - p2->val;
	if( x < 0)
		return -1;
	else if(x > 0)
		return 1;
	else
		return 0;
}

int cmpharddata_point_gauss_cprob(const void *a, const void *b)
{
	harddata_point_t *p1 = (harddata_point_t *) a;
	harddata_point_t *p2 = (harddata_point_t *) b;
	float r = p1->gauss_cprob - p2->gauss_cprob;
	if( r < 0)
		return -1;
	else if(r > 0)
		return 1;
	else
		return 0;
}

int cmpvalue_index(const void *a, const void *b)
{
	value_index_t *v1 = (value_index_t *) a;
	value_index_t *v2 = (value_index_t *) b;

	return v1->index - v2->index;
}

int cmpvalue_index_by_value(const void *a, const void *b)
{
	value_index_t *v1 = (value_index_t *) a;
	value_index_t *v2 = (value_index_t *) b;

	float r =  v1->value - v2->value;
	if( r < 0)
		return -1;
	else if(r > 0)
		return 1;
	else
		return 0;
}

int cmpvalue_point(const void *a, const void *b)
{
	value_point_t *v1 = (value_point_t *) a;
	value_point_t *v2 = (value_point_t *) b;

	float r =  v1->value - v2->value;
	if( r < 0)
		return -1;
	else if(r > 0)
		return 1;
	else
		return 0;
}
