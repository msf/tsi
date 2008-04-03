#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dss.h"
#include "dss_legacy.h"
#include "debug.h"
#include "log.h"
#include "tsi_io.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)


#define IXL		1
#define IYL		2
#define IZL		3
#define IVRL	4
#define IWT		0
#define ISECVR	0
#define ITRANS	1
#define ISMOOTH 0
#define ISVR    1
#define ISWT    2



int load_harddata_file(log_t *l, char *filename, harddata_t *harddata) 
{

    int i, j, m;
	char line[256];
	harddata_point_t point;
	TSI_FILE *fp;

    if ((fp = fopen(filename, "r")) == NULL) {
		ERROR(l, "fopen()", filename);
        return 1;
    }
	
    /* ignore gslib header */
	if( !read_gslib_header(l, fp, 4) ) {
		ERROR(l, "load_harddata_file()", "read_gslib_header()");
		return 1;
	}

	/* first find out how many values */
	j = 0;
	while(fgets(line, 255, fp) != NULL)
		j++;
	printf_dbg("read_harddata_file(): %d lines\n",i);

	harddata->point = (harddata_point_t *) tsi_malloc(sizeof(harddata_point_t) * j);
	if(harddata->point == NULL) {
		ERROR(l, "load_harddata_file", "NOT ENOUGH MEMORY");
		exit(1);
	}

	fseek(fp, 0, SEEK_SET); // back to start of file

    /* ignore gslib header */
	if( !read_gslib_header(l, fp, 4) ) {
		ERROR(l, "load_harddata_file()", "read_gslib_header()");
		return 1;
	}
	
	i = 0;
	m = 0;
	while( (m = fscanf(fp,"%f %f %f %f", &point.x, &point.y, &point.z, &point.val)) != EOF ) {
		if(m != 4)
			log_print(l,"load_harddata_file(): ERROR -  %d can't parse data values, line: %d\n",m, 6+(i/4));

		/* skip points that aren't in the [min, max] interval */
		if(point.val < harddata->min_value || point.val > harddata->max_value)
			continue;

		harddata->point[i++] = point;
	}

	/* trim point array to correct size */
	if( j > i) {
		unsigned int siz = i * sizeof(harddata_point_t);
		harddata_point_t *new = (harddata_point_t *) tsi_malloc(siz);
		if(new == NULL) {
			ERROR(l, "load_harddata_file", "NOT ENOUGH MEMORY");
			exit(1);
		}
		memcpy(new, harddata->point, siz);
		tsi_free(harddata->point);
		harddata->point = new;
	}

	harddata->point_count = i;

    fclose(fp);
	return 0;
} /* load_harddata_file */


int readdata(log_t *l, 
			harddata_t *harddata,
			general_vars_t * general)
{
	int i, j;
	double w;
	float vrg;
	int ierr;
	value_index_t *temp;

    printf_dbg2("readdata(): begin\n");
	/* Function Body */
	if (harddata->point_count <= 1) {
		ERROR(l, "readdata", "no harddata/histogram");
		return -1; /* ERROR */
	}
	/* !Sort data by value: */
	qsort(harddata->point, harddata->point_count, sizeof(harddata_point_t), cmpharddata_point_val);

	temp = (value_index_t *) tsi_malloc(sizeof(value_index_t) * harddata->point_count);
	if(temp == NULL) {
		ERROR(l, "readdata", "NOT ENOUGH MEMORY");
		exit(1);
	}

	/* Compute:
	 * 	average of value.
	 *  the cumulative probabilities (gaussian distribution)
	 */
	double cp = 0;
	double oldcp = 0;
	harddata->average = 0;
	int x, y, z;
	int ig = -1;
	harddata_point_t point;
	for (j = 0; j < harddata->point_count; ++j) {
		point = harddata->point[j];
		harddata->average += point.val;

		/* see if value is inside grid */
		x = getIndex(general->xmn, general->xsiz, point.x); 
		y = getIndex(general->ymn, general->ysiz, point.y);
		z = getIndex(general->zmn, general->zsiz, point.z);
		if((x >= 0 && x < general->nx) &&
		   (y >= 0 && y < general->ny) &&
		   (z >= 0 && z < general->nz) ) {

			ig++;
			temp[ig].index = getPos(x,y,z,general->nx,general->nxy);
			temp[ig].value = point.val;
		}

		cp = (double) (j+1) / (double) harddata->point_count;
		w = (cp + oldcp) / 2;
		oldcp = cp;

		ierr = gauinv( w, &vrg);
		if (ierr == 1) 
			harddata->point[j].gauss_cprob = general->nosim_value;
		else
			harddata->point[j].gauss_cprob = vrg;
		if(j < 100 || harddata->point_count - j < 100)
			printf_dbg2("readdata(): vrtr[%u] = %f,\tvrgtr[%u] = %f\n",
				j, 
				harddata->point[j].val,
				j,
				harddata->point[j].gauss_cprob);
	}
	harddata->average /= harddata->point_count;

	/* sort by grid index, so that candidates to the same cell become adjacent */
	qsort(temp, ig, sizeof(value_index_t), cmpvalue_index);

	harddata->in_grid = (value_index_t *) tsi_malloc(sizeof(value_index_t) * ig);
	if(harddata->in_grid == NULL) {
		ERROR(l, "readdata", "NOT ENOUGH MEMORY");
		exit(1);
	}

	/* temp might have several values candidate for same grid cell.
	 * (when harddata has more resolution that our grid cells this is very common)
	 * make shure in_grid has only one value (the 1st seen ) for each different index
	 */
	 
	j = 0;
	i = 0;
	int last= -1;
	do {
		
		if(last == temp[i++].index)
			continue;
	
		last = temp[i-1].index;
		harddata->in_grid[j].index = last;
		harddata->in_grid[j].value = temp[i].value;
		if(j < 100)
			printf_dbg2("in_grid[%u] idx = %u, value = %f\n", 
				j,
				harddata->in_grid[j].index,
				harddata->in_grid[j].value);
		j++;
	} while( i < ig) ;
	harddata->in_grid_count = j;

	tsi_free(temp);
	/* trim in_grid to only the needed size */
	if(ig > harddata->in_grid_count) {
		unsigned int siz = harddata->in_grid_count * sizeof(value_index_t);
		value_index_t *new = (value_index_t *) tsi_malloc(siz);
		if(new == NULL) {
			ERROR(l, "readdata", "NOT ENOUGH MEMORY");
			exit(1);
		}
		memcpy(new, harddata->in_grid, siz);
		tsi_free(harddata->in_grid);
		harddata->in_grid = new;
	}

	printf_dbg("harddata():\tharddata->in_grid_count: %d\tharddata->point_count: %d\n",
			harddata->in_grid_count, harddata->point_count);


	double r1;
	harddata->variance = 0;
	for (i = 0; i < harddata->point_count; ++i) {
		/* Computing 2nd power */
		r1 = harddata->point[i].val - harddata->average;
		harddata->variance += r1 * r1;
	}
	harddata->variance /= harddata->point_count;

	printf_dbg("readdata(): average: %f, variance: %f\n",
			harddata->average, harddata->variance);

	return 0;
} /* readdata */


