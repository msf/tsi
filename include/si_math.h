/* si_math.h */

#ifndef _SI_MATH_H
#define _SI_MATH_H

#include <stdlib.h>
#include "debug.h"
#include "registry.h"
#include "grid_heap.h"
#include "si.h"

/* build reflections grid */
int make_reflections_grid(si *s, float *AI, float *RG);

/* builds synthetic seismic grid for evaluation */
int make_synthetic_grid(si *s, float *RG, float *SY);

/* builds a correlations grid based on layers configuration */
int make_correlations_grid(si *s, float *seismic, float *SY);

/* expands the correlations grid to a regular grid for DSS */
int expand_correlations_grid(cm_grid *cmg, float *CM);

#endif /* _SI_MATH_H */

/* end of file si_math.h */
