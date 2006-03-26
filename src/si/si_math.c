/* si_math.c */

#include <stdlib.h>
#include "debug.h"
#include "registry.h"
#include "grid_heap.h"
#include "si.h"


/* wavelet functions */
float point_value(si *s, int point);
float index_value(si *s, int index);



int make_reflections_grid(si *s, float *AI, float *RG) {
	printf_dbg("make_reflections_grid(): called but not implemented\n");
    return 1;
} /* make_reflections_grid */



int make_synthetic_grid(si *s, float *RG, float *SY) {
	printf_dbg("make_syntethic_grid(): called but not implemented\n");
    return 1;
} /* make_synthetic_grid */



int make_correlations_grid(si *s, float *seismic, float *SY) {
	printf_dbg("make_correlations_grid(): called but not implemented\n");
    return 1;
} /* make_correlations_grid */



int expand_correlations_grid(cm_grid *cmg, float *CM) {
	printf_dbg("expand_correlations_grid(): called but not implemented\n");
    return 1;
} /* expand_correlations_grid */



float point_value(si* s, int point) {
    printf_dbg2("Point %d is at array position %d.\n", point, (s->wavelet_used_values / 2) + point);
    return s->values[(s->wavelet_used_values / 2) + point];
} /* point_value */



float index_value(si *s, int index) {
    return s->values[index];
} /* index_value */

/* end of file si_math.c */
