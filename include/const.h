
#ifndef _CONST_H
#define _CONST_H

//error constants
#define OK		0
#define ERROR	-1

#define WAVELET_MAX_VALUES 1024		// max number of values used in wavelet
#define WAVELET_USED_VALUES 256		// values actually used in wavelet. MUST BE AN EVEN VALUE.
#define WAVELET_SPOTS 128

//#define SECTIONS_MAX 11

#define DEBUG_LEVEL 0

#define CONFIGS_MAX_PARAMS 128
#define CONFIGS_TOKENS "[]=\n"

#define INTERNAL_CUBE_TOKENS "\n"

#define GLOBAL_NULL_VALUE -99999.25

#define GLOBAL_FILENAME_SIZE 256
#define GLOBAL_LINE_SIZE 256

#define GLOBAL_DUMPTOFILE_CORR 0
#define GLOBAL_DUMPTOFILE_SYNTH 0
#define GLOBAL_DUMPTOFILE_BAI 0
#define GLOBAL_DUMPTOFILE_BCM 0

#define GLOBAL_INIT_SIMULATIONS 1
#define GLOBAL_END_SIMULATIONS 0

#define CONFIGS_FILENAME "./conf/tsi_config.par"

#define INTERNAL_FILENAME_MAXSIZE 1024
#define INTERNAL_ATOI_BUFFER_MAXSIZE 32

#define INTERNAL_DUMP_DOUBLE_VALUE_FORMAT "%.5f\n"

#define BCM_INIT_VALUE 0
#define BAI_INIT_VALUE 3


#endif 


