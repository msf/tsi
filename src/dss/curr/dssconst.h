#ifndef DSSCONST_H_
#define DSSCONST_H_

/* The "C" attribute prevents C++ name mangling.
   Remove it if the file type is .c */
/*
#ifdef __cplusplus
	extern "C" 
#endif
*/

/* This source file can be used with the Fortran routine
   being in a DLL or directly in the project. */
/*
#ifdef USEDLL
	__declspec(dllimport)
#endif
*/

/* defines the position of all the parameters in the ParamsArray */
#define NVARI 0
#define IXL 1
#define IYL 2
#define IZL 3
#define IVRL 4
#define IWT 5
#define ISECVR 6
#define TMIN 7
#define TMAX 8
#define ITRANS 9 
#define ISMOOTH 10
#define ISVR 11
#define ISWT 12
#define ZMIN 13
#define ZMAX 14
#define LTAIL 15
#define LTPAR 16
#define UTAIL 17
#define UTPAR 18
#define NTRY 19
#define ICMEAN 20
#define ICVAR 21
#define NX 22
#define XMN 23
#define XSIZ 24
#define NY 25
#define YMN 26
#define YSIZ 27
#define NZ 28
#define ZMN 29
#define ZSIZ 30
#define NOSVALUE 31
#define IMASK 32
#define SEED 33
#define NDMIN 34
#define NDMAX 35
#define NODMAX 36
#define SSTRAT 37
#define MULTS 38
#define NMULTS 39
#define NOCT 40
#define RADIUS 41
#define RADIUS1 42
#define RADIUS2 43
#define SANG 44
#define SANG1 45
#define SANG2 46
#define KTYPE 47
#define COLOCORR 48
#define VARRED 49
#define NVARIL 50
#define ICOLLVM 51
#define VARNUM 52
#define NUGGET 53

#define DSSDLL_TOTAL_PARAMS_NUM 54

/*
extern "C" { 
	int dssdll(float* DSSParams,
		   float* DSSModels,
           double* HardData,
           int* HardDataSize,
           float* BCMData,
           float* BAIData,
           float* MaskData,
           float* OutputData);
}
*/

#endif 


