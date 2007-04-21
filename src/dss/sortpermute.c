#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dss.h"


#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE (1)
#define FALSE (0)





/* ----------------------------------------------------------------------- */
/*                      Quickersort Subroutine */
/*                      ********************** */

/* This is a subroutine for sorting a real array in ascending order. This */
/* is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen */
/* in collected algorithms of the ACM. */

/* The method used is that of continually splitting the array into parts */
/* such that all elements of one part are less than all elements of the */
/* other, with a third part in the middle consisting of one element.  An */
/* element with value t is chosen arbitrarily (here we choose the middle */
/* element). i and j give the lower and upper limits of the segment being */
/* split.  After the split a value q will have been found such that */
/* a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then */
/* performs operations on the two segments (i,q-1) and (q+1,j) as follows */
/* The smaller segment is split and the position of the larger segment is */
/* stored in the lt and ut arrays.  If the segment to be split contains */
/* two or fewer elements, it is sorted and another segment is obtained */
/* from the lt and ut arrays.  When no more segments remain, the array */
/* is completely sorted. */

/* INPUT PARAMETERS: */
/*   start,end        start and end index of the array to be sorted */
/*   a            array, a portion of which has to be sorted. */
/*   iperm        0 no other array is permuted. */
/*                1 array b is permuted according to array a */

/* OUTPUT PARAMETERS: */
/*    a      = the array, a portion of which has been sorted. */
/*    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm) */



int sort_permute_float(int start,  int end, float *a, float *b)
           
{
	/* System generated locals */
	int i1;

	/* Local variables */
	int i, j, k, m, p, q;
	float  ta = 0, tb = 0,   /* LPL Init values */
                 xa ;
	int lt[64];
	float xb;
	int ut[64], iring;


	/* The dimensions for lt and ut have to be at least log (base 2) n */


	/* Initialize: */

	/* Parameter adjustments */
	--b;
	--a;

	/* Function Body */
	j = end;
	m = 1;
	i = start;
	iring = 2;

	/* If this segment has more than two elements  we split it */

L10:
	if ((i1 = j - i - 1) < 0) {
		goto L100;
	} else if (i1 == 0) {
		goto L90;
	} else {
		goto L15;
	}

	/* p is the position of an arbitrary element in the segment we choose the */
	/* middle element. Under certain circumstances it may be advantageous */
	/* to choose p at random. */

L15:
	p = (j + i) / 2;
	ta = a[p];
	a[p] = a[i];
	
	tb = b[p];
	b[p] = b[i];

	/* Start at the beginning of the segment, search for k such that a(k)>t */

	q = j;
	k = i;
L20:
	++k;
	if (k > q) {
		goto L60;
	}
	if (a[k] <= ta) {
		goto L20;
	}

	/* Such an element has now been found now search for a q such that a(q)<t */
	/* starting at the end of the segment. */

L30:
	if (a[q] < ta) {
		goto L40;
	}
	--q;
	if (q > k) {
		goto L30;
	}
	goto L50;

	/* a(q) has now been found. we interchange a(q) and a(k) */

L40:
	xa = a[k];
	a[k] = a[q];
	a[q] = xa;

	xb = b[k];
	b[k] = b[q];
	b[q] = xb;

	/* Update q and search for another pair to interchange: */

	--q;
	goto L20;
L50:
	q = k - 1;
L60:

	/* The upwards search has now met the downwards search: */

	a[i] = a[q];
	a[q] = ta;

	b[i] = b[q];
	b[q] = tb;

	/* The segment is now divided in three parts: (i,q-1),(q),(q+1,j) */
	/* store the position of the largest segment in lt and ut */

	if (q << 1 <= i + j) {
		goto L70;
	}
	lt[m - 1] = i;
	ut[m - 1] = q - 1;
	i = q + 1;
	goto L80;
L70:
	lt[m - 1] = q + 1;
	ut[m - 1] = j;
	j = q - 1;

	/* Update m and split the new smaller segment */

L80:
	++m;
	goto L10;

	/* We arrive here if the segment has  two elements we test to see if */
	/* the segment is properly ordered if not, we perform an interchange */

L90:
	if (a[i] <= a[j]) {
		goto L100;
	}
	xa = a[i];
	a[i] = a[j];
	a[j] = xa;
	xb = b[i];
	b[i] = b[j];
	b[j] = xb;

	/* If lt and ut contain more segments to be sorted repeat process: */

L100:
	--m;
	if (m <= 0) {
		goto L110;
	}
	i = lt[m - 1];
	j = ut[m - 1];
	goto L10;
L110:


	return 0;
} /* sortpermute */



/* ----------------------------------------------------------------------- */
/*                      Quickersort Subroutine */
/*                      ********************** */

/* This is a subroutine for sorting a real array in ascending order. This */
/* is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen */
/* in collected algorithms of the ACM. */

/* The method used is that of continually splitting the array into parts */
/* such that all elements of one part are less than all elements of the */
/* other, with a third part in the middle consisting of one element.  An */
/* element with value t is chosen arbitrarily (here we choose the middle */
/* element). i and j give the lower and upper limits of the segment being */
/* split.  After the split a value q will have been found such that */
/* a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then */
/* performs operations on the two segments (i,q-1) and (q+1,j) as follows */
/* The smaller segment is split and the position of the larger segment is */
/* stored in the lt and ut arrays.  If the segment to be split contains */
/* two or fewer elements, it is sorted and another segment is obtained */
/* from the lt and ut arrays.  When no more segments remain, the array */
/* is completely sorted. */

/* INPUT PARAMETERS: */
/*   start,end        start and end index of the array to be sorted */
/*   a            array, a portion of which has to be sorted. */
/*   iperm        0 no other array is permuted. */
/*                1 array b is permuted according to array a */
/* NO EXTERNAL ROUTINES REQUIRED: */
/* ----------------------------------------------------------------------- */



int sort_permute_int(int start,  int end,
           float *a, int *b)
{
	/* System generated locals */
	int i1;

	/* Local variables */
	int i, j, k, m, p, q, tb = 0;
	float  xa, ta = 0;  /* LPL Init values */
	int lt[64];
	int xb;
	int ut[64], iring;


	/* The dimensions for lt and ut have to be at least log (base 2) n */


	/* Initialize: */

	/* Parameter adjustments */
	--b;
	--a;

	/* Function Body */
	j = end;
	m = 1;
	i = start;
	iring = 2;

	/* If this segment has more than two elements  we split it */

L10:
	if ((i1 = j - i - 1) < 0) {
		goto L100;
	} else if (i1 == 0) {
		goto L90;
	}

	/* p is the position of an arbitrary element in the segment we choose the */
	/* middle element. Under certain circumstances it may be advantageous */
	/* to choose p at random. */

	p = (j + i) / 2;
	ta = a[p];
	a[p] = a[i];
	switch (iring) {
		case 1:  goto L21;
		case 2:  goto L19;
	}
L19:
	tb = b[p];
	b[p] = b[i];
L21:

	/* Start at the beginning of the segment, search for k such that a(k)>t */

	q = j;
	k = i;
L20:
	++k;
	if (k > q) {
		goto L60;
	}
	if (a[k] <= ta) {
		goto L20;
	}

	/* Such an element has now been found now search for a q such that a(q)<t */
	/* starting at the end of the segment. */

L30:
	if (a[q] < ta) {
		goto L40;
	}
	--q;
	if (q > k) {
		goto L30;
	}
	goto L50;

	/* a(q) has now been found. we interchange a(q) and a(k) */

L40:
	xa = a[k];
	a[k] = a[q];
	a[q] = xa;
	switch (iring) {
		case 1:  goto L45;
		case 2:  goto L44;
	}
L44:
	xb = b[k];
	b[k] = b[q];
	b[q] = xb;
L45:

	/* Update q and search for another pair to interchange: */

	--q;
	goto L20;
L50:
	q = k - 1;
L60:

	/* The upwards search has now met the downwards search: */

	a[i] = a[q];
	a[q] = ta;
	switch (iring) {
		case 1:  goto L65;
		case 2:  goto L64;
	}
L64:
	b[i] = b[q];
	b[q] = tb;
L65:

	/* The segment is now divided in three parts: (i,q-1),(q),(q+1,j) */
	/* store the position of the largest segment in lt and ut */

	if (q << 1 <= i + j) {
		goto L70;
	}
	lt[m - 1] = i;
	ut[m - 1] = q - 1;
	i = q + 1;
	goto L80;
L70:
	lt[m - 1] = q + 1;
	ut[m - 1] = j;
	j = q - 1;

	/* Update m and split the new smaller segment */

L80:
	++m;
	goto L10;

	/* We arrive here if the segment has  two elements we test to see if */
	/* the segment is properly ordered if not, we perform an interchange */

L90:
	if (a[i] <= a[j]) {
		goto L100;
	}
	xa = a[i];
	a[i] = a[j];
	a[j] = xa;
	switch (iring) {
		case 1:  goto L95;
		case 2:  goto L94;
	}
L94:
	xb = b[i];
	b[i] = b[j];
	b[j] = xb;
L95:

	/* If lt and ut contain more segments to be sorted repeat process: */

L100:
	--m;
	if (m <= 0) {
		goto L110;
	}
	i = lt[m - 1];
	j = ut[m - 1];
	goto L10;
L110:


	return 0;
} /* sortem_ */

