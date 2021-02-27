#ifndef _REGDATA_H
#define _REGDATA_H

// data structure
typedef struct regdata_struct {
	int n;
	int p;
	double *w;		// sampling weight
	double *x;		// design matrix (raw and weighted)
	double *wx;
	double *y;		// response vector (raw and weighted)
	double *wy;
} regdata;
#endif
