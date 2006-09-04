
#include <R.h>
#include <Rdefines.h>

/* compute the stress measure based on Moor Neighborhoods, i.e. the 
 * sums of the squared distances of a point to its eight (five at the 
 * margins and three at the corners) adjacent neighbors as defined by 
 * the row and column indexes (or subsets of it).
 *
 * this function counts each edge distance only once! so, if you 
 * prefer the measure from the paper you have to take twice the 
 * value.
 * 
 * note that NAs are omitted. however, the function does not return 
 * NA if there was no legal edge at all.
 */

double stressMoore(double *x, int *r, int *c, int nr, int nc, int nrx) {
        
    double d, v, z;
    int i, j, l, ll, k, kk;
    
    z = 0;
    l = r[0];
    for (i = 0; i < nr-1; i++) {
	ll = r[i+1];
	k  = c[0] * nrx;
	for (j = 0; j < nc-1; j++) {
	    kk = c[j+1] * nrx;    
	    v  = x[l+k];
	    if (!ISNAN(v)) {
	       d = v - x[ll+k];
	       if (!ISNAN(d))
	          z += d * d;
	       d = v - x[ll+kk];
	       if (!ISNAN(d))
	          z += d * d;
	       d = v - x[l+kk];
	       if (!ISNAN(d))
	          z += d * d;
	    }
	    d = x[ll+k] - x[l+kk];
	    k = kk;
	    if (!ISNAN(d))
	       z += d * d;
	}
	d  = x[l+k] - x[ll+k];
	l  = ll;
	if (!ISNAN(d))
	   z += d * d;
	R_CheckUserInterrupt();
    }
    k = c[0] * nrx;
    for (j = 0; j < nc-1; j++) {
	kk = c[j+1] * nrx;
	d  = x[l+k] - x[l+kk];
	k  = kk;
	if (!ISNAN(d))
	   z += d * d;
    }
    
    return z;
}

/* same as above but use a von Neumann neighborhood, i.e. the 
 * neighboring points on the diagonals are excluded.
 */

double stressNeumann(double *x, int *r, int *c, int nr, int nc, int nrx) {
        
    double d, v, z;
    int i, j, l, ll, k, kk;
    
    z = 0;
    l = r[0];
    for (i = 0; i < nr-1; i++) {
	ll = r[i+1];
	k  = c[0] * nrx;
	for (j = 0; j < nc-1; j++) {
	    kk = c[j+1] * nrx;    
	    v  = x[l+k];
	    if (!ISNAN(v)) {
	       d = v - x[ll+k];
	       if (!ISNAN(d))
	          z += d * d;
	       d = v - x[l+kk];
	       if (!ISNAN(d))
	          z += d * d;
	    }
	    k = kk;
	}
	d  = x[l+k] - x[ll+k];
	l  = ll;
	if (!ISNAN(d))
	   z += d * d;
	R_CheckUserInterrupt();
    }
    k = c[0] * nrx;
    for (j = 0; j < nc-1; j++) {
	kk = c[j+1] * nrx;
	d  = x[l+k] - x[l+kk];
	k  = kk;
	if (!ISNAN(d))
	   z += d * d;
    }
    
    return z;
}

/* R wrapper to the stress functions
 */

SEXP stress(SEXP R_x, SEXP R_r, SEXP R_c, SEXP R_type) {

    int nrx, nr, nc;
    int k;
    int *r, *c;
	
    SEXP R_obj;

    PROTECT(R_r = arraySubscript(0, R_r, GET_DIM(R_x), getAttrib, 
				                       (STRING_ELT), R_x));
    PROTECT(R_c = arraySubscript(1, R_c, GET_DIM(R_x), getAttrib, 
						       (STRING_ELT), R_x));
    
    nrx = INTEGER(GET_DIM(R_x))[0];		/* number of rows */
    
    nr = LENGTH(R_r);
    nc = LENGTH(R_c);
    
    /* remap R indexes to C indexes
     * this sucks!
     */
    
    r = Calloc(nr, int);
    c = Calloc(nc, int);

    /* copy and shift indexes */
    
    for (k = 0; k < nr; k++)
	r[k] = INTEGER(R_r)[k]-1;
    for (k = 0; k < nc; k++)
	c[k] = INTEGER(R_c)[k]-1;
   
    PROTECT(R_obj = NEW_NUMERIC(1));

    switch (INTEGER(R_type)[0]) {
	case 1:
	    REAL(R_obj)[0] = stressMoore(REAL(R_x), r, c, nr, nc, nrx);
	    break;
	case 2:
	    REAL(R_obj)[0] = stressNeumann(REAL(R_x), r, c, nr, nc, nrx);
	    break;
	default:
	    Free(r);
	    Free(c);
	    error("stress: type not implemented");
    }
    Free(r);
    Free(c);

    UNPROTECT(3);
       
    return R_obj;
}
/* calculate the Moore distances between all pairs of rows or columns.
 * of a matrix. for a given (fixed) row or column ordering the distances 
 * could be used to search for an optimal column or row ordering using 
 * an alternating scheme.
 *
 * if the calculation are over the rows ncx = 1, otherwise the roles 
 * of rows and columns are swapped and nrx = 1.
 *
 * the caller must provide the result array d and the temporary array t.
 *
 * the distances are arranged in lower triangular column format (compare
 * the R function dist).
 *
 * note that the edge distances are computed only once!
 *
 * (C) ceeboo 2005, 2006
 */

void distMoore(double *x, int *r, int *c, int nr, int nc, int nrx, int ncx, 
		                                  double *d, double *t) {

    double v, w, z;
    int i, ii, j, jj, k, kk, kkk, l;
    
    for (k = 0; k < nr*(nr-1)/2; k++)	    /* initialize distances */
	d[k] = 0;

    for (i = 0; i < nr; i++) {
	z  = 0;
	ii = r[i] * ncx;
	kk = c[0] * nrx;
	for (k = 0; k < nc-1; k++) {
	    kkk = c[k+1] * nrx;
	    w = x[ii+kk] - x[ii+kkk];
	    if (!ISNAN(w))
	       z += w * w;
	    kk = kkk;
	}
	t[i] = z;
	R_CheckUserInterrupt();
    }
    l = 0;
    for (i = 0; i < nr-1; i++) {
	ii = r[i] * ncx;
	for (j = i+1; j < nr; j++) {
	    z  = t[i] + t[j];
	    jj = r[j] * ncx;
	    kk = c[0] * nrx;
	    for (k = 0; k < nc-1; k++) {
		kkk = c[k+1] * nrx;
		v   = x[ii+kk];
		if (!ISNAN(v)) {
		   w = v - x[jj+kk];
		   if (!ISNAN(w))
		      z += w * w;
		   w = v - x[jj+kkk];
		   if (!ISNAN(w))
		      z += w * w;
		}
		w = x[jj+kk] - x[ii+kkk];
		if (!ISNAN(w))
		   z += w * w;
		kk = kkk;
	    }
	    w = x[ii+kk] - x[jj+kk];
	    if (!ISNAN(w))
	       z += w * w;

	    d[l++] = z;
	    R_CheckUserInterrupt();
	}
    }
}

/* calculate the von Neumann distances over the rows or columns of a
 * matrix.
 *
 * compare above.
 */

void distNeumann(double *x, int *r, int *c, int nr, int nc, int nrx, int ncx, 
		                                    double *d, double *t) {

    double w, z;
    int i, ii, j, jj, k, kk, kkk, l;
    
    for (k = 0; k < nr*(nr-1)/2; k++)	    /* initialize distances */
	d[k] = 0;

    for (i = 0; i < nr; i++) {
	z  = 0;
	ii = r[i] * ncx;
	kk = c[0] * nrx;
	for (k = 0; k < nc-1; k++) {
	    kkk = c[k+1] * nrx;
	    w = x[ii+kk] - x[ii+kkk];
	    if (!ISNAN(w))
	       z += w * w;
	    kk = kkk;
	}
	t[i] = z;
	R_CheckUserInterrupt();
    }
    l = 0;
    for (i = 0; i < nr-1; i++) {
	ii = r[i] * ncx;
	for (j = i+1; j < nr; j++) {
	    z  = t[i] + t[j];
	    jj = r[j] * ncx;
	    for (k = 0; k < nc-1; k++) {
		kk = c[k] * nrx;
		w = x[ii+kk]- x[jj+kk];
		if (!ISNAN(w))
		   z += w * w;
	    }
	    kk = c[k] * nrx;
	    w  = x[ii+kk] - x[jj+kk];
	    if (!ISNAN(w))
	       z += w * w;
	    
	    d[l++] = z;
	    R_CheckUserInterrupt();
	}
    }
}

/* R wrapper
 */

SEXP stress_dist(SEXP R_x, SEXP R_r, SEXP R_c, SEXP R_bycol, SEXP R_type) {

    int nrx, nr, nc;
    int k;
    int *r, *c;
	
    double *d, *t;
    
    SEXP R_obj = R_NilValue;	/* compiler hack */

    PROTECT(R_r = arraySubscript(0, R_r, GET_DIM(R_x), getAttrib, 
                                                       (STRING_ELT), R_x));
    PROTECT(R_c = arraySubscript(1, R_c, GET_DIM(R_x), getAttrib, 
                                                       (STRING_ELT), R_x));
    
    nrx = INTEGER(GET_DIM(R_x))[0];		/* number of rows */
    
    nr = LENGTH(R_r);
    nc = LENGTH(R_c);
    
    /* remap R indexes to C indexes
     * this sucks!
     */
    
    r = Calloc(nr, int);
    c = Calloc(nc, int);
    
    /* copy and shift indexes */
    
    for (k = 0; k < nr; k++)
        r[k] = INTEGER(R_r)[k]-1;
    for (k = 0; k < nc; k++)
        c[k] = INTEGER(R_c)[k]-1;
   
    switch(LOGICAL(R_bycol)[0]) {
	case 0:
	    PROTECT(R_obj = NEW_NUMERIC(nr*(nr-1)/2));

	    d = REAL(R_obj);
	    t = Calloc(nr, double);
	    
	    switch(INTEGER(R_type)[0]) {
		case 1:
	            distMoore(REAL(R_x), r, c, nr, nc, nrx, 1, d, t);
		    break;
		case 2:
	            distNeumann(REAL(R_x), r, c, nr, nc, nrx, 1, d, t);
		    break;
		default:
		    Free(r);
		    Free(c);
		    Free(t);
		    error("stress_dist: \"type\" not implemented");
	    }
	    Free(t);
	    break;
	case 1:
	    PROTECT(R_obj = NEW_NUMERIC(nc*(nc-1)/2));

	    d = REAL(R_obj);
	    t = Calloc(nc, double);
	    
	    switch(INTEGER(R_type)[0]) {
		case 1:
		    distMoore(REAL(R_x), c, r, nc, nr, 1, nrx, d, t);
		    break;
		case 2:
		    distNeumann(REAL(R_x), c, r, nc, nr, 1, nrx, d, t);
		    break;
		default:
		    Free(r);
		    Free(c);
		    Free(t);
		    error("stress_dist: type not implemented");
	    }
	    Free(t);
	    break;
	default:
	    Free(r);
	    Free(c);
	    error("stress_dist: \"bycol\" invalid");
    }
    Free(r);
    Free(c);

    UNPROTECT(3);

    return R_obj;
}

/* in ordering problems we want find a path of minimum length through
 * a distance graph. this is a TSP problem with a dummy city that is
 * equally distant from all other cities, i.e. the length of the path
 * i -> 0 -> j, closing the tour is irrelevant. however, the length of
 * the leg i -> j is greater or equal than any remaining leg k -> l on
 * the optimum tour.
 *
 * see: Climer, S. and Zhang, W. (2006) Rearrangement Clustering:
 *      Pitfalls, Remedies, and Applications. Journal of Machine
 *      Learning Research 7, pp. 919-943.
 *
 * orderTSP implements a greedy heuristic that exchanges two edges 
 * immediately if this improves the tour length and stops if no further
 * improvement (over all combinations of edges) is possible. exchanging
 * edges amounts to reversing subpaths.
 *
 * the time complexity is O(n^2) with n the number of cities.
 *
 * note: the algorithm could easily be extended to a simulated 
 *	 annealing algorithm. the code is slightly optimized.
*/

SEXP orderTSP(SEXP x, SEXP t) {
    if (TYPEOF(x) != REALSXP)
	error("'x' invalid storage type");
    if (TYPEOF(t) != INTSXP)
	error("'t' invalid storage type");
    int i, n, f = 0;

    n = 1 + (int) sqrt(2*LENGTH(x));
    if (LENGTH(x) != n*(n-1)/2)
	error("'x' invalid length");
    if (LENGTH(t) != n)
	error("'t' invalid length");

    for (i = 0; i < n; i++)
	if (INTEGER(t)[i] < 1 || INTEGER(t)[i] > n)
	    error("'t' invalid");
    
    PROTECT(t = duplicate(t)); 
    do {
	int i, j, k = 0, l = 0, c1, c2, c3, c4 = n-1;
	double e23, e13, e12, e34, e24, e31, e41;
	
	f = 0;
	c1 = INTEGER(t)[0]-1;
	for (i = 1; i < n-1; i++) {
	    c2 = INTEGER(t)[i]-1;
	    c3 = INTEGER(t)[i+1]-1;
	    if (c2 > c3)
		e23 = REAL(x)[c2+c3*(n-1)-c3*(c3+1)/2-1];
	    else
		e23 = REAL(x)[c3+c2*(n-1)-c2*(c2+1)/2-1];
	    if (c1 > c3)
		e13 = REAL(x)[c1+c3*(n-1)-c3*(c3+1)/2-1];
	    else
		e13 = REAL(x)[c3+c1*(n-1)-c1*(c1+1)/2-1];
	    if (e23 > e13) {
		f++;
		for (k = 0; k < (i+1)/2; k++) {
		    l = INTEGER(t)[i-k];
		    INTEGER(t)[i-k] = INTEGER(t)[k];
		    INTEGER(t)[k] = l;
		}
		c1 = INTEGER(t)[0]-1;
	    }
	}
	for (i = 0; i < n-3; i++) {
	    c1 = INTEGER(t)[i]-1;
	    c2 = INTEGER(t)[i+1]-1;
	    if (c1 > c2)
		e12 = REAL(x)[c1+c2*(n-1)-c2*(c2+1)/2-1];
	    else
		e12 = REAL(x)[c2+c1*(n-1)-c1*(c1+1)/2-1];
	    for (j = i+2; j < n-1; j++) {
		c3 = INTEGER(t)[j]-1;
		c4 = INTEGER(t)[j+1]-1;
		if (c3 > c4)
		    e34 = REAL(x)[c3+c4*(n-1)-c4*(c4+1)/2-1];
		else
		    e34 = REAL(x)[c4+c3*(n-1)-c3*(c3+1)/2-1];
		if (c2 > c4)
		    e24 = REAL(x)[c2+c4*(n-1)-c4*(c4+1)/2-1];
		else
		    e24 = REAL(x)[c4+c2*(n-1)-c2*(c2+1)/2-1];
		if (c3 > c1)
		    e31 = REAL(x)[c3+c1*(n-1)-c1*(c1+1)/2-1];
		else
		    e31 = REAL(x)[c1+c3*(n-1)-c3*(c3+1)/2-1];

		if (e12+e34 > e24+e31) {
		    f++;
		    for (k = 0; k < (j-i)/2; k++) {
			l = INTEGER(t)[j-k];
			INTEGER(t)[j-k] = INTEGER(t)[i+1+k];
			INTEGER(t)[i+1+k] = l;
		    }
		    c2 = INTEGER(t)[i+1]-1;
		    if (c1 > c2)
			e12 = REAL(x)[c1+c2*(n-1)-c2*(c2+1)/2-1];
		    else
			e12 = REAL(x)[c2+c1*(n-1)-c1*(c1+1)/2-1];
		}
	    }
	    if (c4 > c1)
		e41 = REAL(x)[c4+c1*(n-1)-c1*(c1+1)/2-1];
	    else
		e41 = REAL(x)[c1+c4*(n-1)-c4*(c4+1)/2-1];
	    if (e12 > e41) {
		f++;
		for (k = 0; k < (j-i)/2; k++) {
		    l = INTEGER(t)[j-k];
		    INTEGER(t)[j-k] = INTEGER(t)[i+1+k];
		    INTEGER(t)[i+1+k] = l;
		}
	    }
	    R_CheckUserInterrupt();
	}
    } while (f);
   
    UNPROTECT(1);
    return t;
}

/**/
