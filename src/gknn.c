
#include <R.h>
#include <Rdefines.h>

/* gknn.c
 * 
 * generic k-nearest neighbor algorithm which operates on a 
 * user-supplied distance matrix. 
 *
 * implements:
 * 
 * 1) inclusion of tied (equi-distant) kth neighbors (see option
 *    use_all)
 * 2) random selection of tied (equi-distant) kth neighbors (see
 *    option use_all)
 * 3) minimum vote test (otherwise NA is returned for 'doubt')
 * 3) breaking of tied votes
 * 
 * expects as input a cross-distance matrix, a factor of class values, 
 * the number of neighbors to use, the minimum number of votes 
 * (for a definite decision), options for for handling ties in votes 
 * and/or distances (see above), and an option for inclusion of the 
 * proportions of winning votes. note that classe values may be NA
 * but missing distance values are ignored.
 *
 * returns a factor of class predictions and, optionally, a vector
 * of proportions of winning votes as attribute "prob".
 * 
 * ceeboo (2005)
 */

SEXP gknn(SEXP R_x, SEXP R_y, SEXP R_k, SEXP R_l, SEXP R_break_ties,
					SEXP R_use_all, SEXP R_prob) {

    int nr, n, nc, nn, nv;
    int break_ties, use_all;
    int *y, *o, *c;
    int i, j, k, l, m, v;
    
    double *x;

    SEXP R_obj, R_pro, R_str;
    
    nr = INTEGER(GET_DIM(R_x))[0];	    /* number of test instances */
    n  = INTEGER(GET_DIM(R_x))[1];	    /* number of training instances */
    
    if (LENGTH(R_y) != n)
       error("gknn: \"x\" and \"y\" do not conform");
    
    nc = LENGTH(GET_LEVELS(R_y));	    /* number of classes */

    if (nc < 1)
       error("gknn: \"y\" invalid number of levels");

    if (STRING_ELT(GET_LEVELS(R_y), nc-1) == NA_STRING)
       error("gknn: \"y\" invalid level");
    
    y = INTEGER(R_y);			    /* class indexes (R shifted) */

    for (i = 0; i < n; i++)
	if (y[i] == NA_INTEGER || y[i] < 1 || y[i] > nc)
	   error("gknn: \"y\" invalid value");
    
    nn = INTEGER(R_k)[0];		    /* number of neighbors */

    if (nn < 1 || nn > n)
       error("gknn: invalid number of neighbors");
    
    nv = INTEGER(R_l)[0];		    /* number of minimum votes */

    if (nv < 0 || nv > nn)
       error("gknn: invalid minimum number of votes"); 
    
    break_ties = LOGICAL(R_break_ties)[0];  /* tie breaking */
    use_all = LOGICAL(R_use_all)[0];	    /* use all neighbors */
    
    o = Calloc(n, int);			    /* order */
    c = Calloc(nc+1, int);		    /* class counts */
    
    x = Calloc(n, double);		    /* distances */

    PROTECT(R_obj = NEW_INTEGER(nr));
    
    if (LOGICAL(R_prob)[0]) {		    /* return proportions */
       PROTECT(R_pro = NEW_NUMERIC(nr));
       
       setAttrib(R_obj, install("prob"), R_pro);
       UNPROTECT(1);
    }
    else
       R_pro = R_NilValue;
		    
    GetRNGstate();
    
    for (i = 0; i < nr; i++) {
	for (j = 0; j < n; j++) {
	    o[j] = j;
	    x[j] = REAL(R_x)[i+j*nr];	    /* copy distances */
	}

	rsort_with_index(x, o, n);

	for (j = 1; j < nc+1; j++)	    /* R shifted */
	    c[j] = 0;

	k = 0;                              /* invalid class */
	for (j = 0; j < nn; j++) {	    /* count classes */
	    if (ISNAN(x[j]))
	       break;
	    k = y[o[j]];
	    c[k]++;
	}
	if (use_all) {			    /* use tied */
	   while (j < n && x[j] == x[j-1]) {
		 k = y[o[j++]];
		 c[k]++;
	   }
	}
	else {				    /* draw from tied */
	   while (j < n && x[j] == x[j-1])
		 j++;
	   if (j > nn) {
	      l = nn - 1 + (int) (unif_rand() * (j-nn+1));
	      l = y[o[l]];
	      if (l != k) {
		 c[k]--;
		 k = l;
		 c[k]++;
	      }
	   }
	}
	l = 0;				    /* number of ties */
	v = 0;				    /* number of votes */
	m = 0;				    /* max votes */
	for (j = 1; j < nc+1; j++) {
	    v += c[j];
	    if (c[j] > m) {
	       m = c[j];
	       k = j;
	       l = 1;
	    } 
	    else if (l > 0 && c[j] == m) {
	       if (unif_rand() > (double) l/(l+1))
	          k = j;
	       l++;
	    }
	}
	if (R_pro != R_NilValue) {
	   if (v > 0)
	      REAL(R_pro)[i] = (double) m/v;
	   else
	      REAL(R_pro)[i] = NA_REAL;
	}
	if (nv > m) {			    /* below minimum vote */
	   INTEGER(R_obj)[i] = NA_INTEGER;
	}
	else {
	   if (l > 0) {
	      if (break_ties)
	         INTEGER(R_obj)[i] = k;
	      else {
	         if (l > 1)
		    INTEGER(R_obj)[i] = NA_INTEGER;
	         else
	            INTEGER(R_obj)[i] = k;
	      }
	   }
	   else
	      INTEGER(R_obj)[i] = NA_INTEGER;
	}
    }
    Free(o);
    Free(c);
    Free(x);
    
    PutRNGstate();

    SET_LEVELS(R_obj, duplicate(GET_LEVELS(R_y)));

    PROTECT(R_str = NEW_STRING(1));
    SET_STRING_ELT(R_str , 0, mkChar("factor"));
		    
    SET_CLASS(R_obj, R_str);
    UNPROTECT(1);
    
    UNPROTECT(1);
    
    return R_obj;
}

/**/
