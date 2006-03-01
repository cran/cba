
#include <R.h>
#include <Rdefines.h>

/* rock.c
 *
 * implements the paper: 
 *
 * S. Guha, R. Rastogi, and K. Shim. Rock: A Robust Clustering Algorithm 
 * for Categorical Attributes. Information Systems, Vol. 25, No. 5, 2000.
 * 
 * The implementation uses a lower triangular matrix representation and 
 * comes in three parts: a function that computes link counts from 
 * distances, another that constructs a cluster solution by merging, and 
 * a function that classifies samples. Implementation of the clustering 
 * problem by separate functions is slightly inefficient but allows for 
 * reuse and further experimentation.
 *
 * The auxiliary functions for computation of "binary" distances are
 * considerably faster than R dist, and the second usage even does not 
 * seem to be available in R (moved to dists.c).
 * 
 * Release: 0.1-1
 * Release date: 2005-06-30
 * 
 * (C) ceeboo 2005
 */


/* 
 * Compute Rock link counts (of the number of common neighbors of
 * any two data points).
 *
 * As input we exepct a lower triangular matrix organized by columns 
 * but with the diagonal omitted which contains the distance values 
 * for any pair of data points. Fixme: currently, the vector has no 
 * attribute indicating its type.
 * 
 * The same data structure is used as the return value.
 *
 * note that we do not test for NA and NaN values because these are by 
 * default large (see the documentation for the behavior of the various
 * sorting functions), so that they (usually) won't pass the threshold 
 * (except the user sets NA, which seems pretty useless to me). 
 */

SEXP rockLink(SEXP R_x, SEXP R_beta) {

    int m, n;
    int i, j, k, kk, l;
    int *v, *p;
    
    double beta;
    double *x; 

    SEXP R_obj;
    
    m = LENGTH(R_x);

    n = 1 + (int) sqrt(2*m);
    
    if (m < 3 || m != n*(n-1)/2)		    /* logical constraint */
       error("calcRLink: invalid vector length");
	       
    x = REAL(R_x);

    beta = REAL(R_beta)[0];
    
    PROTECT(R_obj = NEW_INTEGER(m));
   
    for (l = 0; l < m; l++)
	INTEGER(R_obj)[l] = 0;				/* this sucks!? */
    
    v = Calloc(n, int);
    p = Calloc(n, int);					/* column offset */
    
    for (k = 0; k < n; k++)
	p[k] = k*(n-1)-k*(k+1)/2-1;
    
    for (i = 0; i < n; i++) {
        l = 0;
	for (k = 0; k < i; k++)
	    if (beta >= x[i+p[k]])		    /* omit NA, NaN test */
	       v[l++] = k;
	kk = p[i];
	for (k = i+1; k < n; k++)
	    if (beta >= x[k+kk])		    /* omit NA, NaN test */
	       v[l++] = k;
	for (j = 1; j < l; j++)
	    for (k = 0; k < j; k++) {
		kk = p[v[k]];
		INTEGER(R_obj)[v[j]+kk]++;
	    }
    }
    
    Free(p);
    Free(v);
    
    UNPROTECT(1);
    
    return R_obj;
}

/* 
 * Successively merge two clusters, either unitl the desired 
 * number of clusters is reached, or stop if there are all but 
 * zero link counts left. 
 *
 * returns a list containing a factor with levels labled 
 * contiguously and starting with "1", and a table of cluster 
 * sizes.
 * 
 * This code is optimized for low memory footprint and computation
 * time. However, currently it does not use sparse representations.
 *
 * The search for the maximum merger among n candiates has time 
 * complexity O(n*(n-1)/2) in the worst case and O(n-1) in the 
 * best case.
 *
 * Note that we do not check the link count matrix for NAs or NaNs 
 * because 1) the above function should not return these values, 
 * and 2) the present function is conceptually internal.
 *
 * The neighborhood paramtere is constrained to the interval [0,1) 
 * because inclusion of one results in divison by zero in all 
 * calculations of the merging criterion. For zero link counts
 * this would result in NaNs and undesirable additional checks. 
 * 
 * Fixme: tie breaking!
 */

SEXP rockMerge(SEXP R_x, SEXP R_n, SEXP R_theta, SEXP R_debug) {

    int debug, m, n, nn;
    int i, j, k, l, ii, jj, kk, ll, iii, jjj, kkk;
    int *x, *o, *c, *f, *p, *w;
	
    double y, z;
    double *t, *v;

    char *s;

    SEXP R_obj, R_tmp, R_str, R_dim;
    
    debug = INTEGER(R_debug)[0];

    m = LENGTH(R_x);
    
    n = 1 + (int) sqrt(2*m);			/* number of samples */
    
    if (m < 3 || m != n*(n-1)/2)
       error("rockMerge: invalid vector length");

    nn = INTEGER(R_n)[0];			/* number of clusters */
    
    if (nn < 1)
       error("rockMerge: invalid number of clusters");

    z = REAL(R_theta)[0];			/* neigborhood parameter */
    
    if (z < 0 || z >= 1)
       error("rockMerge: invalid neigborhood parameter");
    
    z = 1 + 2 * (1-z) / (1+z);
    
    x = Calloc(m, int);				/* link counts */
    
    Memcpy(x, INTEGER(R_x), m);

    o = Calloc(n, int);				/* sample index */
    c = Calloc(n, int);				/* cluster index */
    f = Calloc(n, int);				/* cluster size */
    
    p = Calloc(n, int);				/* column offset in dist */
    
    t = Calloc(n+1, double);			/* table of powers */

    v = Calloc(n-1, double);			/* column maximum */
    w = Calloc(n-1, int);			/* row index */
    
    for (k = 0; k < n; k++) {
	o[k] = k;
	c[k] = -1;
	f[k] = 1;
	p[k] = k*(n-1)-k*(k+1)/2-1;
	t[k+1] = pow(k+1, z);
    }
    
    /* find the maximum of a column (in the lower 
     * triangular part if it) and the corresponding
     * row index.
     */
    
    y = t[2]-2*t[1];				    /* initially constant */

    k = 0;
    for (i = 0; i < n-1; i++) {
	v[i] = -1;
	for (j = i+1; j < n; j++) {
	    z = x[k++] / y;
	    if (z > v[i]) {
	       v[i] = z;
	       w[i] = j;
	    }
	}
    }
    if (debug)
       Rprintf(" #cls     clids       sizes     goodness\n");
    
    m = n;
    while (m > nn) {

	z = -1;					    /* find the maximum */
	for (ii = 0; ii < m-1; ii++)
	    if (v[ii] > z) {
	       z = v[ii];
	       i = ii;
	    }

        if (z == 0) 
	   break;
	
	ii = o[i];
	j  = w[i];
	jj = o[j];
	
	if (debug) {
	   Rprintf(" %4i %4i %4i [%4i,%4i] %12.6f", 
		    m, ii, jj, f[ii], f[jj], z);
	   if (f[ii] > 1 && f[jj] > 1)
	      Rprintf("+\n");
	   else
	      Rprintf("\n");
	}
	   
	/* merge the frequencies and link counts; check 
	 * if the new cluster provides a new column maximum; 
	 * this is slightly inefficient in the worst case.
	 */
	
	f[ii] += f[jj];

	for (k = 0; k < i; k++) {
	    kk = o[k];
	    kkk = p[kk];
	    x[ii+kkk] += x[jj+kkk];
	    
	    z = x[ii+kkk] / (t[f[ii]+f[kk]] - 
			     t[f[ii]] - t[f[kk]]);
	    if (z > v[k]) {
	       v[k] = z;
	       w[k] = i;
	    } else 
	       if (w[k] == i || w[k] == j) {
		  v[k] = -1;			    /* be deterministic */
		  for (l = k+1; l < m; l++) {
		      if (l == j)
			 continue;
		      ll = o[l];
		      z = x[ll+kkk] / (t[f[kk]+f[ll]] -
				       t[f[kk]] - t[f[ll]]);
		      if (z > v[k]) {
			 v[k] = z;
			 w[k] = l;
		      }
		  }
	       }
	}
	v[i] = -1;				    /* column changed */
	iii = p[ii];
	for (k = i+1; k < j; k++) {
	    kk = o[k];
	    kkk = p[kk];
	    x[kk+iii] += x[jj+kkk];
	    
	    z = x[kk+iii] / (t[f[ii]+f[kk]] -
			     t[f[ii]] - t[f[kk]]);
	    if (z > v[i]) {
	       v[i] = z;
	       w[i] = k;
	    }
	    if (w[k] == j) {
	       v[k] = -1;
	       for (l = k+1; l < m; l++) {
		   if (l == j)
		      continue;
		   ll = o[l];
		   z = x[ll+kkk] / (t[f[kk]+f[ll]] -
				    t[f[kk]] - t[f[ll]]);
		   if (z > v[k]) {
		      v[k] = z;
		      w[k] = l;
		   }
	       }
	    }
	}
	jjj = p[jj];
	for (k = j+1; k < m; k++) {
	    kk = o[k];
	    x[kk+iii] += x[kk+jjj];
	    
	    z = x[kk+iii] / (t[f[ii]+f[kk]] - 
			     t[f[ii]] - t[f[kk]]);
	    if (z > v[i]) {
	       v[i] = z;
	       w[i] = k;
	    }
	    if (w[k] == j) {
	       v[k] = -1;
	       kkk = p[kk];
	       for (l = k+1; l < m; l++) {
		   if (l == j)
		      continue;
		   ll = o[l];
		   z = x[ll+kkk] / (t[f[kk]+f[ll]] - 
				    t[f[kk]] - t[f[ll]]);
		   if (z > v[k]) {
		      v[k] = z;
		      w[k] = l;
		   }
	       }
	    }
	}

	/* reorganize the indexes of the clusters,
	 * of the rows corresponding to the maxima,
	 * and shrink the vectors. 
	 */

        if (c[ii] == -1)
	   c[ii] = ii;
	if (c[jj] == -1)
	   c[jj] = c[ii];
	else {
	   iii = c[ii];
	   jjj = c[jj];
	   for (k = 0; k < n; k++)
	       if (c[k] == jjj)
		  c[k] = iii;
	}
	
	for (k = 0; k < j; k++)			    /* for clarity here */
	    if (w[k] > j)
	       w[k]--;
	
	for (k = j+1; k < m; k++) {
	    o[k-1] = o[k];
	    v[k-1] = v[k];
	    w[k-1] = w[k]-1;
	}
	m--;
    }
    Free(x);
    
    Free(p);
    Free(t);
    Free(v);
    Free(w);
    
    if (m > nn)
       Rprintf("rockMerge: terminated with %i clusters\n", m);

    PROTECT(R_obj = NEW_LIST(2));

    PROTECT(R_tmp = NEW_INTEGER(n));
    
    /* reorganize the indexes of the 
     * clusters to be contiguous and 
     * to start with one.
     */

    for (k = 0; k < n; k++)
	o[k] = -1;
    m = 0;
    for (k = 0; k < n; k++) {
	if (c[k] == -1)
	   c[k] = k;
	kk = c[k];
	if (o[kk] == -1)
	   o[kk] = ++m;
	INTEGER(R_tmp)[k] = o[kk];
    }

    s = Calloc(m/10+2, char);		/* stringified integers */
    
    PROTECT(R_str = NEW_STRING(m));
    for (j = 0; j < m; j++) {
	sprintf(s,"%i",j+1);
	SET_ELEMENT(R_str, j, mkChar(s));
    }
    Free(s);
    
    SET_LEVELS(R_tmp, R_str);
    UNPROTECT(1);
    
    PROTECT(R_str = NEW_STRING(1));
    SET_ELEMENT(R_str, 0, mkChar("factor"));
	                    
    SET_CLASS(R_tmp, R_str);
    UNPROTECT(1);
    
    SET_ELEMENT(R_obj, 0, R_tmp);
    UNPROTECT(1);

    PROTECT(R_tmp = NEW_INTEGER(m));
  
    for (k = 0; k < n; k++) {
	kk = c[k];
	if (o[kk] != -1) {
	   INTEGER(R_tmp)[o[kk]-1] = f[kk];
	   o[kk] = -1;
	}
    }
    
    Free(o);
    Free(c);
    Free(f);

    PROTECT(R_dim = NEW_INTEGER(1));
    INTEGER(R_dim)[0] = m;
    
    SET_DIM(R_tmp, R_dim);
    UNPROTECT(1);
   
    PROTECT(R_dim = NEW_LIST(1));
    SET_ELEMENT(R_dim, 0, duplicate(GET_LEVELS(VECTOR_ELT(R_obj, 0))));
    
    SET_DIMNAMES(R_tmp, R_dim);
    UNPROTECT(1);
    
    PROTECT(R_str = NEW_STRING(1));
    SET_ELEMENT(R_str, 0, mkChar("table"));
    
    SET_CLASS(R_tmp, R_str);
    UNPROTECT(1);
    
    SET_ELEMENT(R_obj, 1, R_tmp);
    UNPROTECT(1);
    
    UNPROTECT(1);
    
    return R_obj;
}

/* 
 * compute a classification based on a Rock clustering. since we 
 * use a threshold on distances a data point may be assigned to 
 * more than one cluster, or even none.
 *
 * we expect the cluster indexes to be a factor, i.e. to be contiguous 
 * and to start with one. the supplied distances have to be equal or
 * greater than zero. NAs and NaNs are allowed (see the explanation 
 * above).
 *
 * note:
 *
 * 1) ties are broken at random (this may obfuscate that the data 
 *    actually has no structure).
 * 2) points that are not in any neighborhood are assigned the class 
 *    value NA.
 *
 * 
 */

SEXP rockClass(SEXP R_x, SEXP R_l, SEXP R_beta, SEXP R_theta) {

    int nr, nc, nl, na;
    int i, j, h, k;
    int *l, *c, *cf;

    double beta;
    double t, z, y;
    double *n, *x;

    SEXP R_lev, R_obj, R_tmp, R_str, R_dim;
	 
    nr = INTEGER(GET_DIM(R_x))[0];
    nc = INTEGER(GET_DIM(R_x))[1];

    if (LENGTH(R_l) != nc)
       error("rockClass: invalid vector length or number of columns");
 
    R_lev = GET_LEVELS(R_l);
    
    nl = LENGTH(R_lev);
    
    t = REAL(R_theta)[0];

    if (t < 0 || t > 1)
       error("rockMerge: invalid neigborhood parameter");
    
    t = 1 + 2 * (1-t) / (1+t);
    
    l = INTEGER(R_l);			    /* number of levels */

    n = Calloc(nc, double);		    /* expected neighbors */
    
    /* check the validity of the indexes and
     * compute the expected number of neighbors 
     */
    
    for (j = 0; j < nc; j++) {
	i = l[j];
	if (i == NA_INTEGER || i < 1 || i > nl) {
	   Free(n);
	   error("rockClass: invalid cluster index(es)");
	}
	n[i-1]++;
    }
    for (j = 0; j < nl; j++) {
	z = n[j];
	if (z == 0) {				    /* not contiguous */
	   Free(n);
	   error("rockClass: invalid cluster index(es)");
	}
	n[j] = pow(1+z, t);
    }
    
    x = REAL(R_x);				    /* distances */
    
    beta = REAL(R_beta)[0];			    /* threshold */
   
    c = Calloc(nl, int);
    
    PROTECT(R_obj = NEW_LIST(2));
    
    PROTECT(R_tmp = NEW_INTEGER(nr));		    /* class indexes */

    cf = Calloc(nl+1, int);			    /* class frequencies */
    
    GetRNGstate();
    
    for (j = 0; j < nl; j++)
	cf[j] = 0;
    for (i = 0; i < nr; i++) {
	for (j = 0; j < nl; j++)		    /* initialize */
	    c[j] = 0;
	for (j = 0; j < nc; j++)		    /* count neighbors */
	    if (beta >= x[i+j*nr])
	       c[l[j]-1]++;
	k = nl;					    /* include NAs */
	h = 0;					    /* compiler hack */
	z = 0;
	for (j = 0; j < nl; j++) {		    /* determine maximum */
	    y = c[j] / n[j];
	    if (y > z) {
	       z = y;
	       k = j;
	       h = 1;
	    } 
	    else if (h > 0 && y == z) {		    /* break ties */
	       if (unif_rand() > (double) h/(h+1))
		  k = j;
	       h++;
	    }
	}
	cf[k]++;
	INTEGER(R_tmp)[i] = k+1;	    
    }

    PutRNGstate();
    
    Free(n);
    Free(c);

    na = nl+(cf[nl]>0);
    
    PROTECT(R_str = NEW_STRING(na));
    for (j = 0; j < nl; j++)
	SET_ELEMENT(R_str, j, mkChar(CHAR(VECTOR_ELT(R_lev, j))));
    if (na>nl)
       SET_ELEMENT(R_str, j, NA_STRING);
    
    SET_LEVELS(R_tmp, R_str);
    UNPROTECT(1);
    
    PROTECT(R_str = NEW_STRING(1));
    SET_ELEMENT(R_str, 0, mkChar("factor"));
	                    
    SET_CLASS(R_tmp, R_str);
    UNPROTECT(1);
    
    SET_ELEMENT(R_obj, 0, R_tmp);
    UNPROTECT(1);

    PROTECT(R_tmp = NEW_INTEGER(na));
    Memcpy(INTEGER(R_tmp), cf, na);

    PROTECT(R_dim = NEW_INTEGER(1));
    INTEGER(R_dim)[0] = na;

    SET_DIM(R_tmp, R_dim);
    UNPROTECT(1);
    
    PROTECT(R_dim = NEW_LIST(1));
    SET_ELEMENT(R_dim, 0, duplicate(GET_LEVELS(VECTOR_ELT(R_obj, 0)))); 
    
    SET_DIMNAMES(R_tmp, R_dim);
    UNPROTECT(1);

    PROTECT(R_str = NEW_STRING(1));
    SET_ELEMENT(R_str , 0, mkChar("table"));

    SET_CLASS(R_tmp, R_str);
    UNPROTECT(1);		
    
    SET_ELEMENT(R_obj, 1, R_tmp);
    UNPROTECT(1);
    
    UNPROTECT(1);
    
    return R_obj;
}

/**/
