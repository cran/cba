
#include <R.h>
#include <Rdefines.h>

/* cluster_dist
 *
 * cluster an undirected graph as representable by an R dist object,
 * i.e. find all the disconnected subgraphs of graph.
 * 
 * as input we expect R_x the vector storage representation of the 
 * upper/lower triangle of a distance matrix (see dist) and R_beta
 * the distance threshold.
 *
 * returns a factor of cluster labels (integers 1,2, ..., k, with
 * k the number of clusters).
 *
 * NA or NaN distance values are interpreted as no link! this is a 
 * simplification as we do not want to check for the 2^k possible 
 * clusterings given each indeterminate link is actually either above 
 * or below the threshold.
 * 
 * NA or NaN threshold values result in an error as we do not want 
 * to to check for all the possible clusterings given the threshold 
 * assumes a value in the range of the distances.
 * 
 * fixme: 1) can we do this in less than O(n^2) time?
 *	  2) can we use a strict threshold? 
 * 
 * ceeboo 2006
 */

SEXP cluster_dist(SEXP R_x, SEXP R_beta) {
    if (TYPEOF(R_x) != REALSXP)
	error("cluster_dist: 'x' invalid storage type");
    if (TYPEOF(R_beta) != REALSXP)
	error("cluster_dist: 'beta' invalid storage type");
    int i, j, k, l, n, o, na, *c, *b;
    char *s;
    double beta, *x;

    SEXP R_str, R_obj;

    n = (int) sqrt(2 * length(R_x)) + 1;

    if (n < 3 || n * (n - 1) / 2 != length(R_x))
       error("cluster_dist: 'x' invalid length");
   
    beta = REAL(R_beta)[0];	    /* distance threshold */

    if (ISNAN(beta))
       error("cluster_dist: 'beta' NA or NaN");
	       
    PROTECT(R_obj = NEW_INTEGER(n));
    c = INTEGER(R_obj);

    for (i = 0; i < n; i++)
	c[i] = i;
    
    x = REAL(R_x);
    
    k = na = 0;
    for (i = 0; i < n - 1; i++)
	for (j = i + 1; j < n; j++) {
	    if (ISNAN(x[k])) {
	       na++;
	       continue;
	    }
	    if (beta >= x[k++]) {
	       if (c[j] == c[i])
		  continue;
	       if (c[j] == j)
		  c[j] = c[i];
	       else {
		  o = c[j];
		  for (l = 0; l < n; l++)
		      if (c[l] == o)
		         c[l] = c[i];
	       }
	    }
	}
    if (na)
       warning("cluster_dist: found NA (NaN) distance values, different solutions may be possible.");
    
    /* make indexes contiguous */
    
    b = Calloc(n, int);
    
    k = 0;
    for (i = 0; i < n; i++) {
	j = c[i];
	if (b[j] == 0)
	   b[j] = ++k;
	c[i] = b[j];
    }
    Free(b);
    
    /* make return value a factor */
   
    int sn = k/10+2;
    s = Calloc(sn, char);           /* stringified integers */
    
    PROTECT(R_str = NEW_STRING(k));
    for (j = 0; j < k; j++) {
        snprintf(s,sn,"%i",j+1);
        SET_STRING_ELT(R_str, j, mkChar(s));
    }
    Free(s);
 
    SET_LEVELS(R_obj, R_str);
    UNPROTECT(1);

    PROTECT(R_str = NEW_STRING(1));
    SET_STRING_ELT(R_str, 0, mkChar("factor"));
	                            
    SET_CLASS(R_obj, R_str);
    UNPROTECT(1);
		
    /* we are done */
	    
    UNPROTECT(1);
    
    return R_obj;
}



