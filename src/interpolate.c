

#include <R.h>
#include <Rdefines.h>

/* interpolate a logical matrix to a lower resolution.
 * 
 * note 1) that we currently use the full storage representation 
 * of a binary matrix and 2) that some rows and/or columns at the 
 * lower and left margins may get cut off
 * 
 * (C) ceeboo 2005
 */

SEXP lminter(SEXP R_x, SEXP R_block_size, SEXP R_nbin) {

    int nr, nc, np, zr, zc;
    int i, j;

    int *x, *z;
    
    SEXP R_obj, R_dim;

    nr = INTEGER(GET_DIM(R_x))[0];	    /* number of rows */
    nc = INTEGER(GET_DIM(R_x))[1];	    /* number of columns */

    x = LOGICAL(R_x);
    
    np = INTEGER(R_block_size)[0];	    /* number of pixels */
    zr = nr / np;			    /* reduced number of rows */
    zc = nc / np;			    /* reduced number of columns */
    
    PROTECT(R_obj = NEW_INTEGER(zr * zc));
    
    z = INTEGER(R_obj);
    
    for (j = 0; j < zr * zc; j++)	    /* this sucks! */
	z[j] = 0;	
    
    for (j = 0; j < zc * np; j++)
	for (i = 0; i < zr * np; i++)
	    z[i / np + (j / np) * zr] += x[i + j * nr];

    i =  INTEGER(R_nbin)[0];		    /* number of bins  */
    if (i < 0 || i > np)
       error("lminter: invalid number of bins");

    if (i == 0) {			    /* majority */
       i = np * np / 2 + 1;
       for (j = 0; j < zr * zc; j++)
	   z[j] /= i;
    }
    else {				    /* bins */
       i = np * np / i;
       for (j = 0; j < zr * zc; j++)
	   z[j] = ceil((double) z[j] / i);
    }
    PROTECT(R_dim = NEW_INTEGER(2));
    
    INTEGER(R_dim)[0] = zr;
    INTEGER(R_dim)[1] = zc;

    SET_DIM(R_obj, R_dim);
    
    UNPROTECT(2);
    
    return R_obj;
}

/**/
