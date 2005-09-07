
#include <R.h>
#include <Rdefines.h>

/* dummy code a factor where NAs are ignored,
 * i.e. all indicators are FALSE */

SEXP asDummy(SEXP R_x) {

    int n, l, i, j;

    SEXP R_obj, R_tmp;
	
    n = LENGTH(R_x);
    l = LENGTH(GET_LEVELS(R_x));

    if (l == 0)
       return R_NilValue;

    PROTECT(R_obj = NEW_LOGICAL(n*l));

    for (i = 0; i < n*l; i++)		/* this sucks! */
	LOGICAL(R_obj)[i] = FALSE;

    for (i = 0; i < n; i++) {
	j = INTEGER(R_x)[i];
	if (j == NA_INTEGER)
	   continue;
	LOGICAL(R_obj)[i+(j-1)*n] = TRUE;
    }

    PROTECT(R_tmp = NEW_INTEGER(2));

    INTEGER(R_tmp)[0] = n;
    INTEGER(R_tmp)[1] = l;

    SET_DIM(R_obj, R_tmp);
    UNPROTECT(1);
    
    SET_LEVELS(R_obj, duplicate(GET_LEVELS(R_x)));

    UNPROTECT(1);
    
    return R_obj;
}

/**/
