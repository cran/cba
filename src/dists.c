
#include <R.h>
#include <Rdefines.h>

/* reimplementation of common distance functions so that not only
 * the auto-distances of a matrix can be computed but also the
 * cross-distances between two matrices.
 *
 * the implementation of both behavior in one function incurrs only
 * a minor inefficiency. for extensibility, the implementation uses 
 * a separate function for each distance measure.
 * 
 * fixed 'bug' due to the copy on write mechanism that was sneaked
 * in by R 2.3.0 or 2.3.1. this is a tell tale reminder that we
 * should not rely on anything in R :-( this, of course, makes
 * things less structured, maintainable and, efficient.
 *
 * fixme: 1) conceptually we should return proper dist objects here
 *        instead of piecing it together on the R level. however,
 *        see remark above. 2) complex data are not supported.
 * 
 * ceeboo 2005, 2006
 */

/* calculate Minkowsky distances (including one-norm) between the 
 * rows of two matrices. ignores NAs and NaNs (cf. dist). other than 
 * in this implementation an NA row/column gets not removed from the 
 * result matrix).
 */

SEXP pdist(SEXP R_x, SEXP R_y, SEXP R_p) {

    int nc, nx, ny;
    int i, j, k, l, m, as_matrix = 0;

    double p, d, q, z;
    double *x, *y;
    
    SEXP R_obj;
    
    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
    else
       as_matrix = 1;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("pdist: invalid number of columns");

    p = REAL(R_p)[0];

    if (p < 0 || p == R_PosInf)
       error("pdist: invalid parameter");

    if (p == 0)					    /* one norm */
       p = 1;
    
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    if (!as_matrix && R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = allocMatrix(REALSXP, nx, ny));
		    
  
    q = 1 / p;
 
    m = 0;
    for (j = 0; j < ny; j++) {
	if (!as_matrix && R_x == R_y)
	   i = j + 1;
        else
	   i = 0;
	for (i = i; i < nx; i++) {
	    l = 0;
	    z = 0;
	    for (k = 0; k < nc; k++) {
		if (ISNAN(x[i+k*nx]) || ISNAN(y[j+k*ny]))
		   continue;
		d = fabs(x[i+k*nx]-y[j+k*ny]);
		if (ISNAN(d))
		   continue;
		l++;
		z += pow(d, p);
	    }
	    if (l > 0)
	       REAL(R_obj)[m++] = pow(z, q);
	    else
	       REAL(R_obj)[m++] = NA_REAL;
	}
	R_CheckUserInterrupt();
    }
    
    UNPROTECT(1);

    return R_obj;			    
}

/* calculate Maximum distances (cf. dist) 
 */

SEXP mdist(SEXP R_x, SEXP R_y) {

    int nc, nx, ny; 
    int i, j, k, l, m, as_matrix = 0;
    
    double d, z;
    double *x, *y;
  
    SEXP R_obj;

    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue) 
       R_y = R_x;
    else
       as_matrix = 1;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("mdist: invalid number of columns");
	
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    if (!as_matrix && R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = allocMatrix(REALSXP, nx, ny));
		    
    m = 0;
    for (j = 0; j < ny; j++) {
	if (!as_matrix && R_x == R_y)
	   i = j + 1;
	else
	   i = 0;
	for (i = i; i < nx; i++) {
	    l = 0;
	    z = R_NegInf;
	    for (k = 0; k < nc; k++) {
		if (ISNAN(x[i+k*nx]) || ISNAN(y[j+k*ny]))
		   continue;
		d = fabs(x[i+k*nx]-y[j+k*ny]);
		if (ISNAN(d))
		   continue;
		l++;
		if (d > z) 
		   z = d;
	    }
	    if (l > 0)
	       REAL(R_obj)[m++] = z;
	    else
	       REAL(R_obj)[m++] = NA_REAL;
	}
	R_CheckUserInterrupt();
    }
       
    UNPROTECT(1);

    return R_obj;			    
}

/* calculate Canberra distances (cf. dist)
 */

SEXP cdist(SEXP R_x, SEXP R_y) {
    
    int nc, nx, ny; 
    int i, j, k, l, m, as_matrix = 0;
    
    double d, z;
    double *x, *y;

    SEXP R_obj;
  
    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
    else
       as_matrix = 1;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("cdist: invalid number of columns");
	
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    if (!as_matrix && R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = allocMatrix(REALSXP, nx, ny));
		    
    m = 0; 
    for (j = 0; j < ny; j++) {
	if (!as_matrix && R_x == R_y)
	   i = j + 1;
	else
	   i = 0;
	for (i = i; i < nx; i++) {
	    l = 0;
	    z = 0;
	    for (k = 0; k < nc; k++) {
	        if (ISNAN(x[i+k*nx]) || ISNAN(y[j+k*ny]))
		   continue;
		d = fabs(x[i+k*nx]+y[j+k*ny]);
		if (ISNAN(d))
		   continue;
		l++;
		z += fabs(x[i+k*nx]-y[j+k*ny]) / d;
	    }
	    if (l > 0)
	       REAL(R_obj)[m++] = z;
	    else
	       REAL(R_obj)[m++] = NA_REAL;
	}
	R_CheckUserInterrupt();
    }

    UNPROTECT(1);

    return R_obj;
}

/* 
 * compute binary distances (one minus the Jaccard coefficient)
 * between all pairs of rows from two matrices. in the case the
 * second argument is NULL (R_NilValue) the usual auto-distances 
 * are computed.
 * 
 * returns either a full matrix, or a lower triangular matrix 
 * organized by columns and with the diagonal omitted (cf. dist and 
 * the packages Matrix).
 * 
 * combining both behaviors in one function implies only a slight
 * computational inefficiency.
 *
 * the caller needs to ensure that there are no NAs or NaNs because 
 * zero is already interpreted as no information available. if you want 
 * a different behavior use extended binary.
 *
 * we define the distance of two zero vectors to be zero instead of NA.
 * this is inconsitent with dummy coding of NAs but compatible with dist.
 * 
 * ceeboo 2005
 */

SEXP bdist(SEXP R_x, SEXP R_y) {

    int nc, nx, ny; 
    int i, j, k, l, t, as_matrix = 0;
    
    int *x, *y, *s;

    double z;
 
    SEXP R_obj;
    
    nc = INTEGER(GET_DIM(R_x))[1];
    
    if (R_y == R_NilValue)
       R_y = R_x;
    else
       as_matrix = 1;
    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("bdist: invalid number of columns");

    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = INTEGER(R_x);
    y = INTEGER(R_y);
    
    s = Calloc(nx, int);
    
    if (!as_matrix && R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = allocMatrix(REALSXP, nx, ny));
   
    for (i = 0; i < nx; i++) {
	t = 0;
	for (k = 0; k < nc; k++) {
	    if (x[i+k*nx] == NA_LOGICAL)
	       continue;
	    t += x[i+k*nx];
	}
	s[i] = t;	
    }
    l = 0;
    for (j = 0; j < ny; j++) {
	if (as_matrix || R_x != R_y) {
	   t = 0;
	   for (k = 0; k < nc; k++) {
	       if (y[j+k*ny] == NA_LOGICAL)
		  continue;
	       t += y[j+k*ny];
	   }
	   i = 0;
	}
	else {
	   t = s[j];
	   i = j + 1;
	}
	for (i = i; i < nx; i++) {
	    z = 0;
	    for (k = 0; k < nc; k++) {
		if (x[i+k*nx] == NA_LOGICAL || y[j+k*ny] == NA_LOGICAL)
		   continue;
		z += x[i+k*nx] & y[j+k*ny];
	    }
	    z = 1 -  z / (t + s[i] - z);
	    if (ISNAN(z))			    /* division by zero */
	       REAL(R_obj)[l++] = 0;		    /* but be compatible */
	    else
	       REAL(R_obj)[l++] = z;
	}
	R_CheckUserInterrupt();
    }
   
    Free(s);
   
    UNPROTECT(1);

    return R_obj;
}

/* calculate extended binary distances (one minus extended Jaccard),
 * i.e. the squared Euclidean distance divided by the squared
 * Euclidean distance minus the scalar product. 
 *
 * fixme: the policy of mapping NaNs to NA is not very informative.
 *
 */

SEXP ebdist(SEXP R_x, SEXP R_y) {
    
    int nc, nx, ny; 
    int i, j, k, l, m, as_matrix = 0;
    
    double t, z;
    double *x, *y, *s;
 
    SEXP R_obj;
    
    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
    else
       as_matrix = 1;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("ebdist: invalid number of columns");
	
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    s = Calloc(nx, double);
    
    if (!as_matrix && R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = allocMatrix(REALSXP, nx, ny));
		    
    for (i = 0; i < nx; i++) {
	z = 0;
	l = 0;
	for (k = 0; k < nc; k++) {
	    if (ISNAN(x[i+k*nx]))
	       continue;
	    l++;
	    z+= pow(x[i+k*nx], 2);
	}
	if (l > 0)
	   s[i] = z;
	else
	   s[i] = NA_REAL;
    }
    m = 0; 
    for (j = 0; j < ny; j++) {
	if (!as_matrix && R_x == R_y) {
	   t = s[j];
	   i = j + 1;
	}
	else {
	   z = 0;
	   l = 0;
	   for (k = 0; k < nc; k++) {
	       if (ISNAN(y[j+k*ny]))
		  continue;
	       l++;
	       z+= pow(y[j+k*ny], 2);
	   }
	   if (l > 0)
	      t = z;
	   else
	      t = NA_REAL;
	   i = 0;
	}
	for (i = i; i < nx; i++) {
	    if (ISNAN(t) || ISNAN(s[i])) {
	       REAL(R_obj)[m++] = NA_REAL;
	       continue;
	    }
	    l = 0;
	    z = 0;
	    for (k = 0; k < nc; k++) {
	        if (ISNAN(x[i+k*nx]) || ISNAN(y[j+k*ny]))
		   continue;
		l++;
		z+= x[i+k*nx]*y[j+k*ny];
	    }
	    if (l > 0) {
	       z = 1 - z / (t + s[i] - z);
	       if (ISNAN(z))
		  REAL(R_obj)[m++] = 0;	    /* be compatible */
	       else
		  REAL(R_obj)[m++] = z;
	    }
	    else
	       REAL(R_obj)[m++] = NA_REAL;
	}
	R_CheckUserInterrupt();
    }

    Free(s);
    
    UNPROTECT(1);

    return R_obj;
}

/* calculate the fuzzy generalization of binary Jaccard distances 
 * as proposed by Kurt Hornik, i.e. the ratio of the sums of
 * the component (parallel) minima and maxima, respectively.
 *
 * note that the measure is only reasonably defined for positive
 * real-valued components, e.g. from the interval [0,1].
 */

SEXP fbdist(SEXP R_x, SEXP R_y) {
    
    int nc, nx, ny; 
    int i, j, k, l, m, as_matrix = 0;
    
    double x1, x2, z1, z2;
    double *x, *y;
 
    SEXP R_obj;
	        
    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
    else
       as_matrix = 1;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("fbdist: invalid number of columns");
	
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    if (!as_matrix && R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = allocMatrix(REALSXP, nx, ny));
		    
    m = 0; 
    for (j = 0; j < ny; j++) {
	if (!as_matrix && R_x == R_y) 
	   i = j + 1;
	else 
	   i = 0;
	for (i = i; i < nx; i++) {
	    l = 0;
	    z1 = 0;
	    z2 = 0;
	    for (k = 0; k < nc; k++) {
		x1 = x[i+k*nx];
		x2 = y[j+k*ny];
	        if (ISNAN(x1) || ISNAN(x2))
		   continue;
		l++;
		if (x1 > x2) {
		   z1 += x2;
		   z2 += x1;
		}
		else {
		   z1 += x1;
		   z2 += x2;
		}
	    }
	    if (l > 0) {
	       if (x2 == 0)
		  REAL(R_obj)[m++] = 1;
	       else
		  REAL(R_obj)[m++] = 1 - z1 / z2;
	    }
	    else
	       REAL(R_obj)[m++] = NA_REAL;
	}
	R_CheckUserInterrupt();
    }

    UNPROTECT(1);

    return R_obj;
}
/* calculate angular distances (one minus cosine similarity) */

SEXP adist(SEXP R_x, SEXP R_y) {
    
    int nc, nx, ny; 
    int i, j, k, l, m, as_matrix = 0;
    
    double t, z;
    double *x, *y, *s;
 
    SEXP R_obj;
	        
    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
    else
       as_matrix = 1;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("adist: invalid number of columns");
	
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    s = Calloc(nx, double);
    
    if (!as_matrix && R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = allocMatrix(REALSXP, nx, ny));
		    
    #define EPS 1E-16

    for (i = 0; i < nx; i++) {
	z = 0;
	l = 0;
	for (k = 0; k < nc; k++) {
	    if (ISNAN(x[i+k*nx]))
	       continue;
	    l++;
	    z+= pow(x[i+k*nx], 2);
	}
	if (l > 0)
	   s[i] = sqrt(z);
	else
	   s[i] = NA_REAL;
    }
    m = 0; 
    for (j = 0; j < ny; j++) {
	if (!as_matrix && R_x == R_y) {
	   t = s[j];
	   i = j + 1;
	}
	else {
	   z = 0;
	   l = 0;
	   for (k = 0; k < nc; k++) {
	       if (ISNAN(y[j+k*ny]))
		  continue;
	       l++;
	       z+= pow(y[j+k*ny], 2);
	   }
	   if (l > 0)
	      t = sqrt(z);
	   else
	      t = NA_REAL;
	   i = 0;
	}
	for (i = i; i < nx; i++) {
	    if (ISNAN(t) || ISNAN(s[i])) {
	       REAL(R_obj)[m++] = NA_REAL;
	       continue;
	    }
	    l = 0;
	    z = 0;
	    for (k = 0; k < nc; k++) {
	        if (ISNAN(x[i+k*nx]) || ISNAN(y[j+k*ny]))
		   continue;
		l++;
		z+= x[i+k*nx]*y[j+k*ny];
	    }
	    if (l > 0) {
	       z =  1 - z / t / s[i];
	       if (ISNAN(z)) {
		  if (t < EPS && s[i] < EPS)
		     REAL(R_obj)[m++] = 0;	    /* be compatible */
		  else
		     REAL(R_obj)[m++] = 1;
	       }
	       else
		  REAL(R_obj)[m++] = z;
	    }
	    else
	       REAL(R_obj)[m++] = NA_REAL;
	}
	R_CheckUserInterrupt();
    }

    Free(s);
    
    UNPROTECT(1);

    return R_obj;
}

/* subset a dist object. in order to preserve symmetry 
 * we allow only one subset index. for the subsripting
 * to work we need a proper array (see last argument).
 * 
 * note that subscriptArray returns indexes that start
 * at 1 instead of zero. Yuck!
 * 
 * ceeboo 2005, 2006
 */

SEXP subset_dist(SEXP R_x, SEXP R_s, SEXP R_l) {
    if (TYPEOF(R_x) != REALSXP)
	error("invalid data parameter");
    int n, m;
    int i, ii, j, jj, k;

    int *s;

    double *x, *z;
	
    SEXP R_obj;
    
    PROTECT(R_s = arraySubscript(0, R_s, GET_DIM(R_l), getAttrib, 
						       (STRING_ELT), R_l));
    if (TYPEOF(R_s) != INTSXP)
	error("arraySubscript returned unsupported data type");

    n = 1 + (int) sqrt(2*LENGTH(R_x));  /* number of elements */
    
    if (LENGTH(R_x) != n*(n-1)/2)
	error("invalid length");
    
    m = LENGTH(R_s);			/* size of subset */
    
    if (m < 2)
	error("invalid subscripts");
    for (k = 0; k < m; k++)
	if (INTEGER(R_s)[k] == NA_INTEGER)
	    error("invalid subscripts");
		
    s = INTEGER(R_s);
    x = REAL(R_x);

    PROTECT(R_obj = NEW_NUMERIC(m*(m-1)/2));

    z = REAL(R_obj);
    
    k = 0;
    for (i = 0; i < m-1; i++) {
	ii = s[i]-1;
	for (j = i+1; j < m; j++) {
	    if (s[i] == s[j])
	       z[k++] = 0;              /* by definition */
	    else {
	       jj = s[j]-1;
	       if (s[i] > s[j])
		  z[k++] = x[ii+jj*(n-1)-jj*(jj+1)/2-1];
	       else
		  z[k++] = x[jj+ii*(n-1)-ii*(ii+1)/2-1];
	    }
	}
    }
    
    /* it is more secure to do this here */

    if (!isNull((R_l = GetRowNames(GET_DIMNAMES(R_l))))) {
       SEXP R_str;
       
       PROTECT(R_str = NEW_STRING(m));
       for (k = 0; k < m; k++)
	   SET_STRING_ELT(R_str, k, duplicate(STRING_ELT(R_l, s[k]-1)));

       setAttrib(R_obj, install("Labels"), R_str);
       UNPROTECT(1);
    }
    setAttrib(R_obj, install("Size"), ScalarInteger(m));

    UNPROTECT(2);

    return R_obj;
}

/* compute the rowSums for an R dist object. due to
 * symmetry this is equivalent to colSums. rowMeans
 * are not implemented as this can be trivially
 * obtained from the values of rowSums.
 * 
 * na_rm implements the usual meaning of omitting 
 * NA and NaN values.
 *
 * ceeboo 2006
 */

SEXP rowSums_dist(SEXP x, SEXP na_rm) {
    if (TYPEOF(x) != REALSXP)
        error("invalid data parameter");
    if (TYPEOF(na_rm) != LGLSXP)
        error("invalid option paramter");
    int i, j, k, n;
    SEXP r;

    n = 1 + (int) sqrt(2*LENGTH(x));
    
    if (LENGTH(x) != n*(n-1)/2)
        error("invalid length");
    
    PROTECT(r = allocVector(REALSXP, n));
    
    for (i = 0; i < n; i++)
	REAL(r)[i] = 0;
    
    k = 0;
    for (i = 0; i < n-1; i++) 
        for (j = i+1; j < n; j++) {
            double z = REAL(x)[k++];
            if (!R_FINITE(z)) {
               if (ISNAN(z)) {
                  if (LOGICAL(na_rm)[0] == TRUE)
                     continue;
                  if (ISNA(z))
                     REAL(r)[i] = REAL(r)[j] = NA_REAL;
                  else
                     REAL(r)[i] = REAL(r)[j] = R_NaN;
                } else
                  REAL(r)[i] = REAL(r)[j] = z;
               break;
            }
            REAL(r)[i] += z;
	    REAL(r)[j] += z;
        }
    
    UNPROTECT(1);
    return r;
}

/* compute auto- or cross-distances with a user-supplied
 * function given matrix data. this is experimental code.
 *
 * fixme: eval environment?
 * 
 * ceeboo 2006
 */

SEXP apply_dist(SEXP p) {
    int i, j, k, l, n, nx, ny, as_matrix = 0;
    SEXP r, c, tx, ty;
    SEXP x, y, f;

    p = CDR(p);
    if (length(p) < 3)
	error("invalid number of arguments");
    x = CAR(p); y = CADR(p);
    if (!isMatrix(x) || (!isNull(y) && !isMatrix(y)))
	error("invalid data parameter(s)");
    p = CDDR(p); f = CAR(p); 
    if (!isFunction(f))
	error("invalid function parameter");
    p = CDR(p);

    if (isNull(y))
	y = x;	
    else
       as_matrix = 1;
    
    if ((n = INTEGER(GET_DIM(x))[1]) != INTEGER(GET_DIM(y))[1])
	error("data parameters do not conform");
    
    nx = INTEGER(GET_DIM(x))[0];
    ny = INTEGER(GET_DIM(y))[0];

    if (!as_matrix && x == y)
	PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
    else
	PROTECT(r = allocMatrix(REALSXP, nx, ny));
    
    PROTECT(tx = allocVector(REALSXP, n)); 
    PROTECT(ty = allocVector(REALSXP, n));

    c = LCONS(f, LCONS(tx, LCONS(ty, p)));
    
    l = 0;
    for (j = 0; j < ny; j++) {
	for (k = 0; k < n; k++)
	    REAL(ty)[k] = REAL(y)[j+k*ny];
	for (i = ((!as_matrix && x == y) ? j+1 : 0); i < nx; i++) {
	    for (k = 0; k < n; k++)
		REAL(tx)[k] = REAL(x)[i+k*nx];
	    SEXP s = eval(c, R_GlobalEnv);
	    if (TYPEOF(s) != REALSXP || LENGTH(s) != 1)
		error("invalid return value");
	    REAL(r)[l++] = REAL(s)[0];
	}			
	R_CheckUserInterrupt();
    }

    UNPROTECT(3);
    return r;
}

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
    
    s = Calloc(k/10+2, char);           /* stringified integers */
    
    PROTECT(R_str = NEW_STRING(k));
    for (j = 0; j < k; j++) {
        sprintf(s,"%i",j+1);
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



/* the end */
