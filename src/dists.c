
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
 * ceeboo 2005
 */

/* calculate Minkowsky distances (including one-norm) between the 
 * rows of two matrices. ignores NAs and NaNs (cf. dist). other than 
 * in this implementation an NA row/column gets not removed from the 
 * result matrix).
 */

SEXP pdist(SEXP R_x, SEXP R_y, SEXP R_p) {

    int nc, nx, ny;
    int i, j, k, l, m;

    double p, d, q, z;
    double *x, *y;
    
    SEXP R_obj;
    
    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
		    
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
	
    if (R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = NEW_NUMERIC(nx*ny));
		    
  
    q = 1 / p;
 
    m = 0;
    for (j = 0; j < ny; j++) {
	if (R_x == R_y)
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
    }
    
    if (R_x != R_y) {
	    
       SEXP R_tmp;

       PROTECT(R_tmp = NEW_INTEGER(2));

       INTEGER(R_tmp)[0] = nx;
       INTEGER(R_tmp)[1] = ny;
	
       SET_DIM(R_obj, R_tmp);
       UNPROTECT(1);
    }
    
    UNPROTECT(1);

    return R_obj;			    
}

/* calculate Maximum distances (cf. dist) 
 */

SEXP mdist(SEXP R_x, SEXP R_y) {

    int nc, nx, ny; 
    int i, j, k, l, m;
    
    double d, z;
    double *x, *y;
  
    SEXP R_obj;

    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("mdist: invalid number of columns");
	
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    if (R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = NEW_NUMERIC(nx*ny));
		    
    m = 0;
    for (j = 0; j < ny; j++) {
	if (R_x == R_y)
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
    }
       
    if (R_x != R_y) {
	    
       SEXP R_tmp;

       PROTECT(R_tmp = NEW_INTEGER(2));

       INTEGER(R_tmp)[0] = nx;
       INTEGER(R_tmp)[1] = ny;
	
       SET_DIM(R_obj, R_tmp);
       UNPROTECT(1);
    }
    
    UNPROTECT(1);

    return R_obj;			    
}

/* calculate Canberra distances (cf. dist)
 */

SEXP cdist(SEXP R_x, SEXP R_y) {
    
    int nc, nx, ny; 
    int i, j, k, l, m;
    
    double d, z;
    double *x, *y;

    SEXP R_obj;
  
    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("cdist: invalid number of columns");
	
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    if (R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = NEW_NUMERIC(nx*ny));
		    
    m = 0; 
    for (j = 0; j < ny; j++) {
	if (R_x == R_y)
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
    }

    if (R_x != R_y) {

       SEXP R_tmp;

       PROTECT(R_tmp = NEW_INTEGER(2));

       INTEGER(R_tmp)[0] = nx;
       INTEGER(R_tmp)[1] = ny;
	
       SET_DIM(R_obj, R_tmp);
       UNPROTECT(1);
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
    int i, j, k, l, t;
    
    int *x, *y, *s;

    double z;
 
    SEXP R_obj;
    
    nc = INTEGER(GET_DIM(R_x))[1];
    
    if (R_y == R_NilValue)
       R_y = R_x;
    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("bdist: invalid number of columns");

    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = INTEGER(R_x);
    y = INTEGER(R_y);
    
    s = Calloc(nx, int);
    
    if (R_x != R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*ny));
    else
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
   
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
	if (R_x != R_y) {
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
    }
   
    Free(s);
   
    if (R_x != R_y) {
	    
       SEXP R_tmp;

       PROTECT(R_tmp = NEW_INTEGER(2));

       INTEGER(R_tmp)[0] = nx;
       INTEGER(R_tmp)[1] = ny;

       SET_DIM(R_obj, R_tmp);
       UNPROTECT(1);
    }
    
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
    int i, j, k, l, m;
    
    double t, z;
    double *x, *y, *s;
 
    SEXP R_obj;
    
    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("ebdist: invalid number of columns");
	
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    s = Calloc(nx, double);
    
    if (R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = NEW_NUMERIC(nx*ny));
		    
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
	if (R_x == R_y) {
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
    }

    Free(s);
    
    if (R_x != R_y) {

       SEXP R_tmp;

       PROTECT(R_tmp = NEW_INTEGER(2));

       INTEGER(R_tmp)[0] = nx;
       INTEGER(R_tmp)[1] = ny;
	
       SET_DIM(R_obj, R_tmp);
       UNPROTECT(1);
    }
    
    UNPROTECT(1);

    return R_obj;
}

/* calculate angular distances (one minus cosine similarity) */

SEXP adist(SEXP R_x, SEXP R_y) {
    
    int nc, nx, ny; 
    int i, j, k, l, m;
    
    double t, z;
    double *x, *y, *s;
 
    SEXP R_obj;
	        
    nc = INTEGER(GET_DIM(R_x))[1];

    if (R_y == R_NilValue)
       R_y = R_x;
		    
    if (INTEGER(GET_DIM(R_y))[1] != nc)
       error("adist: invalid number of columns");
	
    nx = INTEGER(GET_DIM(R_x))[0];
    ny = INTEGER(GET_DIM(R_y))[0];

    x = REAL(R_x);
    y = REAL(R_y);
	
    s = Calloc(nx, double);
    
    if (R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = NEW_NUMERIC(nx*ny));
		    
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
	if (R_x == R_y) {
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
    }

    Free(s);
    
    if (R_x != R_y) {

       SEXP R_tmp;

       PROTECT(R_tmp = NEW_INTEGER(2));

       INTEGER(R_tmp)[0] = nx;
       INTEGER(R_tmp)[1] = ny;
	
       SET_DIM(R_obj, R_tmp);
       UNPROTECT(1);
    }
    
    UNPROTECT(1);

    return R_obj;
}

/* subset a dist object. in order to preserve symmetry 
 * we allow only one subset index. for the subsripting
 * to work we need a proper array (see last argument).
 * 
 * ceeboo 2005
 */

SEXP subset_dist(SEXP R_x, SEXP R_s, SEXP R_l) {

    int n, m;
    int i, ii, j, jj, k;

    int *s;

    double *x, *z;
	
    SEXP R_obj;
    
    PROTECT(R_s = arraySubscript(0, R_s, GET_DIM(R_l), getAttrib, 
						       (STRING_ELT), R_l));
	
    n = 1 + (int) sqrt(2*LENGTH(R_x));  /* number of elements */
    
    m = LENGTH(R_s);			/* size of subset */
    
    s = INTEGER(R_s);

    x = REAL(R_x);

    PROTECT(R_obj = NEW_NUMERIC(m*(m-1)/2));

    z = REAL(R_obj);
    
    k = 0;
    for (i = 0; i < m-1; i++) {
	ii = s[i]-1;
	ii = ii*(n-1)-ii*(ii+1)/2-1;
	for (j = i+1; j < m; j++) {
	    if (s[i] == s[j])
	       z[k++] = 0;		/* by definition */
	    else {
	       jj = s[j]-1;
	       z[k++] = x[jj+ii];
	    }
	}
    }
    
    /* setting the labels here is more efficient */

    R_l = GET_DIMNAMES(R_l);
    
    if (R_l != R_NilValue) {

       SEXP R_str;
       
       PROTECT(R_str = NEW_STRING(m));
       R_l = VECTOR_ELT(R_l, 0);
       for (k = 0; k < m; k++)
	   SET_ELEMENT(R_str, k, mkChar(CHAR(VECTOR_ELT(R_l, s[k]-1))));

       setAttrib(R_obj, install("Labels"), R_str);
       UNPROTECT(1);
    }
    
    UNPROTECT(2);

    return R_obj;
}

/* the end */
