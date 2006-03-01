
#include <R.h>
#include <Rdefines.h>

#define min(A,B) ((A)>(B) ? (B) : (A)) 
#define max(A,B) ((A)>(B) ? (A) : (B)) 

/* compute auto-distances, or cross-distances on lists of 
 * sequences, i.e. vectors of integers representing the 
 * alphabet used.
 * 
 * D. Gusfield (1997) Algorithms on Strings, Trees, and
 * Sequences. Cambridge University Press.
 * 
 * todo: interface for two (atomic) sequences which returns 
 *	 the distance (similarity) as well as the complete 
 *	 edit (alignment) traceback.
 *
 * fixme: distance between two empty sequences.
 *
 * (C) ceeboo 2005
 */

/* compute the operation weighted edit distance of two
 * sequences, i.e. insertion, deletion and substitution
 * may have different costs.
 *
 * note that the occurence of missing values results 
 * in NA as a match or mismatch cannot be determined.
 */

double edist_ow(int *x, int *y, double *w, int nx, int ny, int nw, 
		double *z0) {

    int i, j, x0, y0;
    double z1, z2, z;	

    for (i = 0; i <= nx; i++) {
	for (j = 0; j <= ny; j++)
	    if (i == 0) {
	       if (j == 0)
		  z0[j] = z2 = 0;
	       else {
		  if (y[j-1] == NA_INTEGER)
		     return NA_REAL;
		  z0[j] = j;
	       }
	    } else if (j == 0) {
		      x0 = x[i-1];
		      if (x0 == NA_INTEGER)
			 return NA_REAL;
		      z1 = z2 = i;
		   }
		   else {
		      y0 = y[j-1];
		      if (x0 == y0)
			 z = w[1];
		      else
			 z = w[2];
		      z2 = min(z0[j],z1)+w[0];
		      z2 = min(z2,z0[j-1]+z);
		      z0[j-1] = z1;
		      if (j == ny)
			 z0[j] = z2;
		      else
		         z1 = z2;
		   }
    }
    return z2;
}

/* compute the alphabet-weighted distance. actually, 
 * we compute the sequential alignment with maximum 
 * similarity (see Gusfield pp. 225) and return it as 
 * a negative number.
 */

double edist_aw(int *x, int *y, double *w, int nx, int ny, int nw, 
		double *z0) {

    int i, j, x0, y0;
    double z1, z2;	

    for (i = 0; i <= nx; i++) {
	for (j = 0; j <= ny; j++)
	    if (i == 0) {
	       if (j == 0)
		  z0[j] = z2 = w[0];
	       else {
		  y0 = y[j-1];
		  if (y0 == NA_INTEGER)
		     return NA_REAL;
		  z0[j] = w[y0*nw];
	       }
	    } else if (j == 0) {
		      x0 = x[i-1];
		      if (x0 == NA_INTEGER)
			 return NA_REAL;
		      z1 = z2 = w[x0];
		   }
		   else {
		      y0 = y[j-1];
		      z2 = max(z0[j]+w[x0],z1+w[y0*nw]),
		      z2 = max(z2,z0[j-1]+w[x0+y0*nw]);
		      z0[j-1] = z1;
		      if (j == ny)
			 z0[j] = z2;
		      else
		         z1 = z2;
		   }
    }
    return -z2;
}

/* as above but align locally instead of globally.
 */

double edist_awl(int *x, int *y, double *w, int nx, int ny, int nw, 
		 double *z0) {

    int i, j, x0, y0;
    double z1, z2, z;	

    for (i = 0; i <= nx; i++) {
	for (j = 0; j <= ny; j++)
	    if (i == 0) {
	       if (j == 0)
		  z0[j] = z = 0;
	       else {
		  if (y[j-1] == NA_INTEGER)
		     return NA_REAL;
		  z0[j] = 0;
	       }
	    } else if (j == 0) {
		      x0 = x[i-1];
		      if (x0 == NA_INTEGER)
			 return NA_REAL;
		      z1 = 0;
		   }
		   else {
		      y0 = y[j-1];
		      z2 = max(0,z0[j]+w[x0]);
		      z2 = max(z2,z1+w[y0*nw]),
		      z2 = max(z2,z0[j-1]+w[x0+y0*nw]);
		      if (z2 > z)
			 z = z2;
		      z0[j-1] = z1;
		      if (j == ny)
			 z0[j] = z2;
		      else
		         z1 = z2;
		   }
    }
    return -z;
}
/* provide a common interface to all internal functions that
 * compute distances (similarities) of sequences.
 * 
 * we expect two lists of integer vectors, an integer code for
 * the internal function, and a double vector (matrix) for the 
 * weights to use.
 * 
 * internal functions for distance computation have six + one 
 * arguments: the first two are pointers to arrays of integer (the 
 * sequences). the third is a pointer to an array of double (the 
 * weights). the last is a pointer to an array of double long enough 
 * to hold temporary results.
 * 
 * returns either a vector in lower triangular format (see dist) 
 * or a matrix of distances.
 */

SEXP sdists(SEXP R_x, SEXP R_y, SEXP R_method, SEXP R_weight) {

    double (*sdfun)(int *, int *, double *, int, int, int, double *) = NULL;
	
    int nx, ny, nw;
    int i, j, k, l, lx, ly;
    int nt = 32;		/* initial length of temporary storage */
    
    double *t;
    
    SEXP x, y;			/* pointer to sequence vector */

    SEXP R_obj;			/* return value */
	
    nw = length(R_weight);
    
    switch (INTEGER(R_method)[0]) {
	case 1:
	    sdfun = edist_ow;
	    break;
	case 2:
	    sdfun = edist_aw;
	    nw    = INTEGER(GET_DIM(R_weight))[0];
	    break;
	case 3:
	    sdfun = edist_awl;
	    nw    = INTEGER(GET_DIM(R_weight))[0];
	    break;
	default:
	    error("method not implemented");
    }
    
    if (R_y == R_NilValue)
       R_y = R_x;
	
    nx = length(R_x);
    ny = length(R_y);

    if (R_x == R_y)
       PROTECT(R_obj = NEW_NUMERIC(nx*(nx-1)/2));
    else
       PROTECT(R_obj = NEW_NUMERIC(nx*ny));

    t  = Calloc(nt, double);			/* temporary storage */
    
    k = 0;
    for (j = 0; j < ny; j++) {
	if (R_x == R_y)
	   i = j + 1;
	else
	   i = 0;
	y  = VECTOR_ELT(R_y,j);
	ly = length(y);
	for (i = i; i < nx; i++) {
	    x  = VECTOR_ELT(R_x,i);
	    lx = length(x);
	    l  = min(lx,ly);
	    if (l >= nt) {
	       do     nt *= 2; 
	       while (l >= nt);
	       t  = Realloc(t,nt,double);
	    }
	    if (lx < ly)
	       REAL(R_obj)[k++] = (*sdfun)(INTEGER(y),INTEGER(x),
			                   REAL(R_weight), ly,lx,nw,t);
	    else
	       REAL(R_obj)[k++] = (*sdfun)(INTEGER(x),INTEGER(y),
					   REAL(R_weight), lx,ly,nw,t);
	}					
    }
    Free(t);

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

/**/
