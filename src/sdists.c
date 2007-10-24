
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
 * for an interface for two (atomic) sequences which returns 
 * the distance (similarity) as well as the complete edit 
 * (alignment) traceback see below.
 *
 * fixme: 1) distance between two empty sequences.
 *        2) FASTA does not seem to be GPL but we may ask
 *           if we could use parts of it in R.
 *
 * note that we do not implement the most efficient algorithmic 
 * concepts known in the field.
 *
 * (C) ceeboo 2005, 2006
 */

/* compute the operation weighted edit distance of two
 * sequences, i.e. insertion, deletion and substitution
 * may have different costs.
 *
 * note that the occurence of missing values results 
 * in NA as a match or mismatch cannot be determined.
 */

double edist_ow(int *x, int *y, double *w, int nx, int ny, int nw, 
		double *z0, char *b, double *v) {

    int i, j, x0 = 0, y0 = 0;
    double z1 = 0, z2 = 0, s0 = 0, s1 = 0, s2 = 0;

    for (i = 0; i <= nx; i++) {
	for (j = 0; j <= ny; j++)
	    if (i == 0) {
	       if (j == 0) {
		  z0[j] = z2 = 0;
		  if (b)
		     b[0] = 0;
		  if (v)
		     v[0] = 0;
	       }
	       else {
		  if (y[j-1] == NA_INTEGER)
		     return NA_REAL;
		  z2 = z0[j] = j * w[0];
		  if (b)
		     b[j*(nx+1)] = 2;
		  if (v)
		     v[j*(nx+1)] = z2;
	       }
	    } else if (j == 0) {
		      x0 = x[i-1];
		      if (x0 == NA_INTEGER)
			 return NA_REAL;
		      z1 = z2 = i * w[0];
		      if (b)
			 b[i] = 1;
		      if (v)
			 v[i] = z1;
		   }
		   else {
		      y0 = y[j-1];
		      s0 = z0[j] + w[0];
		      s1 = z1 + w[0];
		      s2 = z0[j-1] + ((x0 == y0) ? w[1] : w[2]);
		      z2 = min(s0, s1);
		      z2 = min(z2, s2);
		      if (b)
			 b[i+j*(nx+1)] = (s0 == z2) +
				         (s1 == z2) * 2 + 
				         (s2 == z2) * ((x0 != y0) ? 4 : 8);
		      if (v)
			 v[i+j*(nx+1)] = z2;
		      z0[j-1] = z1;
		      if (j == ny)
			 z0[j] = z2;
		      else
		         z1 = z2;
		   }
    }
    return z2;
}

/* compute the alphabet-weighted distance. actually, we compute 
 * the global sequential alignment with maximum similarity (see 
 * Gusfield pp. 225) and return it as a negative number.
 */

double edist_aw(int *x, int *y, double *w, int nx, int ny, int nw, 
		double *z0, char *b, double *v) {

    int i, j, x0 = 0, y0 = 0;
    double z1 = 0, z2 = 0, z3 = 0, s0 = 0, s1 = 0, s2 = 0;

    for (i = 0; i <= nx; i++) {
	for (j = 0; j <= ny; j++)
	    if (i == 0) {
	       if (j == 0) {
		  z0[j] = z2 = z3 = w[0];
		  if (b)
		     b[0] = 0;
		  if (v)
		     v[0] = z3;
	       }
	       else {
		  y0 = y[j-1];
		  if (y0 == NA_INTEGER)
		     return NA_REAL;
		  z2 = z0[j] = z0[j-1] + w[(y0-1)*nw];
		  if (b)
		     b[j*(nx+1)] = 2;
		  if (v)
		     v[j*(nx+1)] = z2;
	       }
	    } else if (j == 0) {
		      x0 = x[i-1];
		      if (x0 == NA_INTEGER)
			 return NA_REAL;
		      z3 += w[(x0-1)];
		      z1 = z2 = z3;
		      if (b)
			 b[i] = 1;
		      if (v)
			 v[i] = z1;
		   }
		   else {
		      y0 = y[j-1];
		      s0 = z0[j] + w[(x0-1)];
		      s1 = z1 + w[(y0-1)*nw];
		      s2 = z0[j-1] + w[(x0-1)+(y0-1)*nw];
		      z2 = max(s0, s1);
		      z2 = max(z2, s2);
		      if (b)
			 b[i+j*(nx+1)] = (s0 == z2) +
				         (s1 == z2) * 2 + 
				         (s2 == z2) * ((x0 != y0) ? 4 : 8);
		      if (v)
			 v[i+j*(nx+1)] = z2;
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
 *
 * notes: 1) a value of zero indicates the empty sequence.
 *	  2) an optimal non-empty solution is indicated by
 *	     the fifth bit in the traceback table
 */

double edist_awl(int *x, int *y, double *w, int nx, int ny, int nw, 
		 double *z0, char *b, double *v) {

    int i, j, x0 = 0, y0 = 0, k = 0, l = 0, *p = 0;
    double z1 = 0, z2 = 0, z = 0, s0 = 0, s1 = 0, s2 = 0;
    
    if (b)
       p = Calloc(nx*ny, int);

    for (i = 0; i <= nx; i++) {
	for (j = 0; j <= ny; j++)
	    if (i == 0) {
	       if (j == 0) {
		  z0[j] = z = 0;
		  if (b)
		     b[0] = 0;
		  if (v)
		     v[0] = 0;
		}
	       else {
		  if (y[j-1] == NA_INTEGER)
		     return NA_REAL;
		  z0[j] = 0;
		  if (b)
		     b[j*(nx+1)] = 2;
		  if (v)
		     v[j*(nx+1)] = 0;
	       }
	    } else if (j == 0) {
		      x0 = x[i-1];
		      if (x0 == NA_INTEGER)
			 return NA_REAL;
		      z1 = 0;
		      if (b)
			 b[i] = 1;
		      if (v)
			 v[i] = 0;
		   }
		   else {
		      y0 = y[j-1];
		      s0 = z0[j] + w[(x0-1)];
		      s1 = z1 + w[(y0-1)*nw];
		      s2 = z0[j-1] + w[(x0-1)+(y0-1)*nw];
		      z2 = max(0, s0);
		      z2 = max(z2, s1),
		      z2 = max(z2, s2);
		      if (b) {
			 k = i+j*(nx+1);
			 b[k] = (z2 > 0 && s0 == z2) +
				(z2 > 0 && s1 == z2) * 2 + 
				(z2 > 0 && s2 == z2) * ((x0 != y0) ? 4 : 8);
			 if (z2 > z) {
			    l = 0;
			    p[l++] = k;
			 } else
			 if (z2 > 0 && z2 == z)
			    p[l++] = k;
		      }
		      if (v)
			 v[i+j*(nx+1)] = z2;
		      if (z2 > z)
			 z = z2;
		      z0[j-1] = z1;
		      if (j == ny)
			 z0[j] = z2;
		      else
		         z1 = z2;
		   }
    }
    if (b) {
       while (l-- > 0)
	  b[p[l]] = b[p[l]] + 16;
       Free(p);
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
 * internal functions for distance computation have six + three 
 * arguments: the first two are pointers to arrays of integer (the 
 * sequences). the third is a pointer to an array of double (the 
 * weights). the seventh is a pointer to an array of double long 
 * enough to hold temporary results. the eighth is either null or
 * a pointer to a temporary array of character used in the computation
 * of the edit transcripts or alignments (see below). the ninth is 
 * either null or a pointer to an array of double large enough to 
 * hold the values of the dynamic programming table.
 * 
 * returns either a vector in lower triangular format (see dist) 
 * or a matrix of distances.
 *
 * in the case of auto-distances we check for asymmetric weights
 * as these may result in asymmetric distances.
 *
 * todo: warning if NA results are encountered (?)
 */

// test for exact symmetry

int is_symmetric(double *x, int n) {
    int i, j, r = 1;			// true
    for (i = 0; i < n-1; i++)
	for (j = i+1; j < n; j++)
	    if (x[i+j*n] != x[j+i*n]) {
		r = 0;
		break;
	    }
    return r;
}

SEXP sdists(SEXP R_x, SEXP R_y, SEXP R_method, SEXP R_weight) {
    if (TYPEOF(R_x) != VECSXP || (!isNull(R_y) && TYPEOF(R_y) != VECSXP))
	error("invalid sequence parameters");
    if (TYPEOF(R_method) != INTSXP)
	error("invalid method parameter");
    if (TYPEOF(R_weight) != REALSXP)
	error("invalid weight parameter");

    double (*sdfun)(int *, int *, double *, int, int, int, double *, char *, double *v) = NULL;
	
    int nx, ny, nw;
    int i, j, k;
    int as_matrix = 0;		/* fixes copy on write */
    SEXP x, y, t, r;		/* return value */
	
    nw = LENGTH(R_weight);
    
    switch (INTEGER(R_method)[0]) {
	case 1:
	    sdfun = edist_ow;
	    break;
	case 2:
	    if (!isMatrix(R_weight))
		error("invalid weight parameter");
	    sdfun = edist_aw;
	    nw    = INTEGER(GET_DIM(R_weight))[0];
	    break;
	case 3:
	    if (!isMatrix(R_weight))
		error("invalid weight parameter");
	    sdfun = edist_awl;
	    nw    = INTEGER(GET_DIM(R_weight))[0];
	    break;
	default:
	    error("method not implemented");
    }
    
    if (R_y == R_NilValue) {
       if (isMatrix(R_weight) && !is_symmetric(REAL(R_weight), nw))
	  error("auto-similarities for asymmetric weights not implemented");
       R_y = R_x;
    } 
    else
       as_matrix = 1;
	
    nx = LENGTH(R_x);
    ny = LENGTH(R_y);

    if (!as_matrix && R_x == R_y)
       PROTECT(r = allocVector(REALSXP, nx*(nx-1)/2));
    else
       PROTECT(r = allocMatrix(REALSXP, nx, ny));

    PROTECT(t = allocVector(REALSXP, 256));	/* temporary storage */

    k = 0;
    for (j = 0; j < ny; j++) {
	if (!as_matrix && R_x == R_y)
	   i = j + 1;
	else
	   i = 0;
	y  = VECTOR_ELT(R_y, j);
	if (LENGTH(y) >= LENGTH(t)) {		/* more storage */
	    UNPROTECT(1);
	    PROTECT(t = allocVector(REALSXP, LENGTH(y) * 2));
	}
	for (i = i; i < nx; i++) {
	    x  = VECTOR_ELT(R_x, i);
	    REAL(r)[k++] = (*sdfun)(INTEGER(x), INTEGER(y), REAL(R_weight),
				    LENGTH(x), LENGTH(y), nw, REAL(t), 0, 0);
	    R_CheckUserInterrupt();
	}
    }

    UNPROTECT(2);
	    
    return r;
}

/*  get the next edit transcript. the input arguments are a pointer 
 *  to an array of traceback codes, the lenghts of the sequences
 *  compared, and pointers to the transcript and its length. the
 *  possible edit oparations are indicated by four bits. the lowest
 *  bit is decoded and if more than one bit is set it is set to zero.
 *  
 *  returns -1 on error, 0 if no more transcripts are left, and 
 *  otherwise the backtrack position in the code array.
 */

static int next_transcript(char *b, int i, int j, char *s, int *l) {
    int b0 = 0, k0 = 0, k1 = 0, k = 0, n = i+1;
    
    while (i > 0 || j > 0) {
	if (i < 0 || j < 0) {
	    REprintf("next_transcript: coding error\n");
	    return -1;
	}
	k0 = i+j*n;
	b0 = b[k0];
	if (b0 & 1) {
	   s[k++] = 'D';
	   if (b0 & 2 || b0 & 4 || b0 & 8)
	      k1 = k0;
	   i--;
	} else
	if (b0 & 2) {
	   s[k++] = 'I';
	   if (b0 & 4 || b0 & 8)
	      k1 = k0;
	   j--;
	} else {
	   if (b0 == 4)
	      s[k++] = 'R';
	   else
	   if (b0 == 8)
	      s[k++] = 'M';
	   else {
	      REprintf("edit_transcript: coding error\n");
	      return -1;
	   }
	   i--;
	   j--;
	}
    }
    *l = k;
    s[k] = (char)0;
    
    if (k1) {
       b0 = b[k1];
       if (b0 & 1)
	  b[k1] = b0 ^ 1;
       else
       if (b0 & 2)
	  b[k1] = b0 ^ 2;
    }
    return k1;
}

/* get the next transcript for a local alignment. first we have to find 
 * the endpoint of a local alignment (if any). then we proceed as above
 * until we hit zero or either of the two sequences is exhausted. 
 * remaining prefixes or suffixes are aligned by padding with wildcards 
 * where insertions or deletions at the ends are used to acount for 
 * differences in lenghths (shifting the shorter prefix or suffix in the 
 * direction of the local alignment). we use the codes {'i', 'd', '?'} 
 * in order to distinguish these edit operations from those necessary
 * to obtain the local alignment.
 * 
 * endpoints within a local alignment are ignored as we seek only
 * alignments of maximum length. bits 6 and 7 are used as temporary
 * storage for bits 1 and 2 which we restore after all solutions for
 * one endpoint have been generated.
 */

static int next_local_transcript(char *b, int i, int j, char *s, int *l) {
    int b0 = 0, k0 = 0, k1 = 0, k2 = 0, k = 0, n = i+1, m = j+1;
   
    for (; i > 0; i--)
	for (j = m-1; j > 0; j--) {
	    k0 = i+j*n;
	    if (b[k0] & 16) {
	       k2 = k0;
	       goto next;
	    }
	}
    return 0;
next:
    
    while (k < n-i-m+j)
	s[k++] = 'd';
    while (k < m-j-n+i)
	s[k++] = 'i';
    while (k < n-i-1 || k < m-j-1)
	s[k++] = '?';
   
    while (i > 0 && j > 0) {
	k0 = i+j*n;
	b0 = b[k0];
	if (b0 == 0) 
	   break;
	else
	if (b0 & 16)
	    b[k0] = b0 = b0 ^ 16;
	if (b0 & 1) {
	   s[k++] = 'D';
	   if (b0 & 2 || b0 & 4 || b0 & 8)
	      k1 = k0;
	   i--;
	} else
	if (b0 & 2) {
	   s[k++] = 'I';
	   if (b0 & 4 || b0 & 8)
	      k1 = k0;
	   j--;
	} else {
	   if (b0 & 4)
	      s[k++] = 'R';
	   else
	   if (b0 & 8)
	      s[k++] = 'M';
	   else {
	      REprintf("edit_transcript: coding error\n");
	      return -1;
	   }
	   i--;
	   j--;
	}
    }

    for (; i > 0 && j > 0; i--, j--)
	s[k++] = '?';
    for (; i > 0; i--)
	s[k++] = 'd';
    for (; j > 0; j--)
	s[k++] = 'i';
    *l = k;
    s[k] = (char)0;
    
    if (k1) {
       b0 = b[k1];
       if (b0 & 1)
	  b[k1] = (b0 ^ 1) | 32;
       else
       if (b0 & 2)
	  b[k1] = (b0 ^ 2) | 64;
       b[k2] |= 16;
    } else
       for (k = 1; k < k2; k++) {
	    b0 = b[k];
	    if (b0 & 16)
		k1 = k;
	    if (b0 & 32)
		b0 = (b0 ^ 32) | 1;
	    if (b0 & 64)
		b0 = (b0 ^ 64) | 2;
	    b[k] = b0;
	}
    return k1;
}

/* compute the distance for two sequences and the corresponding set 
 * of equivalent edit transcripts. the input arguments are the same
 * as above with the exception that the first two are integer vectors 
 * instead of lists.
 *
 * the transcripts are coded as strings over the alphabet {'I', 'D',
 * 'R', 'M'} indicating an insert, delete, replace, or match operation
 * to be applied to the first (second) sequence (supplied). for the
 * extended symbol set for local alignments see above.
 *
 * the distance is returned as attribute 'value'. the values of the 
 * dynamic programming table may be returned as attribute 'table' for 
 * plotting, etc. For attribute 'pointer' contains an R ''segments''
 * compatible representation of the (back)pointers (see also below).
 */

SEXP sdists_transcript(SEXP R_x, SEXP R_y, SEXP R_method, SEXP R_weight, SEXP R_table) {
    if (TYPEOF(R_x) != INTSXP || TYPEOF(R_y) != INTSXP)
	error("invalid sequence parameters");
    if (TYPEOF(R_method) != INTSXP)
	error("invalid method parameter");
    if (TYPEOF(R_weight) != REALSXP)
	error("invalid weight parameter");
    if (TYPEOF(R_table) != LGLSXP)
	error("invalid option parameter");

    double (*sdfun)(int *, int *, double *, int, int, int, double *, char *, double *) = NULL;

    int (*stfun)(char *, int, int, char *, int *) = NULL;
    
    int i, j, k, n, nx, ny, nw;
    
    double d, *v = 0, *t;	// temporary storage
    char c, *b, *s;		// temporary storage
    
    SEXP r, tv = (SEXP)0, tb = (SEXP)0;

    nw = length(R_weight);
    
    switch (INTEGER(R_method)[0]) {
	case 1:
	    sdfun = edist_ow;
	    stfun = next_transcript;
	    break;
	case 2:
	    if (!isMatrix(R_weight))
		error("invalid weight parameter");
	    sdfun = edist_aw;
	    stfun = next_transcript;
	    nw    = INTEGER(GET_DIM(R_weight))[0];
	    break;
	case 3:
	    if (!isMatrix(R_weight))
		error("invalid weight parameter");
	    sdfun = edist_awl;
	    stfun = next_local_transcript;
	    nw    = INTEGER(GET_DIM(R_weight))[0];
	    break;
	default:
	    error("method not implemented");
    }

    nx = length(R_x);
    ny = length(R_y);

    if (LOGICAL(R_table)[0] == TRUE) {
	PROTECT(tv = allocMatrix(REALSXP, nx+1, ny+1));
	PROTECT(tb = allocVector(VECSXP, 4));
	v = REAL(tv);
    }
    // R.2.6.x
    b = (char *) CHAR(PROTECT(allocVector(CHARSXP, (nx+1)*(ny+1))));

    t = Calloc(ny+1, double);

    d = (*sdfun)(INTEGER(R_x),INTEGER(R_y), REAL(R_weight), nx, ny, nw, t, b, v);
    Free(t);

    if (!R_FINITE(d)) {
	UNPROTECT(1);
	if (LOGICAL(R_table)[0] == TRUE)
	    UNPROTECT(2);
	return ScalarReal(d);
    }

#ifdef TB_DEBUG
    Rprintf("traceback codes: 1 = up, 2 = left, 4 = replace, 8 = match\n\n");
    for (i = 0; i <= nx; i++) {
	Rprintf("[%2i]", i);
	for (j = 0; j <= ny; j++)
	    if (b[i+j*(nx+1)] & 16)
		Rprintf("(%2i)",b[i+j*(nx+1)] ^ 16);
	    else
		Rprintf(" %2i ",b[i+j*(nx+1)]);
	Rprintf("\n");
    }
#endif

    if (LOGICAL(R_table)[0] == TRUE) {
	int b0;
	SEXP x0, y0, x1, y1;
	
	k = 0;
	for (i = 1; i < (nx+1)*(ny+1); i++) {
	    b0 = b[i];
	    k += ((b0 & 1) == 1) + ((b0 & 2) == 2) + (((b0 & 4) == 4) || ((b0 & 8) == 8));
	}
	
	SET_VECTOR_ELT(tb, 0, (x0 = allocVector(INTSXP, k)));
	SET_VECTOR_ELT(tb, 1, (y0 = allocVector(INTSXP, k)));
	SET_VECTOR_ELT(tb, 2, (x1 = allocVector(INTSXP, k)));
	SET_VECTOR_ELT(tb, 3, (y1 = allocVector(INTSXP, k)));

	k = 0;
	for (i = 0; i <= nx; i++)
	    for (j = 0; j <= ny; j++) {
		b0 = b[i+j*(nx+1)];
		if (b0 & 1) {
		    INTEGER(x0)[k] = i-1;
		    INTEGER(y0)[k] = j;
		    INTEGER(x1)[k] = i;
		    INTEGER(y1)[k] = j;
		    k++;   
		}
		if (b0 & 2) {
		    INTEGER(x0)[k] = i;
		    INTEGER(y0)[k] = j-1;
		    INTEGER(x1)[k] = i;
		    INTEGER(y1)[k] = j;
		    k++;
		}
		if (b0 & 4 || b0 & 8) {
		    INTEGER(x0)[k] = i-1;
		    INTEGER(y0)[k] = j-1;
		    INTEGER(x1)[k] = i;
		    INTEGER(y1)[k] = j;
		    k++;
		}   
	    }
    }
    // R.2.6.x
    s = (char *) CHAR(PROTECT(allocVector(CHARSXP, nx+ny+1)));

    r = R_NilValue;
    do {
	n = (*stfun)(b, nx, ny, s, &k);
	for (i = 0; i < k/2; i++) {
	    c = s[i];
	    s[i] = s[k-i-1];
	    s[k-i-1] = c;
	}
	PROTECT(r);
	r = CONS(mkChar(s), r);
	UNPROTECT(1);
	R_CheckUserInterrupt();
    } while (n);

    UNPROTECT(2);

    PROTECT(r = PairToVectorList(r));
    SET_TYPEOF(r, STRSXP);
    setAttrib(r, install("value"), ScalarReal(d));

    if (LOGICAL(R_table)[0] == TRUE) {
	setAttrib(r, install("table"), tv);
	setAttrib(r, install("pointer"), tb);
	UNPROTECT(3);
    } else
	UNPROTECT(1);

    return r;
}

// align two sequences according to an edit transcript

SEXP sdists_align(SEXP R_x, SEXP R_y, SEXP t) {
    if (TYPEOF(R_x) != INTSXP || TYPEOF(R_y) != INTSXP)
	error("invalid sequence parameter(s)");
    if (TYPEOF(t) != STRSXP || LENGTH(t) != 1)
	error("invalid transcript parameter");
    
    int i, j, k, i0, j0;

    SEXP r, x = (SEXP)0, y = (SEXP)0;
    
    t = STRING_ELT(t, 0);

    PROTECT(r = allocVector(VECSXP, 2));

    SET_VECTOR_ELT(r, 0, (x = allocVector(INTSXP, LENGTH(t))));
    SET_VECTOR_ELT(r, 1, (y = allocVector(INTSXP, LENGTH(t))));
    
    if (isFactor(R_x)) {
	SET_LEVELS(x, duplicate(GET_LEVELS(R_x)));
	setAttrib(x, install("class"), mkString("factor"));
    }
    if (isFactor(R_y)) {
	SET_LEVELS(y, duplicate(GET_LEVELS(R_y)));
	setAttrib(y, install("class"), mkString("factor"));
    }

    UNPROTECT(1);
 
    i = j = i0 = j0 = 0;
    for (k = 0; k < LENGTH(t); k++) {
	if (i > LENGTH(R_x) || j > LENGTH(R_y))
	    error("invalid edit transcript");
	switch (CHAR(t)[k]) {
	case 'i':
	case 'I':
	    INTEGER(x)[i0++] = NA_INTEGER;
	    INTEGER(y)[j0++] = INTEGER(R_y)[j++];
	    break;
	case 'd':
	case 'D':
	    INTEGER(x)[i0++] = INTEGER(R_x)[i++];
	    INTEGER(y)[j0++] = NA_INTEGER;
	    break;
	case 'R':
	case 'M':
	case '?':
	    INTEGER(x)[i0++] = INTEGER(R_x)[i++];
	    INTEGER(y)[j0++] = INTEGER(R_y)[j++];
	    break;
	default:
	    error("invalid edit symbol");
	}
    }
    if (i < LENGTH(R_x) || j < LENGTH(R_y))
	error("invalid edit transcript");

    return r;
}

/*
 * transform a vector of transcripts into a graph, i.e. a set of edges
 * with weights the number of times the edge is a member of a path in
 * the dynamic programming table.
 *
 * returns a list of 4 vectors of coordinates for use with 'segments',
 * x0, y0, x1, y1, where x denotes the first and y the second sequence,
 * and a vector of edge frequencies.
 * 
 * notes: we code the edges into scalar integers so that we can sort 
 *	  and thus efficiently count them. the cells of the dynamic
 *	  programming table are numbered column by column. an edit 
 *	  path is therfore transformed into a sequence of indexes and 
 *	  pairs of consecutive indexes indicate entries in the edge 
 *	  table. the latter we number again by columns. the time 
 *	  complexity thus depends on sorting.
 *	  
 * fixme: Calloc may raise an error so that we cannot free memory
 *	  previously allocated with calloc or Calloc.
 *
 * ceeboo 2006
 */

SEXP sdists_graph(SEXP x) {
    if (TYPEOF(x) != STRSXP)
	error("invalid type");
    int i = 0, j = 0, h, k, l, p = 0, q = 0, k0, k1, nx = 0, ny = 0, n = 0;
    int *i0, *i1;
    SEXP r, x0, y0, x1, y1, f;

    k0 = 0;
    for (k = 0; k < LENGTH(x); k++)
	k0 += LENGTH(STRING_ELT(x, k));
    
    i0 = Calloc(k0, int);
    
    k0 = 0;
    for (h = 0; h < LENGTH(x); h++) {
	SEXP c = STRING_ELT(x, h);
	
	if (h == 0) {
	    nx = ny = LENGTH(c);
	    for (k = 0; k < LENGTH(c); k++)
		switch (CHAR(c)[k]) {
		case 'i':
		case 'I':
		    nx--;
		    break;
		case 'd':
		case 'D':
		    ny--;
		}
	    n = (nx+1) * (ny+1);
	}
	
	p = q = LENGTH(c);
	i = l = 0;
	for (k = 0; k < LENGTH(c); k++) {
	    switch (CHAR(c)[k]) {
	    case 'i':
	    case 'I':
		i += nx+1;
		p--;
		break;
	    case 'd':
	    case 'D':
		i += 1;
		q--;
		break;
	    case 'R':
	    case 'M':
	    case '?':
		i += nx+2;
		break;
	    default:
		Free(i0);
		error("invalid symbol");
	    }
	    i0[k0++] =  l + i * n;
	    l = i;
	}
	if (p != nx || q != ny) {
	    Free(i0);
	    error("transcripts do not conform");
	}
    }
    
    R_isort(i0, k0);

    i1 = Calloc(k0, int);

    l = i0[0];
    k1 = 0;
    for (k = 0; k < k0; k++) {
	if (i0[k] != l) {
	    l = i0[k];
	    i0[++k1] = l;
	}
	i1[k1]++;
    }
    k1++;

    PROTECT(r = allocVector(VECSXP, 5));
    
    SET_VECTOR_ELT(r, 0, (x0 = allocVector(INTSXP, k1)));
    SET_VECTOR_ELT(r, 1, (y0 = allocVector(INTSXP, k1)));
    SET_VECTOR_ELT(r, 2, (x1 = allocVector(INTSXP, k1)));
    SET_VECTOR_ELT(r, 3, (y1 = allocVector(INTSXP, k1)));
    SET_VECTOR_ELT(r, 4, (f  = allocVector(INTSXP, k1)));
    
    for (k = 0; k < k1; k++) {
	l = i0[k];
	i = l % n;
	j = (l - i) / n;
	INTEGER(x0)[k] = l = i % (nx+1);
	INTEGER(y0)[k] = (i - l) / (nx+1); 
	INTEGER(x1)[k] = l = j % (nx+1);
	INTEGER(y1)[k] = (j - l) / (nx+1);
	INTEGER(f )[k] = i1[k];
    }
    Free(i0);
    Free(i1);

    UNPROTECT(1);
    return r;
}

/**/
