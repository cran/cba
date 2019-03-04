/*
 * proximus.c
 *
 * implements the paper:
 * 
 * M. Koyutürk, A. Graham, and N. Ramakrishnan. Compression, Clustering, 
 * and Pattern Discovery in Very High-Dimensional Descrete-Attribute Data 
 * Sets. IEEE Transactions On Knowledge and Data Engineering, Vol. 17, 
 * No. 4, (April) 2005
 * 
 * As data mining algorithms are supposed to deal with large amounts of 
 * data this code was optimized for low memory footprint and low execution 
 * times. Fasten your seat belts ;-)
 *
 * 
 * (C) ceeboo 2005
*/

#include <R.h>
#include <Rdefines.h>

static int debug = FALSE;  /* user defined */

/* copy a variable to R */

static SEXP var2R(int v) {

    SEXP R_obj;

    R_obj = NEW_INTEGER(1);

    INTEGER(R_obj)[0] = v;

    return R_obj;
}

/* vector for counting and indexing */

typedef struct { 
    int *v;	/* pointer to array of values */
    int n;	/* number elements */
} VEC;

static VEC *newVec(int n) {

    int *v = Calloc(n, int);
    VEC *p = Calloc(1, VEC);
   
    p->v = v;
    p->n = n;
    
    return p;
}

static VEC *copyVec(VEC *v) {
	
    int i;
    VEC *r = newVec(v->n);
    
    for (i = 0; i < v->n; i++)
	r->v[i] = v->v[i];

    return r;
}

static void freeVec(VEC *v) {

    if (v->v != NULL)
       Free(v->v);
    Free(v);
}

/* copy a vector to R where an offset is added 
 * to each element */

static SEXP vec2R(VEC *v, int o) {

    int j;
    
    SEXP R_obj;

    R_obj = NEW_INTEGER(v->n);
 
    for (j = 0; j < v->n; j++)
   	INTEGER(R_obj)[j] = v->v[j] + o;
    
    return R_obj;
}

/* for debugging

static void VecPrintf(VEC *v, char *s) {

    int j;
    
    Rprintf("%s", s);
    for (j = 0; j < v->n; j++) 
	Rprintf(" (%i,%i)", j, v->v[j]);
    Rprintf("\n");
}

*/

/* matrix for binary data in sparse column format */

typedef struct {
    int *ri;	    /* pointer to array of row indexes */
    int *ci;	    /* pointer to array of column start indexes */
    int nr;	    /* number of rows */
    int nc;	    /* number of columns */
} MAT;

static void freeMat(MAT *m) {
	
    Free(m->ri);
    Free(m->ci);
    Free(m);
}

/* copy and transpose R matrix to sparse matrix */

static MAT *R_mat2mat(SEXP R_mat) {

    extern int debug;

    int nr, nc;
    int i, j, k, n;
    int *x, *ci, *ri;
    
    MAT *m;
    
    x = INTEGER(R_mat);
    
    nr = INTEGER(GET_DIM(R_mat))[0];		/* number of rows */
    nc = INTEGER(GET_DIM(R_mat))[1];		/* number of columns */

    ci = Calloc(nr+1, int);			/* column start */

    n = 1024;					/* initial size */
    ri = Calloc(n, int);			/* row indexes */

    k = 0;
    for (j = 0; j < nr; j++) {
	ci[j] = k;
	for (i = 0; i < nc; i++) 
	    if (x[i * nr + j] == 1) {
	       if (k == n) {			/* used up */
		  n *= 2;			/* double size */
		  ri = Realloc(ri, n, int);
	       }
	       ri[k++] = i;
	    }
    }
    ci[j] = k;					/* length of ri */

    if (n > k)					/* free unused */
       ri = Realloc(ri, k, int);

    
    if (debug) {
       Rprintf("Non-Zero: %i\n", k);
       Rprintf("Sparsity: %4.2f\n",k / (double) (nr * nc));
    }
    
    m = Calloc(1, MAT);

    m->ri = ri;
    m->ci = ci;
    m->nr = nc;
    m->nc = nr;

    return m;
}

/* multiply a matrix in sparse column format (m) with a sparse 
 * vector (v) from the left using a subset of the columns (s). the 
 * caller is reponsible for providing a proper results vector (r). */

static void matLeft(VEC *r, VEC *v, VEC *s, MAT *m) {

    int i, j, k, z;
    
    for (i = 0; i < s->n; i++) {	/* columns */
	z = 0;
	k = 0;
	j = m->ci[s->v[i]];
	do {				/* rows */
	   if (m->ri[j] == v->v[k]) {
	      z++;
	      j++;
	      k++;
	   }
	   else if (m->ri[j] < v->v[k])
		   j++; 
		else
		   k++;
	} while (j < m->ci[s->v[i] + 1] && k < v->n);
	r->v[i] = z;
    }
    r->n = s->n;
}

/* as above but multiply from the right */ 

static void matRight(VEC *r, VEC *v, MAT *m) {

    int i, j;

    for (i = 0; i < m->nr; i++)
	r->v[i] = 0;
    r->n = m->nr;
    
    for (i = 0; i < v->n; i++)
	for (j = m->ci[v->v[i]]; j < m->ci[v->v[i] + 1]; j++)
	    r->v[m->ri[j]]++;
}

/* linked list of approximation results */

typedef struct resNode {
    VEC *x;			/* presence vector (column indexes) */
    VEC *y;			/* dominant pattern vector (row indexes) */
    int n;			/* number of ones ... */
    int c;			/* approximation criterion */
    int r;			/* hamming radius */
    struct resNode *next;	/* pointer to result element */
} RES;

static int  res_cnt;		/* number of result elements */
static RES *res_last;		/* last element of result list */

static int freeRes(RES *r) {
	
    int i;
    RES *p, *q;
    
    i = 0;
    for (p = r; p != NULL; p = q) {
	q = p->next;
	
	freeVec(p->x);
	freeVec(p->y);
	
	Free(p);
	i++;
    }

    return i;
}

/* copy result list to R and clean up
 *
 * fixme: pointer protection should be 
 *	  on the level of the caller
 */

static SEXP res2R(RES *r, MAT *m) {

    int i, nr, nc;
    RES *p, *q;

    SEXP R_ret, R_obj, R_lst, R_res;
    
    nc = m->nr;				    /* transpose */
    nr = m->nc;
    
    PROTECT(R_ret = NEW_LIST(3));	    /* results header */
    
    SET_ELEMENT(R_ret, 0, PROTECT(var2R(nr)));
    SET_ELEMENT(R_ret, 1, PROTECT(var2R(nc)));
    UNPROTECT(2);

    PROTECT(R_obj = NEW_STRING(3));
    
    SET_STRING_ELT(R_obj, 0, mkChar("nr"));
    SET_STRING_ELT(R_obj, 1, mkChar("nc"));
    SET_STRING_ELT(R_obj, 2, mkChar("a"));

    SET_NAMES(R_ret, R_obj);

    UNPROTECT(2);
    
    PROTECT(R_lst = NEW_LIST(res_cnt));	    /* results list */
    
    i = 0;
    for (p = r; p != NULL; p = q) {
	q = p->next;
		
	PROTECT(R_res = NEW_LIST(5));
	
	SET_ELEMENT(R_res, 0, PROTECT(vec2R(p->x,1)));
	SET_ELEMENT(R_res, 1, PROTECT(vec2R(p->y,1)));
	UNPROTECT(2);
	
	SET_ELEMENT(R_res, 2, PROTECT(var2R(p->n)));
	SET_ELEMENT(R_res, 3, PROTECT(var2R(p->c)));
	SET_ELEMENT(R_res, 4, PROTECT(var2R(p->r)));
	UNPROTECT(3);

	freeVec(p->x);
	freeVec(p->y);
	Free(p);
	
        PROTECT(R_obj = NEW_STRING(5));
    
        SET_STRING_ELT(R_obj, 0, mkChar("x"));
        SET_STRING_ELT(R_obj, 1, mkChar("y"));
        SET_STRING_ELT(R_obj, 2, mkChar("n"));
        SET_STRING_ELT(R_obj, 3, mkChar("c"));
        SET_STRING_ELT(R_obj, 4, mkChar("r"));

        SET_NAMES(R_res, R_obj);
	UNPROTECT(1);
	
	if (i == res_cnt) {
	    i += freeRes(q);
		 freeMat(m);
	    error("res2R result count error [%i:%i]", i, res_cnt);
	}
	SET_ELEMENT(R_lst, i++, R_res);
	UNPROTECT(1);
    }
    if (i != res_cnt)
       error("res2R result count error [%i:%i]", i, res_cnt);
    
    SET_ELEMENT(R_ret, 2, R_lst);
    
    UNPROTECT(1);
    
    return R_ret;
}

/* compute the rank-one approximation of a column subset of a 
 * matrix. the code is optimized for minimal memory usage */

static int min_size = 1;	/* user defined */
static int max_iter = 16;	/* user defined */

static RES *approximate(VEC *s, MAT *m) {

    extern int min_size;		/* minimum set size */
    extern int max_iter;		/* maximum iterations */
    extern int debug;
    
    int i, j, l, c, z;

    VEC *x, *y, *v;
    RES *p;
    
    x = newVec(s->n);			/* presence set (column indexes) */
    y = newVec(m->nr);			/* dominant pattern (row indexes) */
    
    v = newVec((s->n > m->nr) ? 
		s->n : m->nr);		/* result vector (counts) */

    if (s->n > min_size) {
       i = (int) (unif_rand() * s->n);	/* sample a column */
       y->n = 0;
       for (j = m->ci[s->v[i]]; j < m->ci[s->v[i] + 1]; j++)
	   y->v[y->n++] = m->ri[j];
    } else {
       for (j = 0; j < s->n; j++)
	   x->v[j] = s->v[j];
    }
    z =  0;				/* number of ones in pattern */
    c = -1;				/* stopping criterion */
    i =  0;
    while (i < max_iter) {
	if (s->n > min_size) {
	   matLeft(v, y, s, m);
	   x->n = 0;
	   for (j = 0; j < v->n; j++)
	       if (2 * v->v[j] >= y->n)	/* holds for at least one */
	          x->v[x->n++] = s->v[j];
	}
	
	matRight(v, x, m);
	z = 0;
	y->n = 0;
	for (j = 0; j < v->n; j++)
	    if (2 * v->v[j] >= x->n) {	/* may not hold for any */
	       z += v->v[j];
	       y->v[y->n++] = j;
	    }

	l = c;
	c = 2 * z - x->n * y->n;
	if (c == l)			/* convergence */
	   break;
	i++;
	if (debug > 1)
	   Rprintf("%2i %6i %i\n", i, x->n, c);	
    }
    if (i == max_iter)			/* no convergence */
       warning("approximation: no convergence");

    /* compute the Hamming radius of the presence set */ 
    
    matLeft(v, y, x, m);
    l = 0;
    for (i = 0; i < x->n; i++) {
	j = m->ci[x->v[i] + 1] - m->ci[x->v[i]];
	j += y->n - 2 * v->v[i];
	if (j > l)
	   l = j;
    }
    
    freeVec(v);
    
    x->v = Realloc(x->v, x->n, int);	/* see above */
    if (y->n)				/* see above */
       y->v = Realloc(y->v, y->n, int);
    else {
       Free(y->v);
       y->v = NULL;
    }

    /* package result */
    
    p = Calloc(1, RES);
    
    p->x = x;
    p->y = y;
    p->n = z;
    p->c = c;
    p->r = l;
    
    p->next = NULL;
   
    return p;
}

/* produce a presence set. for now, draw a pattern 
 * and select additional patterns that are within the user
 * defined radius. this may result in a singular set and has 
 * nothing todo with the approximation idea! this is more
 * like vodoo.
 */

static int max_radius = 1;	/* user defined */

static VEC *presenceSet(VEC *s, MAT *m) {

    extern int debug;
    extern int max_radius;
	
    int i, j, z;
    
    VEC *y, *x;
    
    y = newVec(m->nr);			/* pattern vector */
    x = newVec(s->n);			/* presenece vector */

    i = (int) (unif_rand() * s->n);	/* sample a column */
    y->n = 0;
    for (j = m->ci[s->v[i]]; j < m->ci[s->v[i] + 1]; j++)
	y->v[y->n++] = m->ri[j];

    /* select all rows that are within the
     * radius of the selected pattern */
    
    matLeft(x, y, s, m);
    x->n = 0;
    for (i = 0; i < s->n; i++) {
	z = m->ci[s->v[i] + 1] - m->ci[s->v[i]];
	z += y->n - 2 * x->v[i];
	if (z <= max_radius)
	   x->v[x->n++] = s->v[i];
    }
    if (debug > 1)
       Rprintf(" %i %i\n", s->n, x->n);
    
    freeVec(y);
    
    x->v = Realloc(x->v, x->n, int);
   
    return x;
}

/* remove the set x from set s (column indexes). note: this is
 * not a general implementation of the setminus operation. */

static void remSet(VEC *x, VEC *s) {

    int i, j, k;
    
    j = 0;
    k = 0;
    for (i = 0; i < s->n; i++)
        if (j < x->n && x->v[j] == s->v[i])
           j++;
        else
           s->v[k++] = s->v[i];
    s->n = k;
}

/* partition a binary matrix over the columns. the code is optimized
 * for minimal memory usage. for the sake of algorithmic clarity
 * shortcuts with respect to terminal nodes are not implemented. */

static int min_retry  = 10;	/* user defined */

static RES *partition(VEC *s, MAT *m, int d, int i) {

    extern int max_radius;
    extern int min_retry;
    extern int min_size;			    /* see approximation */
    extern int debug;
    
    extern int  res_cnt;
    extern RES *res_last;

    VEC *xx, *ss;
    RES *z, *zz;

    z = approximate(s, m);
    
    if (debug)
       Rprintf("%3i [%i,%i,%i] %i", d, s->n, z->x->n, z->r, i);
    
    if (z->x->n == s->n) {			    /* pure */
       if (z->r <= max_radius ||		    /* homogenous */
	   z->x->n <= min_size) { 		    /* min size */
	  res_cnt++;
	  
	  if (debug)
	     Rprintf(" * %i\n", res_cnt);
	  
	  return res_last = z;
       }
       else if (min_retry &&
		s->n >= min_retry * i) {	    /* retry */
	       if (debug)
		  Rprintf(" +\n");
	       
	       freeRes(z);

	       return partition(s, m, d, i+1);
	    }
	    else {				    /* vodoo !!! */
	       if (debug)
		  Rprintf(" >>\n");
	  
	       freeRes(z);
	       
	       /* compare below */
	 
	       xx = presenceSet(s, m);
	       ss = copyVec(xx);
	       zz = partition(ss, m, d+1, i);
	       
	       freeVec(ss);
	       
	       remSet(xx, s);
	       
	       freeVec(xx);

	       if (s->n) {
		  z = res_last;
		  z->next = partition(s, m, d+1, i);
	       }

	       if (debug)
	          Rprintf("%3i <<\n", d);
	       
	       return zz;
       }
    }
    if (debug)
       Rprintf(" >\n");
  
    /* in order to prevent excessive memory consumption we reuse the 
     * subset vector for the next zero set. as its contents may get 
     * changed in the recursion, the next one set must be a copy of 
     * the current one set. */
    
    ss = copyVec(z->x);
    zz = partition(ss, m, d+1, i);
   
    freeVec(ss);
    
    remSet(z->x, s);
    
    freeRes(z);

    z = res_last;    
    z->next = partition(s, m, d+1, i);
    
    if (debug)
       Rprintf("%3i <\n", d);
    
    return zz;
}

/* R interface */

SEXP proximus(SEXP R_mat, SEXP R_max_radius, SEXP R_min_size, SEXP R_min_retry,						     SEXP R_max_iter, SEXP R_debug) {

    extern int max_radius;		    /* see partition */
    extern int min_size;
    extern int min_retry;
    extern int max_iter;		    /* see approximation */
    extern int debug;

    extern int res_cnt;

    int j;
    
    VEC *s;
    MAT *m;
    RES *r;

    SEXP R_res;
   
    if (!LENGTH(R_max_radius) || 
	!LENGTH(R_min_size  ) || 
	!LENGTH(R_min_retry ) || 
	!LENGTH(R_max_iter  ) || 
	!LENGTH(R_debug     ))
	error("proximus: missing parameter");

    max_radius = INTEGER(R_max_radius)[0];
    min_size   = INTEGER(R_min_size  )[0];
    min_retry  = INTEGER(R_min_retry )[0];
    max_iter   = INTEGER(R_max_iter  )[0];
    debug      = LOGICAL(R_debug     )[0];

    if (!IS_LOGICAL(R_mat))
       error("proximus: matrix not logical");
    
    m = R_mat2mat(R_mat);
    
    s = newVec(m->nc);			    /* column subset vector */

    for (j = 0; j < s->n; j++)
	s->v[j] = j;

    GetRNGstate();

    res_cnt = 0;			    /* reset results counter */
    
    r = partition(s, m, 0, 1);		    /* recursion */

    PutRNGstate();
    
    freeVec(s);

    R_res = res2R(r, m);

    freeMat(m);
    
    UNPROTECT(1);

    return R_res;
}

/***/
