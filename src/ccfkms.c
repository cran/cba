/*
 *  ccfkms.c
 * 
 *  parameter and logistic k-means based on conjugate covex functions
 *  using sparse data structures and centering (or optionally 
 *  standardizing) of the data.
 *  
 *  For details see:
 * 
 *  Helmut Strasser and Klaus Poetzelberger. Data Compression by 
 *  Unsupervised Classification. SFB Report Series, No. 10, 1997.
 *                                                          
 *  convex  function: f(x) = |x|^q/q
 *
 *  kohonen  k-means: q = 1
 *  ordinary k-means: q = 2
 *
 *  convex  function: f(x) = 2*ln(cosh(|x|))/2
 *
 *  logistic k-means: q = 0
 *  
 *  Sparse data structure means that zero data values in the data are not 
 *  stored because they can be ignored in all vector operations, so that 
 *  computations get a considerable boost if the data are highly sparse, 
 *  i.e. the ratio of the number of non-zero data values to the number of 
 *  all data values is low.
 *
 *  Centering means that the mean of the data is subtracted from the 
 *  samples, and standardizing that we further devide by the standard 
 *  deviation. Note that only for Euclidian distances the solution does 
 *  not depend on centering. 
 *
 *  In the case of degenerate and non-convergent solutions the program 
 *  gives a warning message.
 *
 *  fixme: 1) second, ... winner not implemented.
 *         2) maybe return a flag indicating convergence issues?
 *  
 *  Note that the code is prepared for direct intefacing with sparse 
 *  data structures, such as dgCMatrix from the R package Matrix.
 *  
 *  (C) ceeboo 2003, 2004, 2005  
*/    

#include <R.h>
#include <Rdefines.h>

static int debug = FALSE;	    /* be silent */

/* matrix structure for data in sparse row format */

typedef struct {
    int *ri;        /* pointer to array of column start indexes */
    int *ci;        /* pointer to array of column indexes */
    double *cv;	    /* pointer to array of column values */
    int nr;         /* number of rows */
    int nc;         /* number of columns */
} SMAT;

static void FreeSMat(SMAT *m) {
	
    Free(m->ri);
    Free(m->ci);
    Free(m->cv);
    
    Free(m);
}

/* copy R matrix in full-storage representation to sparse 
 * representation. treat as read only (!) */

static SMAT *R_mat2smat(SEXP R_mat) {

    extern int debug;

    int nr, nc, n;
    int i, j, k;
    int *ri, *ci;

    double z;
    double *x, *cv;

    SMAT *m;
    
    nr = INTEGER(GET_DIM(R_mat))[0];
    nc = INTEGER(GET_DIM(R_mat))[1];

    x = REAL(R_mat);
  
    ri = Calloc(nr+1, int);			    /* row start indexes */
    
    n = 1024;					    /* initial memory */

    ci = Calloc(n, int);			    /* column indexes */
    cv = Calloc(n, double);			    /* column values */

    k = 0;
    for (i = 0; i < nr; i++) {			    /* rows */
	ri[i] = k;
	for (j = 0; j < nc; j++) {		    /* columns */
	    z = x[i+j*nr];
	    if (R_FINITE(z) && z != 0.0) {
	       if (k == n) {
		  n *= 2;			    /* double memory */
		  ci = Realloc(ci, n, int);
		  cv = Realloc(cv, n, double);
	       }
	       ci[k] = j;
	       cv[k++] = z;
	    }
	}
    }
    ri[i] = k;
    
    if (n > k) {
       ci = Realloc(ci, k, int);
       cv = Realloc(cv, k, double);
    }

    if (debug) {
       Rprintf("Non-Zero: %i\n", k);
       Rprintf("Sparsity: %4.2f\n",k / (double) (nr * nc));
    }
    
    m = Calloc(1, SMAT);
    
    m->ri = ri;
    m->ci = ci;
    m->cv = cv;
    m->nr = nr;
    m->nc = nc;

    return m;
}

SEXP ccfkms(SEXP R_x, SEXP R_p, SEXP R_par, SEXP R_max_iter, SEXP R_opt_std, 
							     SEXP R_debug) {
    extern int debug;
    
    int opt_std, np, max_iter;
    int i, j, k, l, iter, ap;
    int *pf, *pm;

    double par;
    double x = 0, y, z, max_var, max_inf, old_inf, inf, var;
    double *am, *as, *p, *pt, *cc, *ct;

    char *s;

    SMAT *m;

    SEXP R_obj, R_tmp;
    
    debug = INTEGER(R_debug)[0];
	
    m = R_mat2smat(R_x);			/* data matrix */
    
    /* compute attribute means. standardization 
     * is optional. if used we transform so that 
     * we do not need to revert to a full-storage
     * representation. */ 

    opt_std = INTEGER(R_opt_std)[0];		/* standardization option */
  
    am = Calloc(m->nc, double);			/* attribute means */
    as = NULL;					/* attribute standard
							     deviations */ 
    if (opt_std)
       as = Calloc(m->nc, double);
    
    for (i = 0; i < m->nr; i++)
	for (j = m->ri[i]; j < m->ri[i + 1]; j++) {
	    am[m->ci[j]] += m->cv[j];
	    if (opt_std)
	       as[m->ci[j]] += pow(m->cv[j], 2);
	}

    for (i = 0; i < m->nc; i++) {
	am[i] /= m->nr;
	if (opt_std) {
	   as[i] = sqrt(as[i] / m->nr - pow(am[i], 2));
           if (as[i] == 0) {
	      Free(am);
	      if (opt_std)
	         Free(as);
	      FreeSMat(m);
              error("ccfkms: zero standard deviation");
	   }
	   am[i] /= as[i];
	}
    }

    if (opt_std)				/* prepere data */
       for (i = 0; i < m->nr; i++)
           for (j = m->ri[i]; j < m->ri[i + 1]; j++)
               m->cv[j] /= as[m->ci[j]];

    /* get initial protoypes and allocate 
     * R result object. */

    np = INTEGER(GET_DIM(R_p))[0];		/* number of prototypes */
    
    if (INTEGER(GET_DIM(R_p))[1] != m->nc) {	/* check */
       Free(am);
       if (opt_std)
          Free(as);
       FreeSMat(m);
       error("ccfkms: \"x\" and \"p\" do not conform");
    }
    
    PROTECT(R_obj = NEW_LIST(4));			/* result object */
    
    PROTECT(R_tmp = allocMatrix(REALSXP, np, m->nc));	/* prototypes */

    Memcpy(REAL(R_tmp), REAL(R_p), np * m->nc);		/* copy prototypes */
    
    p = REAL(R_tmp);
    
    SET_VECTOR_ELT(R_obj, 0, R_tmp);
    UNPROTECT(1);

    /* center (standardize) initial prototypes */
    
    for (i = 0; i < np; i++)
	for (j = 0; j < m->nc; j++) {
	    if (opt_std)
	       p[i + j * np] /= as[j];
	    p[i + j * np] -= am[j];
	}
    
    /* get parameter */

    par = REAL(R_par)[0];
  
    /* get maximum number of iterations */

    max_iter = INTEGER(R_max_iter)[0];
  
    /* compute the maximum information and variance, 
     * i.e., each point is a prototype */

    z = 0;
    for (i = 0; i < m->nc; i++)
	z += pow(am[i], 2);

    max_var = 0;
    max_inf = 0;
    for (i = 0; i < m->nr; i++) {
	y = z;
	for (j = m->ri[i]; j < m->ri[i + 1]; j++)
	    y += m->cv[j] * (m->cv[j] - 2 * am[m->ci[j]]);
	max_var += y;
	max_inf += pow(sqrt(y), par) / par;
    }
    max_var /= m->nr;
    max_inf /= m->nr;

    /* allocate remaining result objects 
     * and iterate to a fixpoint solution  */

    PROTECT(R_tmp = NEW_INTEGER(np));	    /* prototype frequencies */
    
    pf = INTEGER(R_tmp);

    SET_VECTOR_ELT(R_obj, 1, R_tmp);
    UNPROTECT(1);
    
    PROTECT(R_tmp = NEW_INTEGER(m->nr));    /* prototype memberships */
    
    pm = INTEGER(R_tmp);
  
    SET_VECTOR_ELT(R_obj, 2, R_tmp);
    UNPROTECT(1);

    GetRNGstate();
    
    pt = Calloc(np * m->nc, double);	    /* prototype temporary */
    cc = Calloc(np, double);		    /* conjugate convex */
    ct = Calloc(np, double);		    /* conjugate convex temporary */

    if (debug)
       Rprintf("\n %3s %5s %5s %3s\n","#","inf","var","nap");

    old_inf = -1;
    inf = 0;
    iter = 0;
    while(inf > old_inf && iter < max_iter) {

	/* map prototype means into domain of dual problem 
	 * and compute conjugate convex function */

	for (i = 0; i < np; i++) {
	    y = 0;
	    for (j = 0; j < m->nc; j++) 
		y += pow(p[i + j * np], 2);
	    y = sqrt(y);
	    if (par)
	       z = pow(y, par-2);
	    else {				    /* logistic */
	       x = (exp(y) - 1) / (exp(y) + 1);
	       z = x / y;
	    }
	    for (j = 0; j < m->nc; j++) 
		p[i + j * np] *= z;
	    z = 0;
	    for (j = 0; j < m->nc; j++)
		z += p[i + j * np] * am[j];
	    if (par)
	        z += pow(y, par) * (par-1) / par;
	    else
		z += (1 + x) * log(1 + x) + (1 - x) * log(1 - x);
	    cc[i] = z;
	    
	    pf[i] = 0;				    /* initialize */
	    for (j=0; j < m->nc; j++)
		pt[i + j * np] = 0;
	}

	/* determine partition and 
	 * calculate prototype means */

	for (i = 0; i < m->nr; i++) {
	    for (k = 0; k < np; k++) {
		ct[k] =- cc[k];
		for (j = m->ri[i]; j < m->ri[i + 1]; j++)
		    ct[k] += m->cv[j] * p[k + m->ci[j] * np];
	    }

	    /* find the closest prototype.
	     * note that tie breaking is used. */ 

	    l = 1;
	    k = 0;
	    z = ct[0];
	    for (j = 1; j < np; j++)
		if (z < ct[j]) {
		   k = j;
		   z = ct[j];
		}
	        else if (z == ct[j]) {
		        if (unif_rand() > l/(l+1))
		           k = j;
		        l++;
		     }
	    pm[i] = k;

	    /* update prototype frequency and means */

	    pf[k]++;
	    for (j = m->ri[i]; j < m->ri[i + 1]; j++)
		pt[k + m->ci[j] * np] += m->cv[j];
	}

	/* update the stopping criterion. compute the means 
	 * and the information and variance of the partition */
     
	old_inf = inf;

	ap = 0;
	inf = 0;
	var = 0;
	for (i = 0; i < np; i++) {
	    if (pf[i] != 0) {
	       ap++;
	       for (j = 0; j < m->nc; j++)
		   p[i + j * np] = pt[i + j * np] / pf[i] - am[j];
	       z =0;
	       for (j = 0; j < m->nc; j++)
		   z += pow(p[i + j * np], 2);
	       var += z * pf[i] / m->nr;
	       z = sqrt(z);
	       if (par)
		  z = pow(z, par) / par;
	       else {
		  z = 2 * log((exp(z / 2) + exp(-z / 2)) / 2);
	       }
	       inf += z * pf[i] / m->nr;
	    }
	} 
	iter++;
	
	if (debug)
	   Rprintf(" %3i %5.3f %5.3f %3i\n", iter, inf / max_inf, 
					           var / max_var, ap);
	/* degenrate solution */
	if (old_inf > inf)
	   warning("ccfkms: decrease in information");
    }
    Free(pt);
    Free(cc);
    Free(ct);

    PutRNGstate();
    
    if (max_iter > 1 && old_inf != inf)
       warning("ccfkms: no convergence");
  
    /* invert the information */
    
    inf = max_inf - inf;

    PROTECT(R_tmp = NEW_NUMERIC(1));
    
    REAL(R_tmp)[0] = inf;

    SET_VECTOR_ELT(R_obj, 3, R_tmp);
    UNPROTECT(1);
    
    /* decenter (destandardize) the prototype means */

    for (i = 0; i < np; i++) {
	for (j = 0; j < m->nc; j++) {
	    p[i + j * np] += am[j];
	    if (opt_std)
	       p[i + j * np] *= as[i];
	}
    }
    Free(am);
    
    if (opt_std)
       Free(as);
    
     /* offset memberships to R indexing. 
      */
    
    for (i = 0; i < m->nr; i++)
	pm[i]++;

    FreeSMat(m);
    
    /* levels attribute */

    s = Calloc(np/10+2, char);		    /* stringified integers */
        
    PROTECT(R_tmp = NEW_STRING(np));
    for (j = 0; j < np; j++) {
	sprintf(s,"%i",j+1);
	SET_STRING_ELT(R_tmp, j, mkChar(s));
    }
    Free(s);
		    
    SET_LEVELS(VECTOR_ELT(R_obj, 2), R_tmp);
    UNPROTECT(1);
    
    UNPROTECT(1);  

    return R_obj;
}

/**/
