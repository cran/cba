
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP ccfkms(SEXP R_x, SEXP R_p, SEXP R_par, SEXP R_max_iter,
		   SEXP R_opt_std, SEXP R_debug);
extern SEXP cluster_dist(SEXP R_x, SEXP R_beta);
extern SEXP as_dummy(SEXP R_x);
extern SEXP gknn(SEXP R_x, SEXP R_y, SEXP R_k, SEXP R_l, SEXP R_break_ties,
		 SEXP R_use_all, SEXP R_prob);
extern SEXP order_optimal(SEXP R_dist, SEXP R_merge);
extern SEXP order_length(SEXP R_dist, SEXP R_order);
extern SEXP order_greedy(SEXP R_dist);
extern SEXP lminter(SEXP R_x, SEXP R_block_size, SEXP R_nbin);
extern SEXP proximus(SEXP R_mat, SEXP R_max_radius, SEXP R_min_size, 
		     SEXP R_min_retry, SEXP R_max_iter, SEXP R_debug);
extern SEXP rockLink(SEXP R_x, SEXP R_beta);
extern SEXP rockMerge(SEXP R_x, SEXP R_n, SEXP R_theta, SEXP R_debug);
extern SEXP rockClass(SEXP R_x, SEXP R_l, SEXP R_beta, SEXP R_theta);
extern SEXP sdists(SEXP R_x, SEXP R_y, SEXP R_method, SEXP R_weight, 
		   SEXP R_pairwise);
extern SEXP sdists_transcript(SEXP R_x, SEXP R_y, SEXP R_method,
			      SEXP R_weight, SEXP R_table);
extern SEXP sdists_graph(SEXP x);
extern SEXP sdists_align(SEXP R_x, SEXP R_y, SEXP t);
extern SEXP stress(SEXP R_x, SEXP R_r, SEXP R_c, SEXP R_type);
extern SEXP stress_dist(SEXP R_x, SEXP R_r, SEXP R_c, SEXP R_bycol,
			SEXP R_type);
extern SEXP orderTSP(SEXP x, SEXP t);

static const R_CallMethodDef CallEntries[] = {
    {"R_ccfkms",	    (DL_FUNC) ccfkms,		 6},
    {"R_cluster_dist",	    (DL_FUNC) cluster_dist,	 2},
    {"R_as_dummy",	    (DL_FUNC) as_dummy,		 1},
    {"R_gknn",		    (DL_FUNC) gknn,		 7},
    {"R_order_optimal",	    (DL_FUNC) order_optimal,	 2},
    {"R_order_length",	    (DL_FUNC) order_length,	 2},
    {"R_order_greedy",	    (DL_FUNC) order_greedy,	 1},
    {"R_lminter",	    (DL_FUNC) lminter,		 3},
    {"R_proximus",	    (DL_FUNC) proximus,		 6},
    {"R_rockLink",	    (DL_FUNC) rockLink,		 2},
    {"R_rockMerge",	    (DL_FUNC) rockMerge,	 4},
    {"R_rockClass",	    (DL_FUNC) rockClass,	 4},
    {"R_sdists",	    (DL_FUNC) sdists,		 5},
    {"R_sdists_transcript", (DL_FUNC) sdists_transcript, 5},
    {"R_sdists_graph",	    (DL_FUNC) sdists_graph,	 1},
    {"R_sdists_align",	    (DL_FUNC) sdists_align,	 3},
    {"R_stress",	    (DL_FUNC) stress,		 4},
    {"R_stress_dist",	    (DL_FUNC) stress_dist,	 5},
    {"R_orderTSP",	    (DL_FUNC) orderTSP,		 2},
    {NULL, NULL, 0}
};

void R_init_cba(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

