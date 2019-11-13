#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include <iostream>

#include <string>
#include <tuple>  // std::tie
#include <vector>

#include "cxx11.inc"
#include "r_interrupt.inc"


#include "mcs/core/lapack.hh"
#include "mcs/core/matrix.hh"

#include "mcs/subset/aic.hh"
#include "mcs/subset/abba.hh"
#include "mcs/subset/bba.hh"
#include "mcs/subset/dca.hh"
#include "mcs/subset/hbba.hh"
#include "mcs/subset/subset.hh"
#include "mcs/subset/table.hh"

#include "mcs/util/misc.hh"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

using matrix_cspan = mcs::core::matrix<const double&>;

using mcs::subset::subset_all;
using mcs::subset::abba_all;
using mcs::subset::hbba_all;
using mcs::subset::bba_all;
using mcs::subset::dca_all;
using mcs::subset::dca_all_boot;

using mcs::subset::subset_best;
using mcs::subset::abba_best;
using mcs::subset::hbba_best;
using mcs::subset::bba_best;
using mcs::subset::dca_best;

using mcs::subset::aic;

using mcs::subset::table_all;
using mcs::subset::table_best;

using mcs::util::to_ordinal;


namespace {

const std::string algo_dflt = "DFLT";
const std::string algo_abba = "abba";
const std::string algo_bba = "bba";
const std::string algo_dca = "dca";
const std::string algo_hbba = "hbba";

}

// [[Rcpp::export]]
SEXP
lmSubsets(
    SEXP r_algo,
    SEXP r_xy,
    SEXP r_mark,
    SEXP r_tau,
    SEXP r_nbest,
    SEXP r_prad,
    SEXP r_nboot
)
{
    int protect_cnt = 0;


    if (!Rf_isNull(r_algo) && !Rf_isString(r_algo))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'algo' must be a character string");
    }

    const auto algo = [&r_algo]() -> std::string {
        if (!Rf_isNull(r_algo))
            return CHAR(STRING_ELT(r_algo, 0));
        return algo_dflt;
    }();


    if (!Rf_isMatrix(r_xy))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'xy' must be a numeric matrix");
    }

    if (!Rf_isReal(r_xy))
    {
        r_xy = Rf_protect(Rf_coerceVector(r_xy, REALSXP));
        ++protect_cnt;
    }

    const int* const xy_dim =
        INTEGER(Rf_coerceVector(Rf_getAttrib(r_xy, R_DimSymbol), INTSXP));
    const int m = xy_dim[0];
    const int n = xy_dim[1] - 1;

    if (m <= n)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'xy' (%d x %d) must be a tall (or square) matrix", m, n + 1);
    }

    matrix_cspan ay_mat(m, n + 1, REAL(r_xy));

    if (!Rf_isNumeric(r_mark))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'mark' must be numeric");
    }

    const int mark = Rf_asInteger(r_mark);

    if (mark < 0)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'mark' [%d] must be a non-negative integer", mark);
    }


    if (!Rf_isNumeric(r_tau))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'tau' must be numeric vector");
    }

    if (LENGTH(r_tau) != n)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'tau' [%d] must be of length %d", LENGTH(r_tau), n);
    }

    if (!Rf_isReal(r_tau))
    {
        r_tau = Rf_protect(Rf_coerceVector(r_tau, REALSXP));
        ++protect_cnt;
    }

    const std::vector<double> tau(REAL(r_tau), REAL(r_tau) + n);


    if (!Rf_isNumeric(r_nbest))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'nbest' must be numeric");
    }

    const int nbest = Rf_asInteger(r_nbest);

    if (nbest < 1)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'nbest' [%d] must be positive integer", nbest);
    }


    if (!Rf_isNumeric(r_prad))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'prad' must be numeric");
    }

    const int prad = Rf_asInteger(r_prad);

    if (prad < 0)
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'prad' [%d] must be a non-negative integer", prad);
    }

    if (!Rf_isNumeric(r_nboot))
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'nboot' must be numeric");
    }

    const int nboot = Rf_asInteger(r_nboot);

    r_interrupt_setup();

    // Main function 
    // build tree

    table_all<double> tab;
    int node_cnt = -1;

    if (algo == algo_dflt)
    {
        tab = subset_all<double>(ay_mat, mark, tau, nbest, prad);
    }
    else if (algo == algo_abba)
    {
        std::tie(tab, node_cnt) =
            abba_all<double>(ay_mat, mark, tau, nbest, prad);
    }
    else if (algo == algo_hbba)
    {
        std::tie(tab, node_cnt) =
            hbba_all<double>(ay_mat, mark, tau, nbest, prad);
    }
    else if (algo == algo_bba)
    {
        std::tie(tab, node_cnt) =
            bba_all<double>(ay_mat, mark, nbest, prad);
    }
    else if (algo == algo_dca)
    {
    	// init dca_state 
    	// loop through dca_state
    	// cout << "DCA_all" << endl;
        // std::tie(tab, node_cnt) =
        //     dca_all<double>(ay_mat, mark, nbest, prad);
        std::tie(tab, node_cnt) =
            dca_all_boot<double>(ay_mat, mark, nbest, prad, nboot);
    } 
    else
    {
        Rf_unprotect(protect_cnt);
        Rf_error("'algo' [%s]: unexpected value", algo.c_str());
    }

    // print the table of residual
    for (int size = 1; size <= n; ++size)
    {
        for (int best = 1; best <= nbest; ++best)
        {
            const auto& res = tab[size - 1][best - 1];

            // std::cout << res.key() << "\t" ;
        }
        // std::cout << std::endl;
    } 

    // print subset
    // for (int size = 1; size <= n; ++size)
    // {
    // 	cout << "size: " << size << endl;
    //     for (int best = 1; best <= nbest; ++best)
    //     {	
    //     	cout <<  "nbest: " << best << endl;
    //         const auto& res = tab[size - 1][best - 1];

    //         for (auto it=res.cbegin(); it!=res.cend(); ++it) {
    //         	cout << ' ' << *it;
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    const int ncase = n * nbest;
 

    SEXP r_submodel_size = Rf_protect(Rf_allocVector(INTSXP, ncase));
    ++protect_cnt;

    SEXP r_submodel_best = Rf_protect(Rf_allocVector(INTSXP, ncase));
    ++protect_cnt;

    SEXP r_submodel_rss = Rf_protect(Rf_allocVector(REALSXP, ncase));
    ++protect_cnt;


    SEXP r_subset_dim = Rf_protect(Rf_allocVector(INTSXP, 2));
    ++protect_cnt;

    INTEGER(r_subset_dim)[0] = ncase;
    INTEGER(r_subset_dim)[1] = n;

    SEXP r_subset = Rf_protect(Rf_allocArray(LGLSXP, r_subset_dim));
    ++protect_cnt;


    for (int size = 1; size <= n; ++size)
    {
        for (int best = 1; best <= nbest; ++best)
        {
            const int i = (size - 1) * nbest + (best - 1);

            INTEGER(r_submodel_size)[i] = size;
            INTEGER(r_submodel_best)[i] = best;

            const auto& res = tab[size - 1][best - 1];

            if (res)
            {
                REAL(r_submodel_rss)[i] = res.key();

                for (int j = 0; j < n; ++j)
                {
                    LOGICAL(r_subset)[i + j * ncase] = FALSE;
                }

                for (int k = 0; k < res.size(); ++k)
                {
                    const int j = res[k];

                    LOGICAL(r_subset)[i + j * ncase] = TRUE;
                }
            }
            else
            {
                REAL(r_submodel_rss)[i] = NA_REAL;

                for (int j = 0; j < n; ++j)
                {
                    LOGICAL(r_subset)[i + j * ncase] = NA_LOGICAL;
                }
            }
        }
    }
 

    SEXP r_submodel_names = Rf_protect(Rf_allocVector(STRSXP, 3));
    ++protect_cnt;

    SET_STRING_ELT(r_submodel_names, 0, Rf_mkChar("SIZE"));
    SET_STRING_ELT(r_submodel_names, 1, Rf_mkChar("BEST"));
    SET_STRING_ELT(r_submodel_names, 2, Rf_mkChar("RSS"));

    SEXP r_submodel_row_names = Rf_protect(Rf_allocVector(INTSXP, 2));
    ++protect_cnt;

    INTEGER(r_submodel_row_names)[0] = NA_INTEGER;
    INTEGER(r_submodel_row_names)[1] = -ncase;

    SEXP r_submodel = Rf_protect(Rf_allocVector(VECSXP, 3));
    ++protect_cnt;

    Rf_setAttrib(r_submodel, R_ClassSymbol,
                 Rf_ScalarString(Rf_mkChar("data.frame")));
    Rf_setAttrib(r_submodel, R_NamesSymbol, r_submodel_names);
    Rf_setAttrib(r_submodel, R_RowNamesSymbol, r_submodel_row_names);

    SET_VECTOR_ELT(r_submodel, 0, r_submodel_size);
    SET_VECTOR_ELT(r_submodel, 1, r_submodel_best);
    SET_VECTOR_ELT(r_submodel, 2, r_submodel_rss);


    SEXP r_ans_names = Rf_protect(Rf_allocVector(STRSXP, 4));
    ++protect_cnt;

    SET_STRING_ELT(r_ans_names, 0, Rf_mkChar("submodel"));
    SET_STRING_ELT(r_ans_names, 1, Rf_mkChar("subset"));
    SET_STRING_ELT(r_ans_names, 2, Rf_mkChar(".interrupted"));
    SET_STRING_ELT(r_ans_names, 3, Rf_mkChar(".nodes"));

    SEXP r_ans = Rf_protect(Rf_allocVector(VECSXP, 4));
    ++protect_cnt;

    Rf_setAttrib(r_ans, R_NamesSymbol, r_ans_names);

    SET_VECTOR_ELT(r_ans, 0, r_submodel);
    SET_VECTOR_ELT(r_ans, 1, r_subset);
    SET_VECTOR_ELT(r_ans, 2, Rf_ScalarLogical(r_interrupt_flag()));
    SET_VECTOR_ELT(r_ans, 3, Rf_ScalarInteger(node_cnt));


    Rf_unprotect(protect_cnt);

    return r_ans;
}