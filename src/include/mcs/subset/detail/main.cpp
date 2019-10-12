#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>

#include "state.h"
#include "node.h"


using namespace std;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List all_subset (
	NumericMatrix Xr,
	NumericVector yr,
	int mark
)
{
	const int n = Xr.nrow(), k = Xr.ncol();
	arma::mat X(Xr.begin(), n, k, false); 
	arma::colvec y(yr.begin(), yr.size(), false);

	state state(X, y, mark);

	return Rcpp::List::create(_["y"] = y);
}