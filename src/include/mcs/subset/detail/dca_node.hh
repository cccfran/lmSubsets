// Copyright 2018  Marc Hofmann
//
// This file is part of the 'mcs' library (see
// <https://github.com/marc-hofmann/mcs/>).
//
// 'mcs' is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// 'mcs' is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with 'mcs'.  If not, see <http://www.gnu.org/licenses/>.



#ifndef MCS_SUBSET_DETAIL_DCA_NODE_HH
#define MCS_SUBSET_DETAIL_DCA_NODE_HH



#include <algorithm>  // std::max_element, std::sort
#include <cmath>  // std::pow
#include <iterator>  // std::distance
#include <numeric>  // std::iota
#include <utility>  // std::move, std::swap
#include <vector>



#include "gsl/gsl"  // gsl::span



#include "mcs/core/matrix.hh"

#include "mcs/subset/detail/dca_qrz.hh"
#include "mcs/subset/detail/dca_subset.hh"

#include "mcs/core/lapack.hh"


namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class dca_node
{

    friend void
    swap(dca_node& a, dca_node& b) noexcept
    {
        a.swap(b);
    }

    using matrix = mcs::core::matrix<Scalar>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;

    using dca_qrz = detail::dca_qrz<Scalar>;

    using lapack = mcs::core::lapack<Scalar>;



private:

    std::vector<int> subset_;

    int mark_;

    // matrix rz_mat_;

    std::vector<Scalar> qty_;

    dca_qrz* qrz_;

    int m;


public:

    matrix rz_mat_;    

    matrix X;

    matrix y;

    dca_node(const int root_size) noexcept :
        rz_mat_(root_size + 1, root_size + 1)
    {
        subset_.reserve(root_size);
    }

    // pass in root size and n
    dca_node(const int root_size, const int ay_nrow) noexcept :
        rz_mat_(root_size + 1, root_size + 1),
        qty_(ay_nrow)
    {
        subset_.reserve(root_size);
    }

    // pass in whole qrz object 
    dca_node(const int root_size, 
             const int ay_nrow, 
             dca_qrz* qrz, 
             matrix X_,
             matrix y_) noexcept :
        rz_mat_(root_size + 1, root_size + 1),
        qty_(ay_nrow),
        qrz_(qrz), 
        X(X_),
        y(y_),
        m(ay_nrow)
    {
        subset_.reserve(root_size);
    }



public:

    void
    swap(dca_node& other) noexcept
    {
        subset_.swap(other.subset_);
        std::swap(mark_, other.mark_);
        rz_mat_.swap(other.rz_mat_);
    }



public:

    void
    root(matrix_cspan rz_mat) noexcept
    {
        const int n = rz_mat.ncol() - 1;

        for (int j = 0; j < n; ++j)  subset_.push_back(j);
        mark_ = 0;
        rz_mat_ = rz_mat;

        // std::cout << "root of node mark: " << mark_ << std::endl;
        // std::cout << "root of node cols: " << n << std::endl;
    }


    // pass tri and qty
    void
    root(matrix_cspan rz_mat, std::vector<Scalar> qty) noexcept
    {
        const int n = rz_mat.ncol() - 1;

        for (int j = 0; j < n; ++j)  subset_.push_back(j);
        mark_ = 0;
        rz_mat_ = rz_mat;
        qty_ = qty;


        // std::cout << "root of node mark: " << mark_ << std::endl;
        // std::cout << "root of node cols: " << n << std::endl;
    }



    const std::vector<int>&
    subset() const noexcept
    {
        return subset_;
    }



    int
    size() const noexcept
    {
        return subset_.size();
    }



    int
    mark() const noexcept
    {
        return mark_;
    }



    int
    rank() const noexcept
    {
        return size() - mark();
    }



    Scalar
    rss() const noexcept
    {
        const int n = size();

        return std::pow(rz_mat_(n, n), 2);
    }



    template<typename Function>
    void
    for_each(Function f) const noexcept
    {
        // get all rss of submodel in this node
        const int n = size();
        const int k = mark_;
        // const int m = qty_.size();
        const std::string side = "L";
        const std::string trans = "N";

        gsl::span<const int> s = subset_;
        std::vector<Scalar> aux_work_(n);

        const Scalar* z_ptr = rz_mat_.ptr(n, n);
        Scalar rss = 0;

        for (int j = n; j > k; --j, --z_ptr)
        {
            rss += std::pow(*z_ptr, 2);

            f(s.first(j), rss);
        }

        std::cout << "dca_node: for_each()" << std::endl;
        std::cout << "size of rz: (" << rz_mat_.nrow() << ", " << rz_mat_.ncol() << std::endl;
        // std::cout << "calculating rss: of (size " << n << "):  ";
        // for(auto x : s) std::cout << x << " ";
        // std::cout << std::endl;

        // lapack::trtrs(rz_mat_({0, n-1}, {0, n-1}), rz_mat_({0, n-1}, {n, n}));
        // for(auto x : rz_mat_({0, n-1}, {0, n-1})) std::cout << x << " ";
        //     std::cout << std::endl;
        // matrix_cspan rz_span = rz_mat_({0, n}, {0, n});
        // matrix_cspan qyc = rz_mat_({0, n}, {n, 1});

        matrix_cspan qyc = rz_mat_({0, n}, {n, 1});

        // beta hat
        matrix betahat(qyc);
        lapack::trtrs(rz_mat_({0, n}, {0, n}), betahat);

        std::cout << "BETA" << std::endl;
        for(int i = 0; i < n; ++i) {
            std::cout << betahat(i, 0) << " ";
        }
        std::cout << std::endl;

        // residual y - Xb
        matrix residual(y);
        int beta_front  = 0, prev_front = subset_.front(), prev = subset_.front() - 1, ctr = 0;
        const char transno = 'N';
        for(auto it = subset_.cbegin(); it != subset_.cend(); ++it) {
            std::cout << "it: " << *it << "; prev: " << prev << std::endl;
            if(*it != (prev + 1)) {
                std::cout << "X: " << std::endl;
                for(int i = 0; i < m; i++) {
                    for(int j = prev_front; j < prev_front + ctr; j++) {
                        std::cout << X(i,j) << "\t";
                    }
                    std::cout << std::endl;
                }
                lapack::gemm(&transno, &transno, m, 1, ctr, -1.0, 
                    X.ptr(0, prev_front), m, betahat.ptr(beta_front, 0), ctr, 1.0, residual.base(), m);
                std::cout << std::endl << "AFTER: " << ctr << std::endl;
                for (int j = 0; j < m; j++) {
                    std::cout << residual(j,0) << " ";
                }
                std::cout << std::endl;

                beta_front += ctr; prev_front = prev = *it; ctr = 1;
            } else {++prev; ++ctr;}

            if(next(it) == subset_.cend()) {
                std::cout << "!X: " << std::endl;
                for(int i = 0; i < m; i++) {
                    for(int j = prev_front; j < prev_front + ctr; j++) {
                        std::cout << X(i,j) << "\t";
                    }
                    std::cout << std::endl;
                }
                std::cout << "Beta: " << beta_front << betahat(beta_front, 0)  << std::endl;
                lapack::gemm(&transno, &transno, m, 1, ctr, -1.0, 
                    X.ptr(0, prev_front), m, betahat.ptr(beta_front, 0), ctr, 1.0, residual.base(), m);
                std::cout << std::endl << "AFTER: " << ctr << std::endl;
                for (int j = 0; j < m; j++) {
                    std::cout << residual(j,0) << " ";
                }
                std::cout << std::endl;
            }
        }

        std::cout << std::endl << "res" << std::endl;
        for (int j = 0; j < m; j++) {
            std::cout << residual(j,0) << " ";
        }
        std::cout << std::endl;
    }



    void
    drop_column(
        const int mark,
        dca_node& result,
        const dca_qrz& qrz
    ) const noexcept
    {
        const int n = size();
        const int k = mark;

        matrix_cspan rz_span = rz_mat_({0, n+1}, {0, n+1});

        // init next node (result)
        // subset_, mark_, rz_mat_
        // drop k-th
        // take care of the subset
        dca_subset::drop_column(subset_, k, result.subset_);
        result.mark_ = k;
        // drop k-th column of R and Givens rotation
        std::cout << "## NODE: drop column " << mark << std::endl;
        // std::cout << "result.rz_mat_: " << result.rz_mat_.ldim() << std::endl;
        qrz.drop_column(rz_span, k, result.rz_mat_);
    }



    void
    preorder_complete(
        dca_node& result,
        const dca_qrz& qrz,
        std::vector<Scalar>& aux_1,
        std::vector<int>& aux_2
    ) const noexcept
    {
        const int n = size();
        const int k = mark_;
        const int p = n - k;

        matrix_cspan rz_span = rz_mat_({0, n+1}, {0, n+1});

        qrz.column_bounds(rz_span, k, aux_1);

        std::iota(aux_2.begin(), aux_2.begin() + p, 0);
        std::sort(aux_2.begin(), aux_2.begin() + p,
                  [&aux_1](const int i, const int j) -> bool {
                      return aux_1[i] > aux_1[j];
                  });

        dca_subset::permute_complete(subset_, k, aux_2, result.subset_);
        result.mark_ = k;
        qrz.permute_complete(rz_span, k, aux_2, result.rz_mat_);
    }



    void
    preorder_partial_1(
        dca_node& result,
        const dca_qrz& qrz,
        std::vector<Scalar>& aux_1
    ) noexcept
    {
        const int n = size();
        const int k = mark_;
        const int p = n - k;

        matrix_cspan rz_span = rz_mat_({0, n+1}, {0, n+1});

        qrz.column_bounds(rz_span, k, aux_1);

        const auto max = std::max_element(aux_1.begin(), aux_1.begin() + p);
        const int j = std::distance(aux_1.begin(), max);

        dca_subset::permute_partial_1(subset_, k, j);
        qrz.permute_partial_1(rz_span, k, j);

        swap(result);
    }



    void
    preorder_partial_2(
        dca_node& result,
        const dca_qrz& qrz,
        std::vector<Scalar>& aux_1
    ) noexcept
    {
        const int n = size();
        const int k = mark_;
        const int p = n - k;

        matrix_cspan rz_span = rz_mat_({0, n+1}, {0, n+1});

        qrz.column_bounds(rz_span, k, aux_1);

        const auto max = std::max_element(aux_1.begin(), aux_1.begin() + p);
        const int j = std::distance(aux_1.begin(), max);

        dca_subset::permute_partial_1(subset_, k, j);
        qrz.permute_partial_1(rz_span, k, j);

        swap(result);
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
