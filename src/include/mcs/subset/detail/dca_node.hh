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

// #include <ctime>


#include "gsl/gsl"  // gsl::span



#include "mcs/core/matrix.hh"

#include "mcs/subset/detail/dca_qrz.hh"
#include "mcs/subset/detail/dca_subset.hh"

#include "mcs/core/lapack.hh"

#include "mcs/core/model.hh"

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

    using model = mcs::core::model<Scalar>;



private:

    

    int mark_;

    // matrix rz_mat_;

    std::vector<Scalar> qty_;

    dca_qrz* qrz_;

    int m;

    std::vector<int> subset_;

    std::vector<model> models_;

public:

    typename decltype(models_)::iterator cur_model_;

    matrix rz_mat_;

    matrix qt_mat_;

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
        qt_mat_(root_size + 1, ay_nrow),
        qty_(ay_nrow),
        qrz_(qrz), 
        X(X_),
        y(y_),
        m(ay_nrow)
    {
        subset_.reserve(root_size);    
        models_.reserve(root_size);
    }



public:

    void
    swap(dca_node& other) noexcept
    {
        subset_.swap(other.subset_);
        std::swap(mark_, other.mark_);
        rz_mat_.swap(other.rz_mat_);
        qt_mat_.swap(other.qt_mat_);
        std::swap(models_, other.models_);
        std::swap(cur_model_, other.cur_model_);
        // std::cout << "SWAPPPPPPPP" << std::endl;
    }



public:

    void
    root(matrix_cspan rz_mat) noexcept
    {
        const int n = rz_mat.ncol() - 1;

        for (int j = 0; j < n; ++j)  subset_.push_back(j);
        mark_ = 0;
        rz_mat_ = rz_mat;
        qrz_->qt(qt_mat_);
        // qt_mat_ = qrz_->get_qt();
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

        gsl::span<const int> s = subset_;

        const Scalar* z_ptr = rz_mat_.ptr(n, n);
        Scalar rss = 0;
        // typename std::vector<model>::const_iterator model_ptr = models_.begin();
        auto model_ptr = models_.begin();
        
        for (int j = n; j > k; --j, model_ptr++)
        // for (int j = n; j > k; --j, --z_ptr)
        {
            // rss += std::pow(*z_ptr, 2);
            // f(s.first(j), rss);
            // std::cout << "model size: " << model_ptr->size() << std::endl;
            // std::cout << "maxt: " << model_ptr->maxt() << std::endl;

            f(s.first(j), model_ptr->maxt());
        }

        // std::cout << "dca_node: for_each()" << std::endl;
        // std::cout << "size of rz: (" << rz_mat_.nrow() << ", " << rz_mat_.ncol() << ")" << std::endl;
        // std::cout << "mark_ in node: " << k << std::endl;

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
        const int m = qt_mat_.ncol();

        matrix_cspan rz_span = rz_mat_({0, n+1}, {0, n+1});
        matrix_cspan qt_span = qt_mat_({0, n}, {0, m});

        // init next node (result)
        // subset_, mark_, rz_mat_
        // drop k-th
        // take care of the subset_
        dca_subset::drop_column(subset_, k, result.subset_);
        result.mark_ = k;
        // drop k-th column of R and Givens rotation
        // std::cout << "## NODE: drop column " << mark << std::endl;
        // std::cout << "qt_mat_: " << qt_span.ldim() << std::endl;
        // std::cout << "qt_mat_ nrow: " << qt_span.nrow() << std::endl;
        qrz.drop_column(rz_span, k, result.rz_mat_, qt_span, result.qt_mat_);

        // qrz.drop_column(rz_span, k, result.rz_mat_);

    }


    void
    get_t(
        const int mark,
        const dca_qrz& qrz,
        const matrix& X,
        const matrix& y 
    ) noexcept
    {
        const int n = size();

        // clock_t begin = std::clock();
        // clock_t tmp;
        
        // const int m = y.nrow();

        // matrix residual(m, 1);

        // models_.reserve(n);
        models_.clear();
        for(int j = n; j > mark_; j--) {
            models_.emplace_back(j);
        }
        cur_model_ = models_.begin();

        for(int j = n; j > mark_; j--, cur_model_++) {
            get_beta(j);
            // tmp = std::clock();
            // std::cout << "inner_get_beta: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
            // begin = tmp;
            get_sds(j);

            // tmp = std::clock();
            // std::cout << "inner_get_sd: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
            // begin = tmp;
        }

        // tmp = std::clock();
        // std::cout << "inner_get_t_end: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
        // begin = tmp;

    }


    void
    get_beta(const int model_size) noexcept
    {   
        int n = size();
        // int cur_size = 0;

        matrix betahat(rz_mat_({0, model_size}, {n, 1}));
        
        lapack::trtrs(rz_mat_({0, model_size}, {0, model_size}), betahat);
        
        cur_model_->set_beta(betahat);
        // std::cout << "In get_beta() 3: " << std::endl;
        // for(auto i = subset_.begin(); cur_size < model_size && i != subset_.end(); i++, ++cur_size) {
        //     std::cout << *i << "\t";
        // }
        // std::cout << std::endl;
        // for(int i = 0; i < model_size; i++) {
        //     std::cout <<  cur_model_->beta()(i,0) << "\t";
        // }
        // std::cout << std::endl;
    }


    void
    get_residual(
        // const matrix& X,
        // const matrix& y,
        double* residual_mat,
        const int model_size
    )
    {

        matrix residual(y);
        int beta_front  = 0, prev_front = subset_.front(), prev = subset_.front() - 1;
        int cur_size = 0, ctr = 0;
        const char transno = 'N';

        for(auto it = subset_.cbegin(); cur_size < model_size; ++it, ++cur_size) {
            if(*it != (prev + 1)) {  
                lapack::gemm(lapack::no_trans, lapack::no_trans, m, 1, ctr, -1.0, 
                    X.ptr(0, prev_front), m, cur_model_->beta().ptr(beta_front, 0), ctr, 1.0, residual.base(), m);
                beta_front += ctr; prev_front = prev = *it; ctr = 1;
            } else {++prev; ++ctr;}

            if((cur_size + 1) == model_size) {
                lapack::gemm(lapack::no_trans, lapack::no_trans, m, 1, ctr, -1.0, 
                    X.ptr(0, prev_front), m, cur_model_->beta().ptr(beta_front, 0), ctr, 1.0, residual.base(), m);
            }
        }

        for(int i = 0; i < m; i++) {
            // residual_mat(i, i) = residual(i, 0);
            *(residual_mat + i * (m + 1)) = residual(i, 0);

        }

        // std::cout << std::endl << "res" << std::endl;
        // for(int i = 0; i < m; i++) {
        //     for (int j = 0; j < m; j++) {
        //         std::cout <<  *(residual_mat + i + j*m) << "\t";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;
    }

    void
    get_sds(
        // const matrix& X,
        // const matrix& y,
        const int model_size
    ) noexcept
    {   
        // const int m = y.nrow();

        clock_t begin = std::clock();
        clock_t tmp;

        // matrix residual_mat_2(m, m);
        // matrix qtr(m, m);
        // std::cout << "In get_sds() 1: " << std::endl;
        double residual_mat[m*m], qtr[model_size*m];
        // std::cout << "In get_beta() 2: " << std::endl;
        // tmp = std::clock();
        // std::cout << "\tget_t_init: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
        // begin = tmp;

        std::fill(residual_mat, residual_mat + m*m, 0.0);

        // tmp = std::clock();
        // std::cout << "\tget_t_fill: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
        // begin = tmp;

        matrix sds(model_size, 1);
        

        // tmp = std::clock();
        // std::cout << "\tget_t_init_sds: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
        // begin = tmp;


        // qtr = 0;        
        // residual_mat = 0;
        // get_residual(X, y, residual_mat, model_size);
        get_residual(residual_mat, model_size);

        // tmp = std::clock();
        // std::cout << "\tget_t_res: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
        // begin = tmp;

        // std::cout << "Q^T" << std::endl;
        // for (int i = 0; i < qt_mat_.ldim(); i++) {
        //     for (int j = 0; j < m; j++) {
        //         std::cout << qt_mat_(i,j) << "\t";
        //     }
        //     std::cout << std::endl;
        // }

        // multiply Q^T R
        // lapack::ormqr(lapack::left, lapack::trans, m, qrz_->get_qrr(), qrz_->get_tau(), residual_mat, aux_work_);
        // lapack::gemm(lapack::no_trans, lapack::no_trans, model_size, m, m, 1.0, 
        //             qt_mat_.base(), m, residual_mat.base(), m, 0.0, qtr.base(), m);
        // lapack::gemm(lapack::no_trans, lapack::no_trans, 1.0, qt_mat_, residual_mat, 0.0, qtr);
        lapack::gemm(lapack::no_trans, lapack::no_trans, model_size, m, m, 1.0,
            qt_mat_.base(), qt_mat_.ldim(), residual_mat, m, 0.0, qtr, m);
        // tmp = std::clock();
        // std::cout << "\tget_t_qtR: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
        // begin = tmp;

        // std::cout << "Q^T residual_mat" << std::endl;
        // for (int i = 0; i < model_size; i++) {
        //     for (int j = 0; j < m; j++) {
        //         std::cout << *(qtr + i + j * m) << "\t";
        //     }
        //     std::cout << std::endl;
        // }
        // lapack::trtrs(rz_mat_({0, model_size}, {0, model_size}), qtr);
        lapack::trtrs(lapack::upper, lapack::no_trans, lapack::no_trans, 
            model_size, m,
            rz_mat_({0, model_size}, {0, model_size}),
            qtr, m);

        // tmp = std::clock();
        // std::cout << "\tget_t_tri: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
        // begin = tmp;

        for (int i = 0; i < model_size; i++) {
            double tmp = 0; 
            for (int j = 0; j < m; j++) {
                tmp += std::pow(*(qtr + i + j * m), 2.0);
            }    
            sds(i, 0) = std::sqrt(tmp);
            // std::cout << tmp << "\t";
        }

        // std::cout << std::endl;

        // tmp = std::clock();
        // std::cout << "\tget_t_pow: " << double(tmp - begin) / CLOCKS_PER_SEC << std::endl;
        // begin = tmp;


        // std::cout << "sd" << std::endl;
        // for (int i = 0; i < model_size; i++) {
        //     std::cout << sds(i, 0) << "\t";
        // }
        // std::cout << std::endl;

        cur_model_->set_sds(sds);

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
