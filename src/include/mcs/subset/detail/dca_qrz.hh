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



#ifndef MCS_SUBSET_DETAIL_DCA_QRZ_HH
#define MCS_SUBSET_DETAIL_DCA_QRZ_HH



#include <algorithm>  // std::copy_n, std::fill_n
#include <cstdlib>  // std::abs
#include <utility>  // std::swap
#include <vector>



#include "mcs/core/givens.hh"
#include "mcs/core/lapack.hh"
#include "mcs/core/matrix.hh"



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class dca_qrz
{

    using matrix = mcs::core::matrix<Scalar>;

    using matrix_span = mcs::core::matrix<Scalar&>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;

    using givens = mcs::core::givens<Scalar>;

    using lapack = mcs::core::lapack<Scalar>;



private:

    mutable std::vector<Scalar> aux_work_;

    mutable std::vector<Scalar> aux_tau_;

    mutable std::vector<givens> aux_givens_;

    mutable matrix qty_;

    mutable matrix qrq;


public:

    dca_qrz(const int root_size) noexcept :
        aux_work_(root_size + 1),
        aux_tau_(root_size + 1),
        aux_givens_(root_size + 1)
    {
    }

    dca_qrz(const int root_size, const int ay_nrow) noexcept :
        aux_work_(root_size + 1),
        aux_tau_(root_size + 1),
        aux_givens_(root_size + 1),
        qty_(ay_nrow, 1),
        qrq(ay_nrow, root_size)
    {
    }

public:

    // QR ay_mat and get tau and qty_
    matrix
    rz(matrix_cspan ay_mat) noexcept
    {
        const int n = ay_mat.ncol() - 1;
        const int m = ay_mat.nrow() ;
        const int dy = 1;
        const std::string side = "L";
        const std::string trans = "T";

        matrix rz_tmp = ay_mat;
        lapack::geqr2(rz_tmp, aux_tau_, aux_work_);
        qrq = rz_tmp({0, m}, {0, n});

        // std::cout << "dca_qrz: rz()" << std::endl;
        // std::cout << "Checkthings ********************" << std::endl;
        for(int i = 0; i < m; ++i) {
            qty_(i,0) = ay_mat(i, n);
            // if (i < 10) std::cout << qty_(i,0) << " ";
        }
        // std::cout << std::endl;

        // std::cout << "QR matrix ********************" << std::endl;
        // for(int i = 0; i < m; i++) {
        //     for (int j = 0; j < n; j++) {
        //         std::cout << qrq(i, j) << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;

        lapack::orgqr(n, qrq, aux_tau_, aux_work_);

        // std::cout << "Q matrix ********************" << std::endl;
        // for(int i = 0; i < m; i++) {
        //     for (int j = 0; j < n; j++) {
        //         std::cout << qrq(i, j) << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;

        // lapack::ormqr(side, trans, m, qrq, aux_tau_, qty_, aux_work_);
        // lapack::gemm(lapack::trans, lapack::no_trans, 1.0, qrq, qty_, 0.0, qty_);

        // std::cout << "QTy ********************" << std::endl;
        // for(int i = 0; i < m; ++i) {
        //     std::cout << qty_(i,0) << " ";
        // }
        // std::cout << std::endl;


        return rz_tmp({0, n + 1}, {0, n + 1});
    }

    std::vector<Scalar>
    get_tau() {
        return aux_tau_;
    }

    matrix
    get_qty() {
        return qty_;
    }

    matrix
    get_qrq() {
        return qrq;
    }



    void
    column_bounds(
        matrix_cspan rz_mat,
        const int mark,
        std::vector<Scalar>& out
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;
        const int p = n - k;

        column_bounds(p, rz_mat.ptr(k, k), rz_mat.ldim(), out.data(),
                      aux_givens_.data());
    }



    void
    permute_complete(
        matrix_cspan rz_mat,
        const int mark,
        const std::vector<int>& pos,
        matrix_span out_mat
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;
        const int p = n - k;

        permute_complete(p, rz_mat.ptr(k, k), rz_mat.ldim(), pos.data(),
                         out_mat.ptr(k, k), out_mat.ldim(), aux_tau_.data(),
                         aux_work_.data());
    }


    void
    permute_partial_1(
        matrix_span rz_mat,
        const int mark,
        const int pos
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;

        permute_partial_1(n - k, rz_mat.ptr(k, k), rz_mat.ldim(), pos - k);
    }



    void
    permute_partial_2(
        matrix_span rz_mat,
        const int mark,
        const int pos
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;

        permute_partial_2(n - k, rz_mat.ptr(k, k), rz_mat.ldim(), pos - k);
    }


    void
    drop_column(
        matrix_cspan rz_mat,
        const int mark,
        matrix_span out_mat
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;

        
        // std::cout << "### QRZ_: drop column public " << k << std::endl;
        // // std::cout << "### rz_mat: (" << k << ", " << k << "): " << rz_mat(k, k) << std::endl;
        // for (int i = 0; i < rz_mat.nrow(); ++i) {
        //     for (int j = 0; j < rz_mat.ncol(); ++j ) {
        //         std::cout << *(rz_mat.ptr(i, j)) << "\t" ;
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "### result.rz_mat_: " << out_mat.ldim() << std::endl;
        // for (int i = 0; i < out_mat.nrow(); ++i) {
        //     for (int j = 0; j < out_mat.ncol(); ++j ) {
        //         std::cout << *(out_mat.ptr(i, j)) << "\t" ;
        //     }
        //     std::cout << std::endl;
        // }

        drop_column(n - k, rz_mat.ptr(k, k), rz_mat.ldim(), out_mat.ptr(k, k),
                    out_mat.ldim());
        // fill in firs˜t (k-1) row
        for(int i = 0; i < k; i++) {
            for(int j = 0; j < n; j++) {
                if(j < k) out_mat.elem(i,j) = rz_mat.elem(i, j);
                else out_mat.elem(i,j) = rz_mat.elem(i, j+1);
            }
        }

        // std::cout << "### result.rz_mat_: " << out_mat.ldim() << std::endl;
        // for (int i = 0; i < out_mat.nrow(); ++i) {
        //     for (int j = 0; j < out_mat.ncol(); ++j ) {
        //         std::cout << *(out_mat.ptr(i, j)) << "\t" ;
        //     }
        //     std::cout << std::endl;
        // }
    }


    void
    drop_column(
        matrix_cspan rz_mat,
        const int mark,
        matrix_span out_mat,
        matrix_cspan qt_mat,
        matrix_span qt_out_mat
    ) const noexcept
    {
        const int n = rz_mat.ncol() - 1;
        const int k = mark;

        
        // std::cout << "### QRZ_: drop column public " << k << " ;n: " << n << std::endl;
        // std::cout << "### rz_mat: (" << k << ", " << k << "): " << rz_mat(k, k) << std::endl;
        // for (int i = 0; i < rz_mat.nrow(); ++i) {
        //     for (int j = 0; j < rz_mat.ncol(); ++j ) {
        //         std::cout << *(rz_mat.ptr(i, j)) << "\t" ;
        //     }
        //     std::cout << std::endl;
        // }

        // std::cout << "### qt_mat ldim: " << qt_mat.ldim() << std::endl;
        // std::cout << "### qt_mat nrow: " << qt_mat.nrow() << std::endl;
        // std::cout << "### qt_out_mat ldim: " << qt_out_mat.ldim() << std::endl;
        // std::cout << "### qt_out_mat nrow: " << qt_out_mat.nrow() << std::endl;

        // std::cout << "### qt_mat_: " << std::endl;
        // for (int i = 0; i < qt_mat.nrow(); ++i) {
        //     for (int j = 0; j < qt_mat.ncol(); ++j ) {
        //         std::cout << *(qt_mat.ptr(i, j)) << "\t" ;
        //     }
        //     std::cout << std::endl;
        // }

        drop_column(n - k, rz_mat.ptr(k, k), rz_mat.ldim(), 
                    out_mat.ptr(k, k), out_mat.ldim(), 
                    qt_mat.ptr(k, 0), qt_out_mat.ptr(k, 0), qt_mat.ldim(), qt_mat.ncol());
        // fill in firs˜t (k-1) row
        for(int i = 0; i < k; i++) {
            for(int j = 0; j < n; j++) {
                if(j < k) out_mat.elem(i,j) = rz_mat.elem(i, j);
                else out_mat.elem(i,j) = rz_mat.elem(i, j+1);
            }
        }

        for(int i = 0; i < k; i++) {
            for(int j = 0; j < qt_mat.ncol(); j++) {
                qt_out_mat(i, j) = qt_mat(i, j);
            }
        }

        // std::cout << "### result.qt_out_mat_: " << std::endl;
        // for (int i = 0; i < qt_out_mat.nrow(); ++i) {
        //     for (int j = 0; j < qt_out_mat.ncol(); ++j ) {
        //         std::cout << *(qt_out_mat.ptr(i, j)) << "\t" ;
        //     }
        //     std::cout << std::endl;
        // }

        // std::cout << "### result.rz_mat_: " << out_mat.ldim() << std::endl;
        // for (int i = 0; i < out_mat.nrow(); ++i) {
        //     for (int j = 0; j < out_mat.ncol(); ++j ) {
        //         std::cout << *(out_mat.ptr(i, j)) << "\t" ;
        //     }
        //     std::cout << std::endl;
        // }
    }


    matrix& 
    qt()
    {
        const int m = qrq.ncol();
        const int n = qrq.nrow();

        matrix qrqt(m, n);
        for(int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                qrqt(i, j) = qrq(j, i);
            }
        }

        return qrqt;
    }


private:

    static void
    drop_column(
        int n,
        const Scalar* rz,
        const int ldrz,
        Scalar* out,
        const int ldout
    ) noexcept
    {
        // std::cout << "rz matrix" << std::endl;
        // for (int i = 0; i < n; ++i) {
        //     for (int j = 0; j < n; ++j ) {
        //         std::cout << *(rz + j*n + i) << "\t" ;
        //     }
        //     std::cout << std::endl;
        // }

        // std::cout << "### QRZ_: drop column private" << std::endl;
        // std::cout << "ldrz: " << ldrz << std::endl;
        // std::cout << "ldout: " << ldout << std::endl;
        // // column-order matrix
        // std::cout << "rz+ldrz: " << *(rz+ldrz) << std::endl;
        // std::cout << "rz+ldrz+1: " << *(rz+ldrz+1) << std::endl;
        // std::cout << "out: " << *(out) << std::endl;
        // std::cout << "out+1: " << *(out+1) << std::endl;
        // std::cout <<"########### inside drop col:\t" << n << " " << std::endl;

        givens::zero(n, rz + ldrz, ldrz, rz + ldrz + 1, ldrz, out, ldout,
                     out + 1, ldout);

        // std::cout << n << " " << std::endl;
        //     for (int i = 0; i < n; ++i) {
        //         for (int j = 0; j < n; ++j ) {
        //             std::cout << *(out + j*n + i) << "\t" ;
        //         }
        //         std::cout << std::endl;
        //     }

        while (--n > 0)
        {
            // increasing element
            rz += ldrz + 1;
            out += ldout + 1;

            givens::zero(n, out, ldout, rz + ldrz + 1, ldrz, out, ldout,
                         out + 1, ldout);
            // std::cout << "inside drop col:\t" << n << " " << std::endl;
            // for (int i = 0; i < n; ++i) {
            //     for (int j = 0; j < n; ++j ) {
            //         std::cout << *(out + j*n + i) << "\t" ;
            //     }
            //     std::cout << std::endl;
            // }
        }


    }


    static void
    drop_column(
        int n,
        const Scalar* rz,
        const int ldrz,
        Scalar* out,
        const int ldout, 
        const Scalar* qt,
        Scalar* qt_out,
        const int ldqt,
        const int m
    ) noexcept
    {
        // std::cout << n << std::endl;
        givens::zero(n, rz + ldrz, ldrz, rz + ldrz + 1, ldrz, out, ldout,
                     out + 1, ldout, qt, qt + 1, qt_out, qt_out + 1, ldqt, m);

        while (--n > 0)
        {
            // increasing element
            rz += ldrz + 1;
            out += ldout + 1;

            qt++;
            qt_out++;

            // std::cout << n << std::endl;
            if(n == 1) {
                givens::zero(n, out, ldout, rz + ldrz + 1, ldrz, out, ldout,
                         out + 1, ldout);
            } else {
                givens::zero(n, out, ldout, rz + ldrz + 1, ldrz, out, ldout,
                         out + 1, ldout, qt_out, qt + 1, qt_out, qt_out + 1, ldqt, m);
            }
        }

    }



    static void
    column_bounds(
        const int n,
        const Scalar* const rz,
        int ldrz,
        Scalar* out,
        givens* const aux_givens
    ) noexcept
    {
        const Scalar* colj = rz;

        for (int j = 0; j < n; ++j, colj += ldrz)
        {
            const Scalar* coli = colj + ldrz;

            for (int i = j + 1; i <= n; ++i, coli += ldrz)
            {
                Scalar t = coli[j];

                for (int g = j + 1; g < i; ++g)
                {
                    t = -aux_givens[g].s() * t + aux_givens[g].c() * coli[g];
                }

                aux_givens[i].gen(t, coli[i]);
            }

            *(out++) = std::abs(aux_givens[n].r());
        }
    }



    static void
    permute_complete(
        const int n,
        const Scalar* rz,
        const int ldrz,
        const int* const pos,
        Scalar* const out,
        const int ldout,
        Scalar* const aux_tau,
        Scalar* const aux_work
    ) noexcept
    {
        for (int i = 0; i < n; ++i)
        {
            const int j = pos[i];

            std::copy_n(rz + j * ldrz, j + 1, out + i * ldout);
            std::fill_n(out + i * ldout + j + 1, n - j, 0);
        }

        std::copy_n(rz + n * ldrz, n + 1, out + n * ldout);
        lapack::geqr2(n + 1, n + 1, out, ldout, aux_tau, aux_work);
    }



    static void
    permute_partial_1(
        const int n,
        Scalar* const rz,
        const int ldrz,
        const int pos
    ) noexcept
    {
        const int j = pos;

        Scalar* const col0 = rz;
        Scalar* const colj = rz + j * ldrz;

        for (int i = j - 1; i >= 0; --i)
        {
            const givens g(colj[i], colj[i + 1]);
            g.rot(j - i, colj - ldrz + i, -ldrz, colj - ldrz + i + 1, -ldrz);
            g.rot(n - j, colj + ldrz + i, ldrz, colj + ldrz + i + 1, ldrz);

            colj[i] = g.r();
            colj[i + 1] = 0;
        }

        std::swap(col0[0], colj[0]);
        std::swap(col0[1], colj[1]);

        Scalar* coli = col0;

        for (int i = 1; i < j; ++i)
        {
            coli += ldrz;

            const givens g(coli[i], coli[i + 1]);
            g.rot(n - i, coli + ldrz + i, ldrz, coli + ldrz + i + 1, ldrz);

            coli[i] = g.r();
            coli[i + 1] = 0;
        }
    }



    static void
    permute_partial_2(
        const int n,
        Scalar* const rz,
        const int ldrz,
        const int pos
    ) noexcept
    {
        if (pos < 1)
        {
            return;
        }

        const int j = pos;

        Scalar* const col0 = rz;
        Scalar* const colj = rz + j * ldrz;

        *(colj - ldrz + j) = 0;

        givens g(colj[j - 1], colj[j]);
        g.rot(1, colj - ldrz + j - 1, -ldrz, colj - ldrz + j, -ldrz,
              colj + j - 1, -ldrz, colj + j, -ldrz);
        g.rot(n - j, colj + ldrz + j - 1, ldrz, colj + ldrz + j, ldrz);

        Scalar* rii = colj - ldrz + j - 1;
        for (int i = j - 1; i > 0; --i, rii -= ldrz + 1)
        {
            *rii = 0;

            g.gen(colj[i - 1], g.r());
            g.rot(j - i + 1, colj - ldrz + i - 1, -ldrz, colj + i, -ldrz,
                  colj + i - 1, -ldrz, colj + i, -ldrz);
            g.rot(n - j, colj + ldrz + i - 1, ldrz, colj + ldrz + i, ldrz);
        }

        col0[0] = g.r();
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
