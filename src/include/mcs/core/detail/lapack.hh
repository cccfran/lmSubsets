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



#ifndef MCS_CORE_DETAIL_LAPACK_HH
#define MCS_CORE_DETAIL_LAPACK_HH



#include <string>



#include "gsl/gsl"  // gsl::span



#include "mcs/core/detail/matrix.hh"
#include "mcs/core/detail/vector.hh"



namespace mcs    {
namespace core   {
namespace detail {



struct lapack_base
{

    static const std::string upper;

    static const std::string lower;

    static const std::string full;

    static const std::string trans;

    static const std::string no_trans;

    static const std::string left;

    static const std::string right;

};



template<typename Scalar>
class lapack : private lapack_base
{

    using vector_span = vector<Scalar&>;

    using vector_cspan = vector<const Scalar&>;

    using matrix_span = matrix<Scalar&>;

    using matrix_cspan = matrix<const Scalar&>;



public:

    using lapack_base::upper;

    using lapack_base::lower;

    using lapack_base::full;

    using lapack_base::trans;

    using lapack_base::no_trans;

    using lapack_base::left;

    using lapack_base::right;



public:

    static void
    lacpy(
        matrix_cspan a,
        matrix_span b
    ) noexcept
    {
        lacpy(full, a, b);
    }



    static void
    lacpy(
        const std::string& uplo,
        matrix_cspan a,
        matrix_span b
    ) noexcept
    {
        const int m = a.nrow();
        const int n = a.ncol();

        lacpy(uplo.c_str(), m, n, a.base(), a.ldim(), b.base(), b.ldim());
    }



    static void
    lacpy(
        const char* uplo,
        int m,
        int n,
        const Scalar* a,
        int lda,
        Scalar* b,
        int ldb
    ) noexcept;



    static int
    geqrf(
        matrix_span a,
        gsl::span<Scalar> tau,
        gsl::span<Scalar> work
    ) noexcept
    {
        const int m = a.nrow();
        const int n = a.ncol();

        return geqrf(m, n, a.base(), a.ldim(), tau.data(), work.data(),
                     work.size());
    }



    static int
    geqrf(
        int m,
        int n,
        Scalar* a,
        int lda,
        Scalar* tau,
        Scalar* work,
        int lwork
    ) noexcept;



    static int
    geqr2(
        matrix_span a,
        gsl::span<Scalar> tau,
        gsl::span<Scalar> work
    ) noexcept
    {
        const int m = a.nrow();
        const int n = a.ncol();

        return geqr2(m, n, a.base(), a.ldim(), tau.data(), work.data());
    }



    static int
    geqr2(
        int m,
        int n,
        Scalar* a,
        int lda,
        Scalar* tau,
        Scalar* work
    ) noexcept;



    static int
    orgqr(
        const int k,
        matrix_span a,
        vector_cspan tau,
        gsl::span<Scalar> work
    ) noexcept
    {
        const int m = a.nrow();
        const int n = a.ncol();

        return orgqr(m, n, k, a.base(), a.ldim(), tau.base(), work.data(),
                     work.size());
    }



    static int
    orgqr(
        int m,
        int n,
        int k,
        Scalar* a,
        int lda,
        const Scalar* tau,
        Scalar* work,
        int lwork
    ) noexcept;



    // static int
    // ormqr(
    //     const std::string& side,
    //     const std::string& trans,
    //     const int k,
    //     matrix_cspan a,
    //     vector_cspan tau,
    //     matrix_span c,
    //     gsl::span<Scalar> work
    // ) noexcept
    // {
    //     const int m = c.nrow();
    //     const int n = c.ncol();

    //     return ormqr(side.c_str(), trans.c_str(), m, n, k, a.base(), a.ldim(),
    //                  tau.base(), c.base(), c.ldim(), work.data(), work.size());
    // }

    static int
    ormqr(
        const std::string& side,
        const std::string& trans,
        const int k,
        matrix_cspan a,
        std::vector<Scalar> tau,
        matrix_span c,
        gsl::span<Scalar> work
    ) noexcept
    {
        const int m = c.nrow();
        const int n = c.ncol();

        return ormqr(side.c_str(), trans.c_str(), m, n, k, a.base(), a.ldim(),
                     tau.data(), c.base(), c.ldim(), work.data(), work.size());
    }


    static int
    ormqr(
        const char* side,
        const char* trans,
        int m,
        int n,
        int k,
        const Scalar* a,
        int lda,
        const Scalar* tau,
        Scalar* c,
        int ldc,
        Scalar* work,
        int lwork
    ) noexcept;


    static int
    trtrs(
        matrix_cspan a,
        matrix_span b
    ) noexcept
    {
        const int n = a.nrow();
        const int d = a.ncol();
        const int dy = b.ncol();
        const char uplo = 'U';
        const char trans = 'N';
        const char diag = 'N';
        
        return trtrs(&uplo, &trans, &diag, d, dy, a.base(), a.ldim(), b.base(), b.ldim());
    }


    static int
    trtrs(
        const char* uplo,
        const char* trans,
        const char* diag,
        int n,
        int nrhs,
        const Scalar* a,
        int lda,
        Scalar* b,
        int ldb
    ) noexcept;


    static int
    gemm(
        const char* transA,
        const char* transB,
        int m,
        int n,
        int k,
        double alpha,
        const Scalar* A,
        int lda,
        const Scalar* B,
        int ldb,
        double beta, 
        Scalar* C, 
        int ldc
    ) noexcept;
    // C = alpha * A * B + beta*C
    // m: row of A,C; n: col of B; k: col of A

public:

    lapack() = delete;

};



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs



#include "mcs/core/detail/lapack.inc"
#endif
