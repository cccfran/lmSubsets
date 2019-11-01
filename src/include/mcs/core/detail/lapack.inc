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
#error Do not include this file!
#endif



#include "mcs/detail/fortran.hh"  // MCS_F77_NAME, MCS_F77_CALL



namespace mcs    {
namespace core   {
namespace detail {



extern "C" {



void
MCS_F77_NAME(dlacpy)(
    const char* uplo,
    const int* m,
    const int* n,
    const double* a,
    const int* lda,
    double* b,
    const int* ldb
);



void
MCS_F77_NAME(dgeqrf)(
    const int* m,
    const int* n,
    double* a,
    const int* lda,
    double* tau,
    double* work,
    const int* lwork,
    int* info
);



void
MCS_F77_NAME(dgeqr2)(
    const int* m,
    const int* n,
    double* a,
    const int* lda,
    double* tau,
    double* work,
    int* info
);



void
MCS_F77_NAME(dorgqr)(
    const int* m,
    const int* n,
    const int* k,
    double* a,
    const int* lda,
    const double* tau,
    double* work,
    const int* lwork,
    int* info
);



void
MCS_F77_NAME(dormqr)(
    const char* side,
    const char* trans,
    const int* m,
    const int* n,
    const int* k,
    const double* a,
    const int* lda,
    const double* tau,
    double* c,
    const int* ldc,
    double* work,
    const int* lwork,
    int* info
);


void
MCS_F77_NAME(dtrtrs)(
    const char* uplo,
    const char* trans,
    const char* diag,
    const int* n,
    const int* nrhs,
    const double* a,
    const int* lda,
    double* b,
    const int* ldb,
    int* info
);



void
MCS_F77_NAME(dgemm)(
    const char* transA,
    const char* transB,
    const int* m,
    const int* n,
    const int* k,
    const double* alpha,
    const double* A,
    const int* lda,
    const double* B,
    const int* ldb,
    const double* beta, 
    double* C, 
    const int* ldc
);


double
MCS_F77_NAME(dnrm2)(
    const int* n,
    const double* X,
    const int* incx
);


}  // end extern "C"



const std::string lapack_base::upper = "U";

const std::string lapack_base::lower = "L";

const std::string lapack_base::full = "_";

const std::string lapack_base::trans = "T";

const std::string lapack_base::no_trans = "N";

const std::string lapack_base::left = "L";

const std::string lapack_base::right = "R";



template<>
inline void
lapack<double>::lacpy(
    const char* const uplo,
    const int m,
    const int n,
    const double* const a,
    const int lda,
    double* const b,
    const int ldb
) noexcept
{
    MCS_F77_CALL(dlacpy)(uplo, &m, &n, a, &lda, b, &ldb);
}



template<>
inline int
lapack<double>::geqrf(
    const int m,
    const int n,
    double* const a,
    const int lda,
    double* const tau,
    double* const work,
    const int lwork
) noexcept
{
    int info = 0;

    MCS_F77_CALL(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);

    return info;
}



template<>
inline int
lapack<double>::geqr2(
    const int m,
    const int n,
    double* const a,
    const int lda,
    double* const tau,
    double* const work
) noexcept
{
    int info = 0;

    MCS_F77_CALL(dgeqr2)(&m, &n, a, &lda, tau, work, &info);

    return info;
}



template<>
inline int
lapack<double>::orgqr(
    const int m,
    const int n,
    const int k,
    double* const a,
    const int lda,
    const double* const tau,
    double* const work,
    const int lwork
) noexcept
{
    int info = 0;

    MCS_F77_CALL(dorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);

    return info;
}



template<>
inline int
lapack<double>::ormqr(
    const char* const side,
    const char* const trans,
    const int m,
    const int n,
    const int k,
    const double* const a,
    const int lda,
    const double* const tau,
    double* const c,
    const int ldc,
    double* const work,
    const int lwork
) noexcept
{
    int info = 0;

    MCS_F77_CALL(dormqr)(side, trans, &m, &n, &k, a, &lda, tau, c, &ldc, work,
                         &lwork, &info);

    return info;
}


template<>
inline int
lapack<double>::trtrs(
    const char* const uplo,
    const char* const trans,
    const char* const diag,
    const int n,
    const int nrhs,
    const double* const a,
    const int lda,
    double* const b,
    const int ldb
) noexcept
{
    int info = 0;

    MCS_F77_CALL(dtrtrs)(uplo, trans, diag, &n, &nrhs, a, &lda, b, &ldb, &info);

    return info;
}


template<>
inline int
lapack<double>::gemm(
    const char* const transA,
    const char* const transB,
    const int m,
    const int n,
    const int k,
    const double alpha,
    const double* const A,
    const int lda,
    const double* const B,
    const int ldb,
    double const beta, 
    double* const C, 
    const int ldc
) noexcept
{
    int info = 0;

    MCS_F77_CALL(dgemm)(transA, transB, &m, &n, &k, 
                        &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
    info = 1;

    return info;
}


template<>
inline double
lapack<double>::nrm2(
    const int n,
    const double* const X,
    const int incx
) noexcept
{
    double ret = MCS_F77_CALL(dnrm2)(&n, X, &incx);
    
    return ret;
}



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs
