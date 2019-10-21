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



#ifndef MCS_CORE_DETAIL_MODEL_HH
#define MCS_CORE_DETAIL_MODEL_HH



#include <cmath>  // std::abs, std::sqrt, std::copysign



#include "mcs/core/detail/matrix.hh"
// #include "mcs/core/detail/lapack.hh"



namespace mcs    {
namespace core   {
namespace detail {



template<typename Scalar>
class model
{

    using matrix = detail::matrix<Scalar>;

    using matrix_cspan = detail::matrix<Scalar&>;

    // using lapack = mcs::core::lapack<Scalar>;



public:

    matrix beta_;

    matrix sds_;

    matrix tstat_;

    int model_size;


public:

    model(int subset_size) noexcept :
        beta_(subset_size, 1),
        sds_(subset_size, 1),
        tstat_(subset_size, 1),
        model_size(subset_size)
    {
    }

public:

    const int
    size() const noexcept
    {
        return model_size;
    }
    
    const matrix&
    beta() const noexcept
    {
        return beta_;
    }

    matrix&
    beta() noexcept
    {
        return beta_;
    }

    const matrix&
    sds() const noexcept
    {
        return sds_;
    }

    const matrix&
    tstat() const noexcept
    {
        return tstat_;
    }

    void 
    set_beta(
        matrix beta
    ) noexcept
    {
        beta_ = beta;
    }

    void 
    set_sds(
        matrix& sds
    ) noexcept
    {
        sds_ = sds;
    }

    void 
    set_tstat(
        matrix& tstat
    ) noexcept
    {
        tstat_ = tstat;
    }

};



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs



#endif
