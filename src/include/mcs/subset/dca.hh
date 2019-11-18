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



#ifndef MCS_SUBSET_DCA_HH
#define MCS_SUBSET_DCA_HH



#include <utility>  // std::pair
#include <vector>



#include "mcs/core/matrix.hh"

#include "mcs/subset/table.hh"

#include "mcs/subset/detail/dca.hh"
#include "mcs/subset/detail/dca_preo.hh"
#include "mcs/subset/detail/dca_state.hh"



namespace mcs    {
namespace subset {



template<typename Scalar,
         typename CostFunc>
std::pair<table_best<Scalar>, int>
dca_best(
    mcs::core::matrix<const Scalar&> ay_mat,
    const int mark,
    const CostFunc& cost_func,
    const int nbest,
    const int prad = 0
) noexcept
{
    using namespace mcs::subset::detail;

    using preo_complete = dca_preo::complete<Scalar>;
    using preo_radius = dca_preo::radius<Scalar, preo_complete>;
    using dca_state = dca_state_best<Scalar, CostFunc, preo_radius>;

    dca_state state(ay_mat, mark, cost_func, nbest, preo_radius(prad));
    const int nodes = detail::dca_best<Scalar, dca_state>(state);

    return { state.table(), nodes };
}



template<typename Scalar>
std::pair<table_all<Scalar>, int>
dca_all(
    mcs::core::matrix<const Scalar&> ay_mat,
    const int mark,
    const int nbest,
    const int prad = 0
) noexcept
{
    using namespace mcs::subset::detail;

    using preo_complete = dca_preo::complete<Scalar>;
    using preo_radius = dca_preo::radius<Scalar, preo_complete>;
    using dca_state = dca_state_all<Scalar, preo_radius>;

    // std::cout << "Inside dca" << std::endl;
    dca_state state(ay_mat, mark, nbest, preo_radius(prad));
    // std::cout << "FInish init" << std::endl;
    const int nodes = detail::dca_all<Scalar, dca_state>(state);

    // build table and return number of nodes
    return { state.table(), nodes };
}


template<typename Scalar>
std::pair<table_all<Scalar>, int>
dca_all_boot(
    mcs::core::matrix<const Scalar&> ay_mat,
    const int mark,
    const int nbest,
    const int prad = 0,
    const int nboot = 100,
    const int boot_type = 1
) noexcept
{
    using namespace mcs::subset::detail;

    using preo_complete = dca_preo::complete<Scalar>;
    using preo_radius = dca_preo::radius<Scalar, preo_complete>;
    using dca_state = dca_state_all_boot<Scalar, preo_radius>;

    // std::cout << "Inside dca" << std::endl;
    dca_state state(ay_mat, mark, nbest, preo_radius(prad), nboot, boot_type);
    // std::cout << "FInish init" << std::endl;
    const int nodes = detail::dca_all<Scalar, dca_state>(state);

    // build table and return number of nodes
    return { state.table(), nodes };
}



}  // namespace subset
}  // namespace mcs



#endif
