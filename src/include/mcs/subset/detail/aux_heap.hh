// Copyright 2018  Marc Hofmann
//
// This file is part of the 'mcs' library (see
// <https://github.com/marc-hofmann/mcs.cc/>).
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



#ifndef MCS_SUBSET_DETAIL_AUX_HEAP_HH
#define MCS_SUBSET_DETAIL_AUX_HEAP_HH



#include <algorithm>  // std::pop_heap, std::push_heap, std::sort_heap
#include <functional>  // std::function
#include <limits>  // std::numeric_limits
#include <vector>



#include "gsl/gsl"  // gsl::span



#include "mcs/util/algo.hh"  // util::sort_heap, util::transform

#include "mcs/subset/detail/dca_result.hh"



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class aux_heap
{

    using dca_result = detail::dca_result<Scalar>;



private:

    Scalar max_key_;

    Scalar min_key_;

    std::vector<int> heap_;

    std::function<bool(int,int)> heap_comp_;

    // keys is the RSS
    std::vector<Scalar> keys_;

    std::vector<std::vector<int>> subsets_;



public:

    aux_heap(
        const int root_size,
        const int nbest
    ) noexcept :
        max_key_(std::numeric_limits<Scalar>::min()),
        min_key_(max_key_),
        heap_comp_(
            [this](
                const int i,
                const int j
            ) -> bool {
                return keys_[i] > keys_[j];
            }
        )
    {
        heap_.reserve(nbest);
        keys_.reserve(nbest);
        subsets_.resize(nbest);

        for (int i = 0; i < nbest; ++i)
        {
            heap_.push_back(i);
            keys_.push_back(max_key_);
            subsets_[i].reserve(root_size);
        }
    }



public:

    Scalar
    min_key() const noexcept
    {
        return min_key_;
    }



    Scalar
    max_key() const noexcept
    {
        return max_key_;
    }



    void
    insert(
        gsl::span<const int> subset,
        const Scalar key
    ) noexcept
    {
        std::pop_heap(heap_.begin(), heap_.end(), heap_comp_);

        const int pos = heap_.back();
        keys_[pos] = key;
        subsets_[pos].assign(subset.begin(), subset.end());

        std::push_heap(heap_.begin(), heap_.end(), heap_comp_);

        if (key > max_key_)
        {
            max_key_ = key;
        }

        min_key_ = keys_[heap_.front()];
    }



    std::vector<dca_result>
    results() const noexcept
    {
        return util::transform(
            util::sort_heap(heap_, heap_comp_),
            [this](const int i) -> dca_result {
                return { subsets_[i], keys_[i] };
            }
        );
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
