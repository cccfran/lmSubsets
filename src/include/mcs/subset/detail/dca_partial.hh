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



#ifndef MCS_SUBSET_DETAIL_DCA_PARTIAL_HH
#define MCS_SUBSET_DETAIL_DCA_PARTIAL_HH



#include <vector>



#include "gsl/gsl"  // gsl::span



#include "mcs/util/algo.hh"  // util::transform

#include "mcs/subset/detail/aux_heap.hh"
#include "mcs/subset/detail/dca_node.hh"
#include "mcs/subset/detail/dca_result.hh"



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar>
class dca_partial_all
{

    using aux_heap = detail::aux_heap<Scalar>;

    using dca_node = detail::dca_node<Scalar>;

    using dca_result = detail::dca_result<Scalar>;



private:

    std::vector<aux_heap> heaps_;



public:

    dca_partial_all(
        const int root_size,
        const int nbest
    ) noexcept
    {
        heaps_.reserve(root_size);
        std::cout << "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHeap" << root_size << std::endl;
        for (int size = 1; size <= root_size; ++size)
        {
            heaps_.emplace_back(size, nbest);
        }
    }



public:

    Scalar
    min_rss(const int size) const noexcept
    {
        return heaps_[size - 1].max_key();
    }



    void
    update(const dca_node& node) noexcept
    {
        node.for_each(
            [this](
                gsl::span<const int> subset,
                const Scalar rss
            ) -> void {
                const int size = subset.size();

                auto& heap = heaps_[size - 1];
                if (rss < heap.max_key())
                {
                    heap.insert(subset, rss);
                }
            });
    }


    // void
    // update(const dca_node& node, const dca_qrz& qrz_) noexcept
    // {
    //     node.for_each(
    //         [this](
    //             gsl::span<const int> subset,
    //             const Scalar rss
    //         ) -> void {
    //             const int size = subset.size();

    //             auto& heap = heaps_[size - 1];
    //             if (rss < heap.max_key())
    //             {
    //                 heap.insert(subset, rss);
    //             }
    //         });
    // }


    std::vector<std::vector<dca_result>>
    results() const noexcept
    {
        return util::transform(heaps_, [](const aux_heap& h) {
                return h.results();
            });
    }

};



template<typename Scalar>
class dca_partial_best
{

    using aux_heap = detail::aux_heap<Scalar>;

    using dca_node = detail::dca_node<Scalar>;

    using dca_result = detail::dca_result<Scalar>;



private:

    aux_heap heap_;



public:

    dca_partial_best(
        const int root_size,
        const int nbest
    ) noexcept :
        heap_(root_size, nbest)
    {
    }




public:

    Scalar
    min_cost() const noexcept
    {
        return heap_.max_key();
    }



    template<typename CostFunc>
    void
    update(
        const dca_node& node,
        const CostFunc& cost_func
    ) noexcept
    {
        node.for_each(
            [this, &cost_func](
                gsl::span<const int> subset,
                const Scalar rss
            ) {
                const int size = subset.size();
                const Scalar cost = cost_func(size, rss);

                if (cost < heap_.max_key())
                {
                    heap_.insert(subset, cost);
                }
            }
        );
    }



    std::vector<dca_result>
    results() const noexcept
    {
        return heap_.results();
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
