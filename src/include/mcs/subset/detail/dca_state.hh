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



#ifndef MCS_SUBSET_DETAIL_DCA_STATE_HH
#define MCS_SUBSET_DETAIL_DCA_STATE_HH



#include <type_traits>  // std::is_same
#include <utility>  // std::declval
#include <vector>

#include <random>


#include "mcs/core/matrix.hh"

#include "mcs/subset/detail/dca_node.hh"
#include "mcs/subset/detail/dca_partial.hh"
#include "mcs/subset/detail/dca_qrz.hh"

#include "mcs/util/algo.hh"  // algo::concat, algo::iota, algo::map, algo::plus,
                             // algo::repeat, algo::transform
#include "mcs/util/function_traits.hh"

#include "mcs/core/ziggurat.hpp"  // Ziggurat normal generator



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar,
         typename NodeXfer>
class dca_state_base
{

    using dca_node = detail::dca_node<Scalar>;

    using dca_qrz = detail::dca_qrz<Scalar>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;

    // return type of NodeXfer.make(dca_state_base)
    // where NodeXfer is of type dca_state
    // and make return an instance with R and root size
    using node_xfer_inst = decltype(
        std::declval<NodeXfer>().
        make(std::declval<dca_state_base>())
    );



public:

    // node stack 
    std::vector<dca_node> node_stk_;

    typename decltype(node_stk_)::iterator cur_node_;

    typename decltype(node_stk_)::iterator nxt_node_;

    node_xfer_inst node_xfer_;

    dca_qrz qrz_;

    int root_size_;

    int root_mark_;

    int root_rank_;

    Scalar root_rss_;

    matrix_cspan X_;

    matrix_cspan y_;

public:

    dca_state_base(
        matrix_cspan ay_mat,
        const int mark,
        const NodeXfer& node_xfer
    ) noexcept :
        qrz_(ay_mat.ncol() - 1, ay_mat.nrow()),
        root_size_(ay_mat.ncol() - 1),
        root_mark_(mark),
        root_rank_(root_size_ - root_mark_),
        X_(ay_mat({0, ay_mat.nrow()}, {0, ay_mat.ncol()-1})),
        y_(ay_mat({0, ay_mat.nrow()}, {ay_mat.ncol()-1, 1}))
    {
        const int n = root_size_;
        const int k = root_mark_;
        const int p = root_rank_;

        // std::cout << "root size: " << n << std::endl;
        // std::cout << "root mark: " << k << std::endl;
        // std::cout << "root rank: " << p << std::endl;

        // size of the stack
        node_stk_.reserve(p);
        for (int i = 0; i < p; ++i)
        {   
            // emplace_back: append new container (dca_node) at the end
            node_stk_.emplace_back(p, ay_mat.nrow(), 
                &qrz_, X_, y_);
            // node_stk_.emplace_back(p, ay_mat.nrow(), &qrz_);
        }

        cur_node_ = node_stk_.begin();

        nxt_node_ = cur_node_ + 1;
        // pass the tri
        // p+1 because of y
        nxt_node_->root(qrz_.rz(ay_mat)({k, p+1}, {k, p+1}));
        // BOOTSTRAP
        nxt_node_->get_t(root_mark_, qrz_, X_, y_);

        root_rss_ = nxt_node_->rss();

        
        // a return rank_instance<Scalar, complete_ins, null_inst>
        node_xfer_ = node_xfer.make(*this);


    }


    dca_state_base(
        matrix_cspan ay_mat,
        const int mark,
        const NodeXfer& node_xfer,
        const int nboot
    ) noexcept :
        qrz_(ay_mat.ncol() - 1, ay_mat.nrow()),
        root_size_(ay_mat.ncol() - 1),
        root_mark_(mark),
        root_rank_(root_size_ - root_mark_),
        X_(ay_mat({0, ay_mat.nrow()}, {0, ay_mat.ncol()-1})),
        y_(ay_mat({0, ay_mat.nrow()}, {ay_mat.ncol()-1, 1}))
    {
        const int n = root_size_;
        const int k = root_mark_;
        const int p = root_rank_;

        // std::cout << "root size: " << n << std::endl;
        // std::cout << "root mark: " << k << std::endl;
        // std::cout << "root rank: " << p << std::endl;

        // size of the stack
        node_stk_.reserve(p);
        for (int i = 0; i < p; ++i)
        {   
            // emplace_back: append new container (dca_node) at the end
            node_stk_.emplace_back(p, ay_mat.nrow(), 
                &qrz_, X_, y_, nboot);
            // node_stk_.emplace_back(p, ay_mat.nrow(), &qrz_);
        }

        cur_node_ = node_stk_.begin();

        nxt_node_ = cur_node_ + 1;
        // pass the tri
        // p+1 because of y
        nxt_node_->root(qrz_.rz(ay_mat)({k, p+1}, {k, p+1}));
        // BOOTSTRAP
        // nxt_node_->get_t(root_mark_, qrz_, X_, y_, nboot);

        root_rss_ = nxt_node_->rss();

        
        // a return rank_instance<Scalar, complete_ins, null_inst>
        node_xfer_ = node_xfer.make(*this);


    }



public:

    bool
    is_final() const noexcept
    {
        return cur_node_ == nxt_node_;
    }



    void
    next_node() noexcept
    {
        // preorder: rank_inst: L116
        node_xfer_(*nxt_node_, *cur_node_);
        --nxt_node_;
    }



    int
    node_size() const  noexcept
    {
        return root_mark_ + cur_node_->size();
    }



    int
    node_mark() const noexcept
    {
        return root_mark_ + cur_node_->mark();
    }



    Scalar
    node_rss() const noexcept
    {
        return cur_node_->rss();
    }



    void
    drop_column(const int mark) noexcept
    {
        ++nxt_node_;
        cur_node_->drop_column(mark - root_mark_, *nxt_node_, qrz_);

        // std::cout << "\t" << "mark" << mark-root_mark_ << std::endl;
        // std::cout << "\t";
        // for(auto it=nxt_node_->subset().cbegin();
        //     it != nxt_node_->subset().cend(); ++it) {
        //         std::cout << ' ' << *it;
        //     }
        // std::cout << std::endl;
    }

    void
    get_t(const int mark) noexcept
    {
        // BOOTSTRAP
        nxt_node_->get_t(mark - root_mark_, qrz_, X_, y_);
    }

    int
    root_size() const noexcept
    {
        return root_size_;
    }



    int
    root_mark() const noexcept
    {
        return root_mark_;
    }



    int
    root_rank() const noexcept
    {
        return root_rank_;
    }



    Scalar
    root_rss() const noexcept
    {
        return root_rss_;
    }



    const dca_qrz&
    qrz() const noexcept
    {
        return qrz_;
    }

    typename decltype(node_stk_)::iterator get_cur_node() 
    {
        return cur_node_;
    }
};



template<typename Scalar,
         typename NodeXfer>
class dca_state_all : private dca_state_base<Scalar, NodeXfer>
{

    using base = dca_state_base<Scalar, NodeXfer>;

    using dca_partial = detail::dca_partial_all<Scalar>;

    using dca_result = detail::dca_result<Scalar>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;



private:

    dca_partial partial_;

    int nbest_;



public:

    dca_state_all(
        matrix_cspan ay_mat,
        const int mark,
        const int nbest,
        const NodeXfer& node_xfer
    ) noexcept :
        base(ay_mat, mark, node_xfer),
        partial_(base::root_rank(), nbest),
        nbest_(nbest)
    {
    }


public:

    using base::is_final;

    using base::node_size;

    using base::node_mark;

    using base::node_rss;

    using base::drop_column;

    using base::root_rss;

    using base::get_cur_node;

    using base::get_t;



    Scalar
    rss_inf() const noexcept
    {
        return root_rss();
    }



    void
    next_node() noexcept
    {
        base::next_node();
        // partial get update!
        partial_.update(*base::cur_node_);
        // partial_.update(*base::cur_node_, base::qrz_);
    }



    Scalar
    rss_bound() const noexcept
    {
        return base::cur_node_->rss();
    }



    Scalar
    min_rss(const int size) const noexcept
    {
        return partial_.min_rss(size - base::root_mark());
    }



    std::vector<std::vector<dca_result>>
    table() const noexcept
    {
        const int root_mark = this->root_mark();

        const auto prefix = util::iota(0, root_mark);

        // return dca_result type
        const auto xform = [&prefix, &root_mark](
            const dca_result& r
        ) -> dca_result {
            if (!r)  return {};

            // auto subset_ret = util::transform(r, [&root_mark](int i) {
            //                 return i + root_mark;
            //             });

            // std::cout << "*** table: " << r.key() << std::endl;
            // std::cout << "prefix: " << std::endl;
            // for(auto x : subset_ret) std::cout << x << ' ';
            // std::cout << std::endl;
            // for(auto x : r.subset()) std::cout << x  << ' ';
            // std::cout << std::endl;

            return {
                util::concat(
                    prefix,
                    // [] returns the subset[i] of this state
                    // for each subset[i], add root_mark
                    util::transform(r, [&root_mark](int i) {
                            return i + root_mark;
                        })
                ),
                r.key()
            };
        };
        

        return util::concat(
            util::repeat(util::repeat(dca_result(), nbest_), root_mark),
            util::transform(
                // partial_ restore results of rss in a heap
                partial_.results(),
                [&xform](const std::vector<dca_result>& r) {
                    return util::transform(r, xform);
                }
            )
        );
    }

};




template<typename Scalar,
         typename NodeXfer>
class dca_state_all_boot : private dca_state_base<Scalar, NodeXfer>
{

    using base = dca_state_base<Scalar, NodeXfer>;

    using dca_partial = detail::dca_partial_all_boot<Scalar>;

    using dca_result = detail::dca_result<Scalar>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;



private:

    dca_partial partial_;

    int nbest_;

    // BOOTSTRAP
    int nboot_;

public: 

    std::vector<double> multiplier_mat;

public:


    dca_state_all_boot(
        matrix_cspan ay_mat,
        const int mark,
        const int nbest,
        const NodeXfer& node_xfer,
        // BOOTSTRAP
        const int nboot
    ) noexcept :
        // BOOTSTRAP
        base(ay_mat, mark, node_xfer, nboot),
        partial_(base::root_rank(), nbest, nboot),
        nbest_(nbest),
        nboot_(nboot)
    {   

        multiplier_mat.reserve(nboot * ay_mat.nrow());
        std::mt19937_64 random;
        cxx::ziggurat_normal_distribution<double> normal{0.0, 1.0};

        for(int i = 0; i < nboot * ay_mat.nrow(); ++i) {
            multiplier_mat.push_back(normal(random));
        }

        base::nxt_node_->get_t(base::root_mark_, base::qrz_, base::X_, base::y_, nboot, multiplier_mat);
    }




public:

    using base::is_final;

    using base::node_size;

    using base::node_mark;

    using base::node_rss;

    using base::drop_column;

    using base::root_rss;

    using base::get_cur_node;

    void
    get_t(const int mark) noexcept
    {
        // BOOTSTRAP
        base::nxt_node_->get_t(mark - base::root_mark_, base::qrz_, base::X_, base::y_, nboot_, multiplier_mat);
    }



    Scalar
    rss_inf() const noexcept
    {
        return root_rss();
    }



    void
    next_node() noexcept
    {
        base::next_node();
        // partial get update!
        partial_.update(*base::cur_node_);
        // partial_.update(*base::cur_node_, base::qrz_);
    }



    Scalar
    rss_bound() const noexcept
    {
        return base::cur_node_->rss();
    }



    Scalar
    min_rss(const int size) const noexcept
    {
        return partial_.min_rss(size - base::root_mark());
    }



    std::vector<std::vector<dca_result>>
    table() const noexcept
    {
        const int root_mark = this->root_mark();

        const auto prefix = util::iota(0, root_mark);

        // return dca_result type
        const auto xform = [&prefix, &root_mark](
            const dca_result& r 
        ) -> dca_result {
            if (!r)  return {};

            // auto subset_ret = util::transform(r, [&root_mark](int i) {
            //                 return i + root_mark;
            //             });

            // std::cout << "*** table: " << r.key() << std::endl;
            // std::cout << "prefix: " << std::endl;
            // for(auto x : subset_ret) std::cout << x << ' ';
            // std::cout << std::endl;
            // for(auto x : r.subset()) std::cout << x  << ' ';
            // std::cout << std::endl;

            return {
                util::concat(
                    prefix,
                    // [] returns the subset[i] of this state
                    // for each subset[i], add root_mark
                    util::transform(r, [&root_mark](int i) {
                            return i + root_mark;
                        })
                ),
                r.key()
            };
        };
        

        return util::concat(
            util::repeat(util::repeat(dca_result(), nbest_), root_mark),
            util::transform(
                // partial_ restore results of rss in a heap
                partial_.results(),
                [&xform](const std::vector<dca_result>& r) {
                    return util::transform(r, xform);
                }
            )
        );
    }

};



template<typename Scalar,
         typename CostFunc,
         typename NodeXfer>
class dca_state_best : private dca_state_base<Scalar, NodeXfer>
{

    using cost_func_traits = mcs::util::function_traits<CostFunc>;

    static_assert(
        std::is_same<
            typename cost_func_traits::signature,
            double(int,double)
        >::value,
        "cost function must be 'double(int,double)'"
    );



    using base = dca_state_base<Scalar, NodeXfer>;

    using dca_partial = detail::dca_partial_best<Scalar>;

    using dca_result = detail::dca_result<Scalar>;

    using matrix_cspan = mcs::core::matrix<const Scalar&>;



private:

    dca_partial partial_;

    CostFunc cost_func_;

    Scalar cost_inf_;



public:

    dca_state_best(
        matrix_cspan ay_mat,
        const int mark,
        const CostFunc& cost_func,
        const int nbest,
        const NodeXfer& node_xfer,
        // BOOTSTRAP
        const int nboot
    ) noexcept :
        base(ay_mat, mark, node_xfer, nboot),
        partial_(base::root_rank(), nbest),
        cost_func_(cost_func)
    {
        cost_inf_ = cost_func_(mark + 1, base::root_rss());
    }



public:

    using base::is_final;

    using base::node_size;

    using base::node_mark;

    using base::node_rss;

    using base::drop_column;

    using base::get_t;



    Scalar
    cost_inf() const noexcept
    {
        return cost_inf_;
    }



    void
    next_node() noexcept
    {
        base::next_node();
        partial_.update(
            *base::cur_node_,
            [this](
                const int size,
                const Scalar rss
            ) -> Scalar {
                return cost_func_(base::root_mark() + size, rss);
            }
        );
    }



    Scalar
    cost_bound(const int mark) const noexcept
    {
        return cost_func_(mark + 1, base::cur_node_->rss());
    }



    Scalar
    min_cost() const noexcept
    {
        return partial_.min_cost();
    }



    std::vector<dca_result>
    table() const noexcept
    {
        const int root_mark = this->root_mark();

        const auto prefix = util::iota(0, root_mark);

        const auto xform = [&prefix, &root_mark](
            const dca_result& r
        ) -> dca_result {
            if (!r)  return {};

            return {
                util::concat(
                    prefix,
                    util::transform(r, [&root_mark](int i) {
                            return i + root_mark;
                        })
                ),
                r.key()
            };
        };

        return util::transform(partial_.results(), xform);
    }

};



}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
