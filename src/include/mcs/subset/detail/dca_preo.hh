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



#ifndef MCS_SUBSET_DETAIL_DCA_PREO_HH
#define MCS_SUBSET_DETAIL_DCA_PREO_HH



#include <vector>



#include "mcs/subset/detail/dca_node.hh"
#include "mcs/subset/detail/dca_qrz.hh"
#include "mcs/subset/detail/dca_state.hh"



namespace mcs      {
namespace subset   {
namespace detail   {
namespace dca_preo {



template<typename Scalar>
class null_inst
{

    using dca_node = detail::dca_node<Scalar>;



public:

    null_inst() noexcept
    {
    }



public:

    void
    operator ()(
        dca_node& node,
        dca_node& result
    ) const noexcept
    {
        std::cout << "\t\tSWAPPPING:" << std::endl;
        result.swap(node);
    }

};



template<typename Scalar,
         typename DcaPreoAInst,
         typename DcaPreoBInst>
class rank_inst
{

    using dca_node = detail::dca_node<Scalar>;



private:

    int lim_;

    DcaPreoAInst a_inst_;

    DcaPreoBInst b_inst_;



public:

    rank_inst() noexcept
    {
    }



    rank_inst(
        const int lim,
        const DcaPreoAInst& a_inst,
        const DcaPreoBInst& b_inst
    ) noexcept :
        lim_(lim),
        a_inst_(a_inst),
        b_inst_(b_inst)
    {
    }



public:

    void
    operator ()(
        dca_node& node,
        dca_node& result
    ) const noexcept
    {
        const int rank = node.rank();

        std::cout << "\tlimit: " << lim_ << std::endl;
        std::cout << "\trank: " << rank << std::endl;
        if (rank > lim_)
        {
            std::cout << "\t complete operator " << std::endl;
            a_inst_(node, result);
        }
        else
        {
            std::cout << "\t null operator " << std::endl;
            b_inst_(node, result);
        }
    }

};



template<typename Scalar>
class complete_inst
{

    using dca_node = detail::dca_node<Scalar>;

    using dca_qrz = detail::dca_qrz<Scalar>;



private:

    const dca_qrz* qrz_;

    int aux_size_;

    mutable std::vector<Scalar> aux_1_;

    mutable std::vector<int> aux_2_;



public:

    complete_inst() noexcept
    {
    }



    complete_inst(
        const dca_qrz* const qrz,
        const int aux_size
    ) noexcept :
        qrz_(qrz),
        aux_size_(aux_size),
        aux_1_(aux_size_),
        aux_2_(aux_size_)
    {
    }



public:

    void
    operator ()(
        const dca_node& node,
        dca_node& result
    ) const noexcept
    {
        std::cout << "\t\tComplelet inst operator" << std::endl;
        node.preorder_complete(result, *qrz_, aux_1_, aux_2_);
    }

};



template<typename Scalar>
class partial_1_inst
{

    using dca_node = detail::dca_node<Scalar>;

    using dca_qrz = detail::dca_qrz<Scalar>;



private:

    const dca_qrz* qrz_;

    int aux_size_;

    mutable std::vector<Scalar> aux_1_;



public:

    partial_1_inst() noexcept
    {
    }



    partial_1_inst(
        const dca_qrz* const qrz,
        const int aux_size
    ) noexcept :
        qrz_(qrz),
        aux_size_(aux_size),
        aux_1_(aux_size_)
    {
    }



public:

    void
    operator ()(
        dca_node& node,
        dca_node& result
    ) const noexcept
    {
        std::cout << "\t\t Partial1 inst operator" << std::endl;
        node.preorder_partial_1(result, *qrz_, aux_1_);
    }

};



template<typename Scalar>
class partial_2_inst
{

    using dca_node = detail::dca_node<Scalar>;

    using dca_qrz = detail::dca_qrz<Scalar>;



private:

    dca_qrz* const qrz_;

    int aux_size_;

    mutable std::vector<Scalar> aux_1_;



public:

    partial_2_inst() noexcept
    {
    }



    partial_2_inst(
        const dca_qrz* const qrz,
        const int aux_size
    ) noexcept :
        qrz_(qrz),
        aux_size_(aux_size),
        aux_1_(aux_size_)
    {
    }



public:

    void
    operator ()(
        dca_node& node,
        dca_node& result
    ) const noexcept
    {
        std::cout << "\t\t Partial2 inst operator" << std::endl;
        node.preorder_partial_2(result, *qrz_, aux_1_);
    }

};



template<typename Scalar>
class null
{

public:

    using instance = null_inst<Scalar>;



public:

    template<typename DcaState>
    instance
    make(const DcaState& state) const noexcept
    {
        std::cout << "\t\t000000000000000000000 NULL MAKE" << std::endl;
        return instance();
    }

};



template<typename Scalar,
         typename DcaPreoA,
         typename DcaPreoB = null<Scalar>>
class rank
{

public:

    using instance = rank_inst<
        Scalar,
        typename DcaPreoA::instance,
        typename DcaPreoB::instance>;



private:

    int lim_;

    DcaPreoA a_;

    DcaPreoB b_;



public:

    rank(
        const int lim,
        const DcaPreoA& a = DcaPreoA(),
        const DcaPreoB& b = DcaPreoB()
    ) noexcept :
        lim_(lim),
        a_(a),
        b_(b)
    {
    }



    template<typename DcaState>
    instance
    make(const DcaState& state) const noexcept
    {
        std::cout << "000000000000000000000rankkkkkk make" << std::endl;
        return instance(lim_, a_.make(state), b_.make(state));
    }

};



template<typename Scalar,
         typename DcaPreoA,
         typename DcaPreoB = null<Scalar>>
class radius
{

public:

    using instance = rank_inst<
        Scalar,
        typename DcaPreoA::instance,
        typename DcaPreoB::instance>;



private:

    int lim_;

    DcaPreoA a_;

    DcaPreoB b_;



public:

    radius(
        const int lim,
        // ctor of preo_complete
        const DcaPreoA& a = DcaPreoA(),
        const DcaPreoB& b = DcaPreoB()
    ) noexcept :
        lim_(lim),
        a_(a),
        b_(b)
    {
    }



    template<typename DcaState>
    instance
    make(const DcaState& state) const noexcept
    {
        // making radius
        const int lim = state.root_size() - state.root_mark() - lim_;
        std::cout << "\t\t RADIUS MAKEEEEEEE: " << lim << std::endl;

        // rank_inst oeprator
        // complete a_ and null b_
        return instance(lim, a_.make(state), b_.make(state));
    }

};



template<typename Scalar>
class complete
{

public:

    using instance = complete_inst<Scalar>;



public:

    template<typename DcaState>
    instance
    make(const DcaState& state) const noexcept
    {
        std::cout << "COMPELTE 00000000000000 " << std::endl;
        return instance(&state.qrz(), state.root_size());
    }

};



template<typename Scalar>
class partial_1
{

public:

    using instance = partial_1_inst<Scalar>;



public:

    template<typename DcaState>
    instance
    make(const DcaState& state) const noexcept
    {
        return instance(&state.qrz(), state.root_size());
    }

};



template<typename Scalar>
class partial_2
{

public:

    using instance = partial_2_inst<Scalar>;



public:

    template<typename DcaState>
    instance
    make(const DcaState& state) const noexcept
    {
        return instance(&state.qrz(), state.root_size());
    }

};



}  // end namespace dca_preo
}  // end namespace detail
}  // end namespace subset
}  // end namespace mcs



#endif
