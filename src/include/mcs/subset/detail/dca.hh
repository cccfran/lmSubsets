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



#ifndef MCS_SUBSET_DETAIL_DCA_HH
#define MCS_SUBSET_DETAIL_DCA_HH



namespace mcs    {
namespace subset {
namespace detail {



template<typename Scalar,
         typename DcaState>
int
dca_impl(DcaState& state) noexcept
{
    int node_cnt = 0;

    while (!state.is_final())
    {
        // lots of things done
        // first preorder
        // then partial update
        
        state.next_node();

        const int n = state.node_size();
        const int k = state.node_mark();

        // std::cout << "***********************" << std::endl;
        // std::cout << "** node_cnt: " << node_cnt << std::endl;
        // std::cout << "node_size: " << n << std::endl;
        // std::cout << "marker: " << k << std::endl;
        // std::cout << "Subsets: ";
        // for(auto it=state.get_cur_node()->subset().cbegin();
        //     it != state.get_cur_node()->subset().cend(); ++it) {
        //         std::cout << *it << " ";
        //     }
        // std::cout << std::endl;

        // std::cout << "node: QRZ: " << std::endl;
        // for (int i = 0; i < state.get_cur_node()->rz_mat_.nrow(); ++i) {
        //     for (int j = 0; j < state.get_cur_node()->rz_mat_.ncol(); ++j) {
        //         std::cout << state.get_cur_node()->rz_mat_(i, j) << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // drop k+1 to n-1 columns to sub 
        for (int j = k; j < n - 1; ++j)
        {
            // std::cout << "######################" << std::endl;
            // std::cout << "dropping " << j << std::endl;

            state.drop_column(j);
            state.get_t(j);
        }

        // std::cout << std::endl;
        ++node_cnt;
    }

    // std::cout << "***********************" << std::endl;
    
    return node_cnt;
}



template<typename Scalar,
         typename DcaState>
int
dca_best(DcaState& state) noexcept
{
    return dca_impl<Scalar, DcaState>(state);
}



template<typename Scalar,
         typename DcaState>
int
dca_all(DcaState& state) noexcept
{
    return dca_impl<Scalar, DcaState>(state);
}



}  // namespace detail
}  // namespace subset
}  // namespace mcs



#endif
