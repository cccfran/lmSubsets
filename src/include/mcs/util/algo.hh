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



#ifndef MCS_UTIL_ALGO_HH
#define MCS_UTIL_ALGO_HH



#include "mcs/util/detail/algo.hh"



namespace mcs  {
namespace util {



using detail::arrange;
using detail::arrange_n;
using detail::for_each;
using detail::iota;
using detail::repeat;
using detail::reverse;
using detail::concat;
using detail::transform;
using detail::sort_heap;
// using detail::map;
// using detail::plus;



}  // end namespace util
}  // end namespace mcs



#endif
