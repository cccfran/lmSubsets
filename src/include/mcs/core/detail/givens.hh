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



#ifndef MCS_CORE_DETAIL_GIVENS_HH
#define MCS_CORE_DETAIL_GIVENS_HH



#include <cmath>  // std::abs, std::sqrt, std::copysign



#include "mcs/core/detail/vector.hh"



namespace mcs    {
namespace core   {
namespace detail {



template<typename Scalar>
class givens
{

    using vector_span = detail::vector<Scalar&>;

    using vector_cspan = detail::vector<const Scalar&>;



public:

    static void
    zero(
        vector_span x,
        vector_span y
    ) noexcept
    {
        const int n = x.len();

        zero(n, x.base(), x.inc(), y.base(), y.inc());
    }



    static void
    zero(
        const int n,
        Scalar* const x,
        const int incx,
        Scalar* const y,
        const int incy
    ) noexcept
    {
        if (n > 0)
        {
            const givens g(*x, *y);

            g.rot(n - 1, x + incx, incx, y + incx, incy);

            *x = g.r_;
            *y = 0;
        }
    }



    static void
    zero(
        vector_cspan x,
        vector_cspan y,
        vector_span xx,
        vector_span yy
    ) noexcept
    {
        const int n = x.len();

        zero(n, x.base(), x.inc(), y.base(), y.inc(), xx.base(), xx.inc(),
             yy.base(), yy.inc());
    }



    static void
    zero(
        const int n,
        const Scalar* const x,
        const int incx,
        const Scalar* const y,
        const int incy,
        Scalar* const xx,
        const int incxx,
        Scalar* const yy,
        const int incyy
    ) noexcept
    {
        // std::cout << "#### ZERO: " <<std::endl;
        if (n > 0)
        {
            const givens g(*x, *y);
            // update row
            g.rot(n - 1, x + incx, incx, y + incy, incy, xx + incxx, incxx,
                  yy + incyy, incyy);

            *xx = g.r_;
            *yy = 0;
        }
        // std::cout << std::endl;
    }



private:

    Scalar r_;

    Scalar c_;

    Scalar s_;



public:

    givens() noexcept :
        r_(0),
        c_(0),
        s_(0)
    {
    }



    givens(
        const Scalar dx,
        const Scalar dy
    ) noexcept
    {
        gen(dx, dy);
    }



public:

    Scalar
    r() const noexcept
    {
        return r_;
    }



    Scalar
    c() const noexcept
    {
        return c_;
    }



    Scalar
    s() const noexcept
    {
        return s_;
    }



    void
    gen(
        const Scalar dx,
        const Scalar dy
    ) noexcept
    {
        // see also:  https://en.wikipedia.org/wiki/Givens_rotation

        if (dy == 0)
        {
            c_ = std::copysign(1, dx);
            s_ = 0;
            r_ = std::abs(dx);
        }
        else if (dx == 0)
        {
            c_ = 0;
            s_ = std::copysign(1, dy);
            r_ = std::abs(dy);
        }
        else if (std::abs(dy) > std::abs(dx))
        {
            const Scalar t = dx / dy;
            const Scalar u = std::copysign(std::sqrt(1 + t * t), dy);

            s_ = 1 / u;
            c_ = s_ * t;
            r_ = dy * u;
        }
        else
        {
            const Scalar t = dy / dx;
            const Scalar u = std::copysign(std::sqrt(1 + t * t), dx);

            c_ = 1 / u;
            s_ = c_ * t;
            r_ = dx * u;
        }
    }



    void
    rot(
        vector_span x,
        vector_span y
    ) const noexcept
    {
        const int n = x.len();

        rot(n, x.base(), x.inc(), y.base(), y.inc());
    }



    void
    rot(
        const int n,
        Scalar* x,
        const int incx,
        Scalar* y,
        const int incy
    ) const noexcept
    {
        for (int i = 0; i < n; ++i)
        {
            const Scalar t = c_ * (*x) + s_ * (*y);
            *y = -s_ * (*x) + c_ * (*y);
            *x = t;

            x += incx;
            y += incy;
        }
    }



    void
    rot(
        vector_cspan x,
        vector_cspan y,
        vector_span xx,
        vector_span yy
    ) const noexcept
    {
        const int n = x.len();

        rot(n, x.base(), x.inc(), y.base(), y.inc(), xx.base(), xx.inc(),
            yy.base(), yy.inc());
    }



    void
    rot(
        const int n,
        const Scalar* x,
        const int incx,
        const Scalar* y,
        const int incy,
        Scalar* xx,
        const int incxx,
        Scalar* yy,
        const int incyy
    ) const noexcept
    {
        for (int i = 0; i < n; ++i)
        {   
            // std::cout << "###### Givens: " << std::endl;
            // std::cout << "###### (i,j)" << std::endl;
            const Scalar t = c_ * (*x) + s_ * (*y);
            *yy = -s_ * (*x) + c_ * (*y);
            *xx = t;

            x += incx;
            y += incy;
            xx += incxx;
            yy += incyy;
        }
    }

};



}  // end namespace detail
}  // end namespace core
}  // end namespace mcs



#endif
