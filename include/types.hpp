/******************************************************************************
 *                                                                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 * Copyright (C) 2026 Michal Czakon                                           *
 *                                                                            *
 ******************************************************************************/

#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>
#include <complex>
#include <functional>
#include <type_traits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace ublas = boost::numeric::ublas;

template <class Real>
struct diffeqs_types
{
  using real = Real;
  using cplx = std::complex<Real>;

  using real_vec = std::vector<Real>;

  using cplx_vec = ublas::vector<cplx>;
  using cplx_vec_range = ublas::vector_range<cplx_vec>;
  using matrix = ublas::compressed_matrix<cplx>;

  // Calculate the connection matrix A of the canonical equation given a point x
  // on the integration path its tangent dx and special functions f
  using connection_t = std::function<void(const cplx_vec& x,
                                          const cplx_vec& dx,
                                          const cplx_vec_range& f,
                                          matrix& A)>;

  // Calculate the right-hand side of the differential equations for the special
  // functions f given a point x on the integration path its tangent dx and
  // special functions f
  using vector_field_t = std::function<void(const cplx_vec& x,
                                            const cplx_vec& dx,
                                            const cplx_vec_range& f,
                                            cplx_vec_range& df)>;

  // Provide x(t) and dx(t) = x'(t) with x(0) = x0, x(1) = x1
  using path_t = std::function<void(const real_vec& x0,
                                    const real_vec& x1,
                                    real t,
                                    real_vec& x,
                                    real_vec& dx)>;

  // Caculate a total error estimate given estimates of the function values, f,
  // and their error estimates, f_err. This total error estimate is used by the
  // differential-equation solver that tries to achieve a requested error.
  using error_estimator_t = std::function<real(const cplx_vec& f,
                                               const cplx_vec& f_err)>;
};

/******************************************************************************
 *                                                                            *
 * Operations between int and cplx                                            *
 *                                                                            *
 ******************************************************************************/

namespace diffeqs_int_cplx
{

  template <class Int, class Real,
	    typename std::enable_if<std::is_integral<Int>::value,int>::type = 0>
  std::complex<Real> operator+(Int a, std::complex<Real> z)
  {
    z += Real(a);
    return z;
  }

  template <class Real, class Int,
	    typename std::enable_if<std::is_integral<Int>::value,int>::type = 0>
  std::complex<Real> operator+(std::complex<Real> z, Int a)
  {
    z += Real(a);
    return z;
  }

  template <class Int, class Real,
	    typename std::enable_if<std::is_integral<Int>::value,int>::type = 0>
  std::complex<Real> operator-(Int a, const std::complex<Real>& z)
  {
    std::complex<Real> r(Real(a), Real(0));
    r -= z;
    return r;
  }

  template <class Real, class Int,
	    typename std::enable_if<std::is_integral<Int>::value,int>::type = 0>
  std::complex<Real> operator-(std::complex<Real> z, Int a)
  {
    z -= Real(a);
    return z;
  }

  template <class Int, class Real,
	    typename std::enable_if<std::is_integral<Int>::value,int>::type = 0>
  std::complex<Real> operator*(Int a, std::complex<Real> z)
  {
    z *= Real(a);
    return z;
  }

  template <class Real, class Int,
	    typename std::enable_if<std::is_integral<Int>::value,int>::type = 0>
  std::complex<Real> operator*(std::complex<Real> z, Int a)
  {
    z *= Real(a);
    return z;
  }

  template <class Int, class Real,
	    typename std::enable_if<std::is_integral<Int>::value,int>::type = 0>
  std::complex<Real> operator/(Int a, const std::complex<Real>& z)
  {
    std::complex<Real> r(Real(a), Real(0));
    r /= z;
    return r;
  }

  template <class Real, class Int,
	    typename std::enable_if<std::is_integral<Int>::value,int>::type = 0>
  std::complex<Real> operator/(std::complex<Real> z, Int a)
  {
    z /= Real(a);
    return z;
  }

} // namespace diffeqs_int_cplx

#endif
