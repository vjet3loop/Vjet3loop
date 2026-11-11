/******************************************************************************
 *                                                                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 * Copyright (C) 2025 Michal Czakon                                           *
 *                                                                            *
 ******************************************************************************/

#ifndef DIFFEQS_ERROR_CHECKER_HPP
#define DIFFEQS_ERROR_CHECKER_HPP

#include "types.hpp"

template <class Real>
struct diffeqs_error_checker
{
  using types = diffeqs_types<Real>;

  using real = typename types::real;
  using cplx_vec = typename types::cplx_vec;
  using error_estimator_t = typename types::error_estimator_t;

  real error(const cplx_vec& f, const cplx_vec& f_err) const
  {return _error_estimator(f,f_err)/_eps;};

  error_estimator_t _error_estimator;
  real              _eps;
};

#endif
