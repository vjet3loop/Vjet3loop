/******************************************************************************
 *                                                                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 * Copyright (C) 2026 Michal Czakon                                           *
 *                                                                            *
 ******************************************************************************/

#ifndef DIFFEQS_HPP
#define DIFFEQS_HPP

#include <stdexcept>
#include <iostream>
#include <limits>

#include "types.hpp"

/******************************************************************************
 *                                                                            *
 * exceptions                                                                 *
 *                                                                            *
 ******************************************************************************/

struct diffeqs_exception : std::runtime_error
{
  using std::runtime_error::runtime_error;
};

struct diffeqs_max_steps_exception : diffeqs_exception
{
  diffeqs_max_steps_exception() :
    diffeqs_exception("Max number of steps exceeded.") {};
};

struct diffeqs_max_evals_exception : diffeqs_exception
{
  diffeqs_max_evals_exception() :
    diffeqs_exception("Max number of derivative evaluations exceeded.") {};
};

struct diffeqs_max_time_exception : diffeqs_exception
{
  diffeqs_max_time_exception() :
    diffeqs_exception("Max allowed evaluation time exceeded.") {};
};

struct diffeqs_min_step_exception : diffeqs_exception
{
  diffeqs_min_step_exception() :
    diffeqs_exception("Step size fell below min allowed.") {};
};

/******************************************************************************
 *                                                                            *
 * diffeqs                                                                    *
 *                                                                            *
 ******************************************************************************/

template <class Real = double>
class diffeqs
{
  using types = diffeqs_types<Real>;

public:

  using real = typename types::real;
  using cplx = typename types::cplx;
  using real_vec = typename types::real_vec;
  using cplx_vec = typename types::cplx_vec;
  using cplx_vec_range = typename types::cplx_vec_range;
  using matrix = typename types::matrix;
  using connection_t = typename types::connection_t;
  using vector_field_t = typename types::vector_field_t;
  using path_t = typename types::path_t;
  using error_estimator_t = typename types::error_estimator_t;

  // The differential-equation system is defined by the connection matrix of
  // the canonical differential equation and the right-hand side (vector_field)
  // of the differential equations satisfied by the special functions. The
  // boundary condition is read from a stream, which must have the following
  // format.
  //
  // dimx
  // x_1 ... x_dimx
  // order dimbasis nspecial
  // real(f_1_0) imag(f_1_0)
  // real(f_2_0) imag(f_2_0)
  // ...
  // real(f_dimbasis_0) imag(f_dimbasis_0)
  // real(f_1_1) imag(f_1_1)
  // ...
  // real(f_dimbasis_order) imag(f_dimbasis_order)
  // real(special_f_1) imag(special_f_1)
  // ...
  // real(special_f_nspecial) imag(special_f_nspecial)
  //
  // x_i is the i-th component of the boundary point,
  // f_i_j is the coefficient of eps^j in the Laurent expansion of the i-th
  // canonical master integral,
  // special_f_i is the i-th special function.
  diffeqs(connection_t connection, vector_field_t vector_field,
          std::istream& boundary);

  // Evaluate the functions a x. Deformation defines the transformation of a
  // real path point to a complex one during integration. Specifically, x[i]
  // at a given 0 < t < 1 is transformed to
  //
  // x[i] =  x0[i] + (1 + 4*i*(1-t)*deformation[i])*(x[i]-x0[i])
  //
  // where x0 is the boundary point set at construction. For a linear
  // parameterization the deformation is proportional to i*t*(1-t). The result
  // is contained in f, which has a similar ordering to the boundary conditions,
  // but contains only the non-constant functions
  //
  // f_1_1,f_2_1,...,f_dimbasis_1,f_1_2,...,f_dimbasis_2,...,f_dimbasis_order,
  // special_f_1,...,special_f_nspecial
  void evaluate(const real_vec& x, const real_vec& deformation,
                const real& error, cplx_vec& f);

  // The default error estimate is obtained as the maximum of absolute values
  // of function error estimates. This method allows to change it.
  void change_error_estimator(error_estimator_t e) {_error_estimator = e;};

  // The differential equation is integrated along a complex path from the
  // boundary point defined at construction to a requested point during a call
  // to evaluate(). The path is a complex deformation of a real path with
  // parameters specified in the call to evaluate(). The default real path is a
  // line parameterized linearly between the boundary point and endpoint. This
  // method allows to change it.
  void change_path(path_t path) {_path = path;};

  // Integration progress, statistics and timings will be sent to out during
  // a call to evaluate()
  void set_log_stream(std::ostream& out) {_log = &out;};

  // Throw diffeqs_max_steps_exception if set number of steps exceeded
  void set_max_steps(const unsigned n) {_max_steps = n;}

  // Throw diffeqs_max_evals_exception if set number of derivative
  // evaluations exceeded
  void set_max_evals(const unsigned n) {_max_evals = n;}

  // Throw diffeqs_max_time_exception if set time exceeded
  void set_max_time(const unsigned n) {_max_time = n;}

  // Throw diffeqs_min_step_exception if step size below set value
  void set_min_step(const real d) {_min_step = d;}

private:

  static void line(const real_vec& x0,const real_vec& x1, real t,
                   real_vec& x, real_vec& dx);

  static real abs_error(const cplx_vec& f, const cplx_vec& f_err);

private:

  connection_t      _connection;
  vector_field_t    _vector_field;
  path_t            _path = line;
  error_estimator_t _error_estimator = abs_error;

  real_vec _x0;
  unsigned _order;
  cplx_vec _const;
  cplx_vec _f0;

  std::ostream* _log = nullptr;

  unsigned _max_steps = std::numeric_limits<unsigned>::max();
  unsigned _max_evals = std::numeric_limits<unsigned>::max();
  unsigned _max_time =  std::numeric_limits<unsigned>::max();
  real     _min_step =  0;
};

#endif
