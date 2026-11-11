/******************************************************************************
 *                                                                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 * Copyright (C) 2026 Michal Czakon                                           *
 *                                                                            *
 ******************************************************************************/

#include <chrono>

#include "diffeqs.hpp"

// diffeqs version overrides the original boost code
// check bulirsch_stoer.patch for differences
#include "bulirsch_stoer.hpp"

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace std::chrono;
using namespace boost::numeric::odeint;

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>

using ublas::range;
using ublas::element_prod;
using ublas::axpy_prod;

namespace {

  /****************************************************************************
   *                                                                          *
   * derivative                                                               *
   *                                                                          *
   ****************************************************************************/

  template <class Real>
  class derivative
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

    derivative(connection_t connection, vector_field_t vector_field,
               path_t path, unsigned order, const cplx_vec& constants,
               const real_vec& x0, const real_vec& x1,
               const real_vec& deformation) :
      _connection(connection), _vector_field(vector_field), _path(path),
      _order(order), _const(constants), _x0(x0), _x1(x1),
      _deformation(deformation) {};

    void operator()(cplx_vec& f, cplx_vec& df, real t);

    unsigned n_evaluations() const {return _n_evaluations;}

    real evaluation_time() const {return _evaluation_time.count();}

  private:

    connection_t   _connection;
    vector_field_t _vector_field;
    path_t         _path;

    const unsigned _order;
    const cplx_vec _const;
    const real_vec _x0, _x1;
    const real_vec _deformation;

    unsigned       _n_evaluations = 0;
    duration<double> _evaluation_time = duration<double>::zero();
  };

  template<class Real>
  void derivative<Real>::operator()(cplx_vec& f, cplx_vec& df, real t)
  {
    ++_n_evaluations;

    real_vec r, dr;
    _path(_x0,_x1,t,r,dr);

    // complex deformation

    cplx_vec x(_x0.size()), dx(_x0.size());
    for (unsigned i = 0; i < x.size(); ++i)
      {
        const cplx  z(1, 4*(1-t)*_deformation[i]);
        const cplx dz(0,-4*      _deformation[i]);

         x(i) =  _x0[i]+ z*(r[i]-_x0[i]);
        dx(i) = z*dr[i]+dz*(r[i]-_x0[i]);
      }

    // canonical differential equations

    const unsigned n = _const.size();

    matrix A(n,n);
    const range special_functions_range(_order*n,f.size());
    cplx_vec_range special_functions(f,special_functions_range);
    const auto start{high_resolution_clock::now()};
    _connection(x,dx,special_functions,A);
    _evaluation_time += high_resolution_clock::now()-start;

    cplx_vec_range df1(df,range(0,n));
    axpy_prod(A,_const,df1,true);
    for (unsigned i = 1; i < _order; ++i)
      {
        cplx_vec_range f_range(f,range((i-1)*n,i*n));
        cplx_vec_range df_range(df,range(i*n,(i+1)*n));
        axpy_prod(A,f_range,df_range,true);
      }

    // differential equations for the special functions

    cplx_vec_range d_special_functions(df,special_functions_range);
    _vector_field(x,dx,special_functions,d_special_functions);
  }

  /****************************************************************************
   *                                                                          *
   * observer                                                                 *
   *                                                                          *
   ****************************************************************************/

  template <class Real>
  struct observer
  {
    using types = diffeqs_types<Real>;

    using real = typename types::real;
    using cplx = typename types::cplx;
    using cplx_vec = typename types::cplx_vec;

    void operator()(const cplx_vec& f, const real t);

    ostream*                 _log;
    derivative<Real>&        _der;
    unsigned                 _max_steps;
    unsigned                 _max_evals;
    unsigned                 _max_time;
    real                     _min_step;

    real                     _prev_t = 0;
    unsigned                 _n_steps = 0;
    time_point<steady_clock> _start = steady_clock::now();
  };

  template<class Real>
  void observer<Real>::operator()(const cplx_vec& f , const real t)
  {
    if (++_n_steps > _max_steps) throw diffeqs_max_steps_exception();

    if (_der.n_evaluations() > _max_evals) throw diffeqs_max_evals_exception();

    const duration<double> time{steady_clock::now()-_start};
    if (time.count() > _max_time) throw diffeqs_max_time_exception();

    const real delta_t = t-_prev_t;
    if (_n_steps > 1 && delta_t < _min_step) throw diffeqs_min_step_exception();
    _prev_t = t;

    if (_log != nullptr)
      {
        if (_n_steps == 1)
          *_log << "step      t               delta_t\n\n";

        ostringstream dist;
        dist << t;
        *_log << left << setw(10) << _n_steps-1 << setw(16);
        if (t < 0.5 || dist.str()[0] != '1') *_log << t;
        else *_log << setw(2) << "1-" << setw(14) << 1-t;
        *_log << setw(16) << delta_t << endl;
      }
  }

} // unnamed namespace

/******************************************************************************
 *                                                                            *
 * diffeqs::diffeqs                                                           *
 *                                                                            *
 ******************************************************************************/

template<class Real>
diffeqs<Real>::diffeqs(connection_t connection, vector_field_t vector_field,
                       std::istream& boundary) :
  _connection(connection), _vector_field(vector_field)
{
  if (!boundary) throw runtime_error("Empty boundary stream.");

  unsigned dimx;  // dimension of phase-space point vector
  boundary >> dimx;

  for (unsigned i = 0; i < dimx; ++i)
    {
      real x;
      boundary >> x;
      _x0.push_back(x);
    }

  unsigned dimbasis, nspecial;  // dimension of basis and # of special functions
  boundary >> _order >> dimbasis >> nspecial;
  _f0.resize(_order*dimbasis+nspecial);
  _const.resize(dimbasis);

  for (unsigned i = 0; i < dimbasis; ++i)
    {
      real x, y;
      boundary >> x >> y;
      _const(i) = cplx(x,y);
    }
  for (unsigned o = 1; o <= _order; ++o)
    for (unsigned i = 0; i < dimbasis; ++i)
      {
        real x, y;
        boundary >> x >> y;
        _f0((o-1)*dimbasis+i) = cplx(x,y);
      }
  for (unsigned i = 0; i < nspecial; ++i)
    {
      real x, y;
      boundary >> x >> y;
      _f0(_order*dimbasis+i) = cplx(x,y);
    }
}

/******************************************************************************
 *                                                                            *
 * diffeqs::evaluate                                                          *
 *                                                                            *
 ******************************************************************************/

template<class Real>
void diffeqs<Real>::evaluate(const real_vec& x, const real_vec& d,
                             const real& error, cplx_vec& f)
{
  f = _f0;
  const real zero = 0, one = 1;
  derivative<Real> der(_connection,_vector_field,_path,_order,_const,_x0,x,d);
  bulirsch_stoer<cplx_vec,real> bs(_error_estimator,error);
  observer<Real> obs{_log,der,_max_steps,_max_evals,_max_time,_min_step};

  const auto start{steady_clock::now()};
  const unsigned n_steps = integrate_adaptive(bs,ref(der),f,zero,one,one,obs);
  const auto finish{steady_clock::now()};
  const duration<double> time{finish-start};

  const unsigned n_evals = der.n_evaluations();
  const real eval_time = der.evaluation_time()/n_evals;

  if (_log != nullptr)
    *_log << "\nNumber of steps:                  " << n_steps
          << "\nNumber of derivative evaluations: " << n_evals
          << "\nConnection evaluation time:       " << eval_time
          << "\nElapsed time in seconds:          " << time.count()
          << endl;
}

/******************************************************************************
 *                                                                            *
 * diffeqs::line                                                              *
 *                                                                            *
 ******************************************************************************/

template<class Real>
void diffeqs<Real>::line(const real_vec& x0, const real_vec& x1, real t,
                         real_vec& x, real_vec& dx)
{
  x.clear();
  dx.clear();
  for (unsigned i = 0; i < x0.size(); ++i)
    {
      dx.push_back(x1[i]-x0[i]);
       x.push_back(x0[i]+t*dx[i]);
    }
}

/******************************************************************************
 *                                                                            *
 * diffeqs::abs_error                                                         *
 *                                                                            *
 ******************************************************************************/

template<class Real>
typename diffeqs<Real>::real
diffeqs<Real>::abs_error(const cplx_vec& f, const cplx_vec& f_err)
{
  real f_err_max = 0;
  for (const auto& x : f_err) if (abs(x) > f_err_max) f_err_max = abs(x);
  return f_err_max;
}

/******************************************************************************
 *                                                                            *
 * diffeqs instantiation                                                      *
 *                                                                            *
 ******************************************************************************/

template class diffeqs<double>;

#ifdef QD

#include "dd_real.h"
#include "qd_real.h"

template class diffeqs<dd_real>;
template class diffeqs<qd_real>;

#endif

// zgqq: double-double (include/dd.hpp) as the quad-precision type, portable to any architecture
#include "dd.hpp"

template class diffeqs<dd>;
