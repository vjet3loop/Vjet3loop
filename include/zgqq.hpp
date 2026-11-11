/******************************************************************************
 *                                                                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 *                                                                            *
 ******************************************************************************/

// zgqq -- full-colour finite remainders, 0 to 3 loops, for V -> qbar(p1) q(p2) g(p3) in the omega basis:
// omega_1..3 of M = M_{L+L} and of sigma12 M (p1 <-> p2), at loops 0-3 in their colour layers Nc^c Nf^f NfZ^z,
// 72 + 72 complex values per point, in the decay region and the three production regions.
//
// The iterated integrals are integrated with diffeqs<Real> (Czakon, arXiv:2606.30354) along a complex-deformed path;
// the generators are read off and assembled linearly with the rational functions.
//
//   zgqq::Evaluator ev("data");
//   zgqq::Result r = ev.evaluate(-0.5, -0.25);      // physical x = s12/q^2, y = s13/q^2
//   r.value[i]                                              // ev.values()[i]: coefficient, loop, colour layer

#ifndef ZGQQ_HPP
#define ZGQQ_HPP

#include <complex>
#include <memory>
#include <string>
#include <vector>

namespace zgqq {

enum class Region { decay, x, y, z };
enum class Precision { recommended, dbl, quad };
enum class Errors { automatic, all, none };   // automatic: an estimate where the procedure makes one

struct Options {
    Precision precision = Precision::recommended;
    Errors errors = Errors::automatic;
    double tol_double = 1e-10, tol_quad = 1e-16;   // integration tolerances (relative to the largest function)
    double switch_eta = 3.2;                       // switch region: |eta*| > switch_eta
    double switch_sqrts = 4000.0 / 91.1876;        //   or sqrt(s/q^2) > switch_sqrts (4 TeV for q = M_Z)
    double crit_largest = 1e-6;                    // quad if error / largest |value| > crit_largest
    double crit_self = 1e-2;                       //   or error / |value| > crit_self for one of the values
};

// what value[i] is: coefficient (omega1-3 of M, sigma12_omega1-3 of sigma12 M), loop (0-3), colour layer
struct Value {
    std::string coefficient;
    int loop;
    std::string layer;
};

struct Result {
    Region region = Region::decay;
    double x = 0, y = 0;                          // the point
    std::vector<std::complex<double>> value;      // value[i]: Evaluator::values()[i]
    std::vector<std::complex<double>> value_lo;   // quad: the rest of the double-double value, else 0
    std::vector<double> error;                    // absolute error estimate, NaN where not estimated
    bool quad = false;                            // evaluated in quad precision
    double seconds = 0;
};

// one coefficient at one loop summed over its colour layers with weights Nc^c Nf^f NfZ^z, and its error
struct ColourSum {
    std::string coefficient;
    int loop;
    std::complex<double> value;
    double error;
};

class Evaluator {
public:
    explicit Evaluator(const std::string& data_dir);
    ~Evaluator();
    Result evaluate(double x, double y, const Options& opt = Options()) const;   // the region follows from (x, y)
    const std::vector<Value>& values() const;
    std::vector<ColourSum> colour_sum(const Result& r, double Nc, double Nf, double NfZ) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

const char* region_name(Region region);
const char* process_name(Region region);   // e.g. "q qbar -> V g"

}  // namespace zgqq

#endif
