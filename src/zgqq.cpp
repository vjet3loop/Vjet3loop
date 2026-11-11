/******************************************************************************
 *                                                                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 *                                                                            *
 ******************************************************************************/

// zgqq.cpp -- the evaluator: systems (diffeqs<Real> callbacks), read-off, assembly, paths, spurious curves,
// switch.  Per point one integration: M and sigma12 M (p1 <-> p2) are two read-offs of the same functions.

#include "zgqq.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/numeric/odeint/util/odeint_error.hpp>

#include "dd.hpp"
#include "diffeqs.hpp"

namespace zgqq {

namespace {

using cd = std::complex<double>;
const double P0[2] = {0.3, 0.3};                    // reference point of every region (start of the integration)
const std::vector<std::array<double, 2>> DEFORM = {{0.1, 0.2}, {0.15, 0.25}, {0.25, 0.1}, {-0.15, -0.25}, {0.2, 0.3}, {-0.25, -0.1}};
const std::vector<std::array<double, 2>> RETRY = {{0.15, 0.25}, {0.25, 0.1}, {-0.15, -0.25}};
const double MARGIN = 1e-3;                         // minimal distance (in t) of a letter zero from the open path
const double H_DOUBLE = 1e-3, H_QUAD = 1e-3;        // stencil steps (relative to rho in quad)
const double ERR_FACTOR = 10;                       // error estimate: 10 |difference of two paths|

/******************************************************************************
 * small helpers                                                              *
 ******************************************************************************/

std::string slurp(const std::string& fn) {
    std::ifstream f(fn, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + fn);
    std::ostringstream s;
    s << f.rdbuf();
    return s.str();
}

// a whitespace tokenizer over a file (lines starting with '#' are comments)
struct Tokens {
    std::string buf;
    const char* p;
    explicit Tokens(const std::string& fn) : buf(slurp(fn)), p(buf.c_str()) {}
    void skip() {
        for (;;) {
            while (*p && std::isspace(static_cast<unsigned char>(*p))) ++p;
            if (*p != '#') return;
            while (*p && *p != '\n') ++p;
        }
    }
    std::string word() {
        skip();
        const char* s = p;
        while (*p && !std::isspace(static_cast<unsigned char>(*p))) ++p;
        if (s == p) throw std::runtime_error("unexpected end of data");
        return std::string(s, p);
    }
    long integer() { skip(); char* e; long v = std::strtol(p, &e, 10); p = e; return v; }
    double real() { skip(); char* e; double v = std::strtod(p, &e); p = e; return v; }
    void expect(const char* w) { if (word() != w) throw std::runtime_error(std::string("expected '") + w + "'"); }
};

double absd(const cd& z) { return std::abs(z); }
double absd(const std::complex<dd>& z) { return std::hypot(to_double(z.real()), to_double(z.imag())); }

template <class Real> Real ipow(const Real& b, int e) {
    Real r(1);
    for (int i = 0; i < std::abs(e); ++i) r = r * b;
    return e < 0 ? Real(1) / r : r;
}

/******************************************************************************
 * the alphabet                                                               *
 ******************************************************************************/

// one-forms om_k = dlog W_k along the path (point x, y; tangent xp, yp), and dr for the root r = Sqrt[x y z]
template <class Real>
void one_forms(const std::complex<Real>& x, const std::complex<Real>& y, const std::complex<Real>& xp,
               const std::complex<Real>& yp, const std::complex<Real>& r, const std::array<bool, 21>& active,
               std::array<std::complex<Real>, 21>& om, std::complex<Real>& dr) {
    using C = std::complex<Real>;
    const C c1(Real(1)), c2(Real(2)), I(Real(0), Real(1));
    const C z = c1 - x - y, zp = -(xp + yp);
    const C b = x * y * z, bp = y * z * xp + x * z * yp + x * y * zp;
    auto dlog = [&](const C& W, const C& Wx, const C& Wy) { return (Wx * xp + Wy * yp) / W; };
    for (int k = 1; k <= 20; ++k) {
        if (!active[k]) continue;
        switch (k) {
            case 1: om[k] = xp / x; break;
            case 2: om[k] = yp / y; break;
            case 3: om[k] = zp / z; break;
            case 4: om[k] = -xp / (c1 - x); break;
            case 5: om[k] = -yp / (c1 - y); break;
            case 6: om[k] = (xp + yp) / (x + y); break;
            case 7: case 8: {
                C a = (k == 7) ? x * z : x * y, ap = (k == 7) ? xp * z + x * zp : xp * y + x * yp;
                om[k] = I * (c2 * ap * b - a * bp) / (r * (a * a + b));
                break;
            }
            case 9: om[k] = dlog(c1 - y - x * y, -y, -(c1 + x)); break;
            case 10: om[k] = dlog(c1 - x - x * y, -(c1 + y), -x); break;
            case 11: om[k] = dlog(x - y + y * y, c1, c2 * y - c1); break;
            case 12: om[k] = dlog(-x + x * x + y, c2 * x - c1, c1); break;
            case 13: om[k] = dlog((x - c1) * (x - c1) - y, c2 * (x - c1), -c1); break;
            case 14: om[k] = dlog((x + y) * (x + y) - y, c2 * (x + y), c2 * (x + y) - c1); break;
            case 15: om[k] = dlog((x + y) * (x + y) - x, c2 * (x + y) - c1, c2 * (x + y)); break;
            case 16: om[k] = dlog((c1 - y) * (c1 - y) - x, -c1, -c2 * (c1 - y)); break;
            case 17: om[k] = dlog(x + x * y + y * y, c1 + y, x + c2 * y); break;
            case 18: om[k] = dlog(x * x + y + x * y, c2 * x + y, c1 + x); break;
            case 19: om[k] = dlog((y - c1) * (y - c1) + x * y, y, c2 * (y - c1) + x); break;
            case 20: om[k] = dlog((x - c1) * (x - c1) + x * y, c2 * (x - c1) + y, x); break;
        }
    }
    dr = r * bp / (c2 * b);
}

// the polynomials whose zeros are singular points of the forms (letters 7, 8: the zeros of a^2 + b, a z + y, ...)
cd sing(int k, const cd& x, const cd& y) {
    switch (k) {
        case 1: return x;
        case 2: return y;
        case 3: return 1.0 - x - y;
        case 4: return 1.0 - x;
        case 5: return 1.0 - y;
        case 6: return x + y;
        case 7: return x * (1.0 - x - y) + y;
        case 8: return x * y + 1.0 - x - y;
        case 9: return 1.0 - y - x * y;
        case 10: return 1.0 - x - x * y;
        case 11: return x - y + y * y;
        case 12: return -x + x * x + y;
        case 13: return (x - 1.0) * (x - 1.0) - y;
        case 14: return -y + (x + y) * (x + y);
        case 15: return -x + (x + y) * (x + y);
        case 16: return -x + (1.0 - y) * (1.0 - y);
        case 17: return x + x * y + y * y;
        case 18: return x * x + y + x * y;
        case 19: return (y - 1.0) * (y - 1.0) + x * y;
        case 20: return (x - 1.0) * (x - 1.0) + x * y;
    }
    return 1.0;
}

// letters for the path check, spurious letters (11-16 among them) and their thresholds, per region
struct Profile {
    std::vector<int> sing, spurious;
    std::map<int, double> near;   // closer than near * rho: averaged across the curve (double)
};

/******************************************************************************
 * regions                                                                    *
 ******************************************************************************/

// a region's triangle (X, Y), X, Y > 0, X + Y < 1, Z = 1 - X - Y: the decay (x, y) = (X, Y); a production region A
// (the invariant > 1) has A = 1/Z and the other two invariants -X/Z, -Y/Z.  The integration runs in (X, Y).
std::array<double, 2> physical_point(Region region, double X, double Y) {
    double Z = 1 - X - Y;
    switch (region) {
        case Region::decay: return {X, Y};
        case Region::z: return {-X / Z, -Y / Z};
        case Region::x: return {1 / Z, -X / Z};
        case Region::y: return {-X / Z, 1 / Z};
    }
    return {X, Y};
}

std::pair<Region, std::array<double, 2>> region_point(double x, double y) {
    double z = 1 - x - y;
    if (x > 0 && y > 0 && z > 0) return {Region::decay, {x, y}};
    struct { Region R; double a, b, c; } cand[3] = {{Region::z, z, x, y}, {Region::x, x, y, z}, {Region::y, y, x, z}};
    for (const auto& c : cand)
        if (c.a > 1 && c.b < 0 && c.c < 0) return {c.R, {-c.b / c.a, -c.c / c.a}};
    if ((x > 0) + (y > 0) + (z > 0) == 2 && (x < 0) + (y < 0) + (z < 0) == 1)
        throw std::runtime_error("DIS region not supported (two invariants positive, one negative)");
    throw std::runtime_error("point on a region boundary (an invariant is 0)");
}

// sigma12 M at P: sigma12 ((x, y, z) -> (x, z, y)) maps the iterated integrals based at A at sigma12 P into those
// based at sigma12 A at P, with the letters exchanged (sigma12.txt: the image of each function).  Region x: its
// reference point is its own image; y <-> z: each other's (the same triangle point); the decay: the read-off at the
// image of its reference point (readout_decay_sigma12.txt).  So sigma12 M = this read-off on the images.
const char* sigma12_readout(Region region) {
    switch (region) {
        case Region::decay: return "decay_sigma12";
        case Region::x: return "x";
        case Region::y: return "z";
        case Region::z: return "y";
    }
    return "";
}

// |eta*| = |ln(Y/X)| / 2 and sqrt(s/q^2) = 1/sqrt(Z)
bool in_switch_region(Region region, double X, double Y, const Options& o) {
    if (region == Region::decay) return false;
    return std::fabs(std::log(Y / X)) > 2 * o.switch_eta || 1 - X - Y < 1 / (o.switch_sqrts * o.switch_sqrts);
}

/******************************************************************************
 * paths                                                                      *
 ******************************************************************************/

// roots of c[0] t^d + ... + c[d] (Durand-Kerner)
std::vector<cd> roots(std::vector<cd> c) {
    double mx = 0;
    for (auto& v : c) mx = std::max(mx, std::abs(v));
    size_t s = 0;
    while (s < c.size() && !(std::abs(c[s]) > 1e-14 * mx)) ++s;
    c.erase(c.begin(), c.begin() + s);
    int d = static_cast<int>(c.size()) - 1;
    if (d < 1) return {};
    for (auto& v : c) v /= c[0];
    std::vector<cd> z(d);
    for (int k = 0; k < d; ++k) z[k] = std::pow(cd(0.4, 0.9), k);
    for (int it = 0; it < 500; ++it) {
        double change = 0;
        for (int k = 0; k < d; ++k) {
            cd p = c[0], q = 1;
            for (int i = 1; i <= d; ++i) p = p * z[k] + c[i];
            for (int j = 0; j < d; ++j) if (j != k) q *= z[k] - z[j];
            cd dz = p / q;
            z[k] -= dz;
            change = std::max(change, std::abs(dz));
        }
        if (change < 1e-15) break;
    }
    return z;
}

// smallest distance, in the path parameter, of a zero of a letter from the open segment 0 < t < 1
double path_margin(const Profile& pr, double X, double Y, const std::array<double, 2>& d) {
    // along the path each letter polynomial has degree <= 4 in t: its coefficients (highest first) from its values at
    // t = 0, 1/4, 1/2, 3/4, 1 (inverse Vandermonde matrix, times 3)
    static const double ts[5] = {0, 0.25, 0.5, 0.75, 1};
    static const double Vinv[5][5] = {{32, -128, 192, -128, 32}, {-80, 288, -384, 224, -48}, {70, -208, 228, -112, 22},
                                      {-25, 48, -36, 16, -3}, {3, 0, 0, 0, 0}};
    cd zx[5], zy[5];
    for (int i = 0; i < 5; ++i) {
        double t = ts[i];
        zx[i] = P0[0] + cd(t, 4 * d[0] * t * (1 - t)) * (X - P0[0]);
        zy[i] = P0[1] + cd(t, 4 * d[1] * t * (1 - t)) * (Y - P0[1]);
    }
    double best = std::numeric_limits<double>::infinity();
    for (int k : pr.sing) {
        cd v[5];
        for (int i = 0; i < 5; ++i) v[i] = sing(k, zx[i], zy[i]);
        std::vector<cd> c(5);
        for (int i = 0; i < 5; ++i) { c[i] = 0; for (int j = 0; j < 5; ++j) c[i] += Vinv[i][j] / 3 * v[j]; }
        for (const cd& r : roots(c)) {
            if (std::abs(r - 1.0) < 1e-9) continue;
            double dist = (r.real() >= 0 && r.real() <= 1) ? std::fabs(r.imag()) : std::min(std::abs(r), std::abs(r - 1.0));
            best = std::min(best, dist);
        }
    }
    return best;
}

// the default deformation if the path keeps MARGIN from every zero, else the one with the largest margin;
// second: the best of the others (an independent second evaluation)
std::array<double, 2> choose_deformation(const Profile& pr, double X, double Y, bool second) {
    std::vector<double> m;
    for (const auto& d : DEFORM) m.push_back(path_margin(pr, X, Y, d));
    int first = m[0] >= MARGIN ? 0 : static_cast<int>(std::max_element(m.begin(), m.end()) - m.begin());
    if (!second) return DEFORM[first];
    int best = -1;
    for (int i = 0; i < static_cast<int>(m.size()); ++i)
        if (i != first && (best < 0 || m[i] > m[best])) best = i;
    return DEFORM[best];
}

// closest spurious curve if closer than its threshold: (distance, unit normal), else distance < 0
struct Near { double dist = -1; double n[2] = {0, 0}; };
Near spurious_near(const Profile& pr, double X, double Y, bool quad) {
    double rho = std::min({X, Y, 1 - X - Y}), e = 1e-7 * rho;
    Near best;
    for (int k : pr.spurious) {
        auto f = [&](double a, double b) { return sing(k, a, b).real(); };
        double gx = (f(X + e, Y) - f(X - e, Y)) / (2 * e), gy = (f(X, Y + e) - f(X, Y - e)) / (2 * e);
        double g = std::hypot(gx, gy), d = std::fabs(f(X, Y)) / g;
        double nr = pr.near.at(k);
        double thr = quad ? (nr > 1e-6 ? 1e-4 : 1e-12) : nr;
        if (d < thr * rho && (best.dist < 0 || d < best.dist)) { best.dist = d; best.n[0] = gx / g; best.n[1] = gy / g; }
    }
    return best;
}

/******************************************************************************
 * the systems                                                                *
 ******************************************************************************/

struct Matrices {                                   // decay letter matrices: per letter (row, col, value)
    int n = 0;
    std::array<std::vector<std::array<double, 3>>, 21> M;
};
using LetterMap = std::vector<std::array<int, 3>>;  // (k, j, c): dlog W_k = sum_j c dlog W'_j

Matrices read_matrices(Tokens& t) {
    Matrices m;
    t.expect("functions");
    m.n = static_cast<int>(t.integer());
    for (int l = 1; l <= 20; ++l) {
        t.expect("letter");
        int k = static_cast<int>(t.integer()), cnt = static_cast<int>(t.integer());
        for (int e = 0; e < cnt; ++e) {
            double i = static_cast<double>(t.integer()), j = static_cast<double>(t.integer());
            m.M[k].push_back({i, j, t.real()});
        }
    }
    return m;
}

// one system of one region, integrated by diffeqs<Real>: F' = sum_k om_k M_k F, r' = dr; generators = read-off of F
// for M, then for sigma12 M
template <class Real>
struct Engine {
    using T = diffeqs_types<Real>;
    using C = std::complex<Real>;
    using Readout = std::vector<std::tuple<int, int, C>>;   // generator g += c F_j
    struct Entry { int row, col, letter; Real val; };
    int n = 0, ngen = 0;
    std::array<bool, 21> active{};
    std::vector<Entry> M;                           // all letters in one sparse matrix, sorted by row
    Readout readout, readout_s;                     // M, sigma12 M
    std::unique_ptr<diffeqs<Real>> solver;

    Engine(const Matrices& dec, const LetterMap& cmap, const std::string& boundary, int ng, const Readout& ro,
           const Readout& ro_s) : n(dec.n), ngen(ng), readout(ro), readout_s(ro_s) {
        std::map<std::pair<int, int>, std::map<int, double>> acc;   // (row, col) -> letter -> value
        for (int k = 1; k <= 20; ++k)
            for (const auto& e : dec.M[k]) {
                std::pair<int, int> ij(static_cast<int>(e[0]), static_cast<int>(e[1]));
                for (const auto& kjc : cmap)
                    if (kjc[0] == k) acc[ij][kjc[1]] += kjc[2] * e[2];
            }
        for (const auto& a : acc)
            for (const auto& lv : a.second)
                if (lv.second != 0) {
                    M.push_back({a.first.first, a.first.second, lv.first, Real(lv.second)});
                    active[lv.first] = true;
                }

        auto connection = [](const typename T::cplx_vec&, const typename T::cplx_vec&,
                             const typename T::cplx_vec_range&, typename T::matrix&) {};
        auto vector_field = [this](const typename T::cplx_vec& x, const typename T::cplx_vec& dx,
                                   const typename T::cplx_vec_range& f, typename T::cplx_vec_range& df) {
            const C* fp = &f(0);
            C* dfp = &df(0);
            std::array<C, 21> om;
            C dr;
            one_forms<Real>(x(0), x(1), dx(0), dx(1), fp[n], active, om, dr);
            size_t e = 0;
            for (int i = 0; i < n; ++i) {               // the entries are sorted by row
                C s(Real(0));
                for (; e < M.size() && M[e].row == i; ++e) s += (om[M[e].letter] * M[e].val) * fp[M[e].col];
                dfp[i] = s;
            }
            dfp[n] = dr;
        };
        std::ifstream b(boundary);
        if (!b) throw std::runtime_error("cannot open " + boundary);
        solver.reset(new diffeqs<Real>(connection, vector_field, b));
        // local error relative to the largest function
        solver->change_error_estimator([](const typename T::cplx_vec& f, const typename T::cplx_vec& e) {
            double fm = 0, em = 0;
            for (size_t i = 0; i < f.size(); ++i) { fm = std::max(fm, absd(f(i))); em = std::max(em, absd(e(i))); }
            return Real(em / (fm > 0 ? fm : 1));
        });
        solver->set_min_step(Real(1e-12));
        solver->set_max_steps(100000);
    }

    // the generators at (X, Y), on the path with deformation d: those of M (0 .. ngen-1), then of sigma12 M
    std::vector<C> generators(double X, double Y, const std::array<double, 2>& d, double tol) {
        typename T::cplx_vec f;
        typename T::real_vec x{Real(X), Real(Y)}, def{Real(d[0]), Real(d[1])};
        solver->evaluate(x, def, Real(tol), f);
        std::vector<C> g(2 * ngen, C(Real(0)));
        for (const auto& r : readout) g[std::get<0>(r)] += std::get<2>(r) * f(std::get<1>(r));
        for (const auto& r : readout_s) g[ngen + std::get<0>(r)] += std::get<2>(r) * f(std::get<1>(r));
        return g;
    }
};

/******************************************************************************
 * the amplitude                                                              *
 ******************************************************************************/

struct Row { int r, g, k; cd c; std::complex<dd> cq; };   // c * rational_r * T_g * root^k

// a rational function of (t, u): (sum_j p_j/q_j t^a_j u^b_j) * prod_f factor_f^e_f,
// factors t, u, 1-t-u, t+u, 1-u, 1-t, t^2-t+u, u^2-u+t
struct Rational {
    std::array<int, 8> e;
    std::vector<std::array<long, 4>> num;           // p, q, a, b
};

struct Amplitude {
    std::vector<Rational> rational;
    std::vector<Value> values;
    std::vector<std::vector<Row>> rows;             // per value

    explicit Amplitude(const std::string& fn) {
        Tokens t(fn);
        t.expect("rational");
        rational.resize(t.integer());
        for (auto& r : rational) {
            for (int& e : r.e) e = static_cast<int>(t.integer());
            r.num.resize(t.integer());
            for (auto& m : r.num) for (long& v : m) v = t.integer();
        }
        t.expect("values");
        int nv = static_cast<int>(t.integer());
        for (int v = 0; v < nv; ++v) {
            t.expect("value");
            Value val;
            val.coefficient = t.word();
            val.loop = static_cast<int>(t.integer());
            val.layer = t.word();
            values.push_back(val);
            std::vector<Row> rs(t.integer());
            for (auto& row : rs) {
                row.r = static_cast<int>(t.integer()); row.g = static_cast<int>(t.integer()); row.k = static_cast<int>(t.integer());
                dd re = dd_from_string(t.word()), im = dd_from_string(t.word());
                row.cq = std::complex<dd>(re, im);
                row.c = cd(to_double(re), to_double(im));
            }
            rows.push_back(std::move(rs));
        }
    }

    // the values at the physical point (x, y), root r = rsign Sqrt[|x y z|], from the generators T
    template <class Real>
    std::vector<std::complex<Real>> assemble(const Real& x, const Real& y, int rsign,
                                             const std::vector<std::complex<Real>>& T) const {
        using C = std::complex<Real>;
        using std::sqrt;
        const Real one(1), t = y, u = one - x - y;
        const Real F[8] = {t, u, one - t - u, t + u, one - u, one - t, t * t - t + u, u * u - u + t};
        std::vector<Real> RF(rational.size());
        for (size_t i = 0; i < rational.size(); ++i) {
            Real v(0);
            for (const auto& m : rational[i].num)
                v = v + Real(static_cast<double>(m[0])) / Real(static_cast<double>(m[1])) *
                        ipow(t, static_cast<int>(m[2])) * ipow(u, static_cast<int>(m[3]));
            for (int f = 0; f < 8; ++f) if (rational[i].e[f]) v = v * ipow(F[f], rational[i].e[f]);
            RF[i] = v;
        }
        Real b = x * y * u;
        Real r = sqrt(b < Real(0) ? -b : b);
        if (rsign < 0) r = -r;
        std::vector<C> out(rows.size());
        for (size_t v = 0; v < rows.size(); ++v) {
            C s(Real(0));
            for (const Row& row : rows[v]) {
                C c;
                if constexpr (std::is_same<Real, double>::value) c = row.c; else c = row.cq;
                Real pref = RF[row.r] * (row.k == 0 ? one : ipow(r, row.k));
                s += c * (T[row.g] * pref);
            }
            out[v] = s;
        }
        return out;
    }
};

// stencil weights: sum_k w_k (f(P + k h n) + f(P - k h n)) / 2 = f(P) + O(h^(2L))
template <class Real> Real stencil_weight(int L, int k) {
    double num = 2, den = 1;
    for (int i = 1; i <= L; ++i) num *= i * i;
    for (int i = 1; i <= L - k; ++i) den *= i;
    for (int i = 1; i <= L + k; ++i) den *= i;
    return (k % 2 ? Real(1) : Real(-1)) * Real(num) / Real(den);
}

}  // namespace

const char* region_name(Region r) { return r == Region::decay ? "decay" : r == Region::x ? "x" : r == Region::y ? "y" : "z"; }

const char* process_name(Region r) {
    switch (r) {
        case Region::decay: return "V -> qbar q g";
        case Region::x: return "q qbar -> V g";
        case Region::y: return "q g -> V q";
        case Region::z: return "qbar g -> V qbar";
    }
    return "";
}

/******************************************************************************
 * the evaluator                                                              *
 ******************************************************************************/

struct Evaluator::Impl {
    std::string dir;
    Amplitude amp;
    std::vector<Value> all;                         // the values of M, then those of sigma12 M
    std::map<Region, LetterMap> cmaps;
    std::map<Region, Profile> profiles;
    std::unique_ptr<Matrices> mat;                  // the word system (decay letter matrices), for both precisions
    std::vector<int> sigma;                         // the image of each function under sigma12
    std::map<Region, std::unique_ptr<Engine<double>>> eng_double;
    std::map<Region, std::unique_ptr<Engine<dd>>> eng_quad;

    explicit Impl(const std::string& d) : dir(d), amp(d + "/amplitude.txt") {
        all = amp.values;
        for (Value v : amp.values) { v.coefficient = "sigma12_" + v.coefficient; all.push_back(v); }
        Tokens t(dir + "/letter_maps.txt");
        for (int i = 0; i < 3; ++i) {
            t.expect("region");
            std::string name = t.word();
            Region R = name == "x" ? Region::x : name == "y" ? Region::y : Region::z;
            int cnt = static_cast<int>(t.integer());
            for (int e = 0; e < cnt; ++e) {
                int k = static_cast<int>(t.integer()), j = static_cast<int>(t.integer()), c = static_cast<int>(t.integer());
                cmaps[R].push_back({k, j, c});
            }
        }
        for (int k = 1; k <= 20; ++k) cmaps[Region::decay].push_back({k, k, 1});   // the decay: the identity
        // profiles: the images of the decay letters in use (1-8, 11-16); spurious: 11-16 among them, with the poles of
        // the rational functions on the images of 15, 16
        for (const auto& cm : cmaps) {
            std::array<bool, 21> img{}, pole{};
            for (const auto& kjc : cm.second) {
                int k = kjc[0];
                if ((k >= 1 && k <= 8) || (k >= 11 && k <= 16)) img[kjc[1]] = true;
                if (k == 15 || k == 16) pole[kjc[1]] = true;
            }
            Profile p;
            for (int j = 1; j <= 20; ++j)
                if (img[j]) {
                    p.sing.push_back(j);
                    if (j >= 11 && j <= 16) { p.spurious.push_back(j); p.near[j] = pole[j] ? 1e-5 : 1e-8; }
                }
            profiles[cm.first] = p;
        }
    }

    // a read-off file readout_<name>.txt; with image: each function j replaced by its image under sigma12
    template <class Real>
    typename Engine<Real>::Readout read_readout(const std::string& name, int& ngen, bool image) const {
        Tokens t(dir + "/readout_" + name + ".txt");
        t.expect("readout");
        ngen = static_cast<int>(t.integer());
        typename Engine<Real>::Readout ro(t.integer());
        for (auto& r : ro) {
            int g = static_cast<int>(t.integer()), j = static_cast<int>(t.integer());
            if (image) j = sigma.at(j);
            dd re_hi(t.real()), re_lo(t.real()), im_hi(t.real()), im_lo(t.real());
            dd re = re_hi + re_lo, im = im_hi + im_lo;
            if constexpr (std::is_same<Real, double>::value) r = std::make_tuple(g, j, cd(to_double(re), to_double(im)));
            else r = std::make_tuple(g, j, std::complex<dd>(re, im));
        }
        return ro;
    }

    // the system of region R in precision Real: the same matrices, start values and read-offs for double and quad
    template <class Real> Engine<Real>& engine(Region R) {
        auto& e = std::get<std::map<Region, std::unique_ptr<Engine<Real>>>&>(std::tie(eng_double, eng_quad))[R];
        if (!e) {
            if (!mat) {
                Tokens t(dir + "/system.txt");
                mat.reset(new Matrices(read_matrices(t)));
                Tokens s(dir + "/sigma12.txt");
                s.expect("sigma12");
                sigma.resize(s.integer());
                for (int& j : sigma) j = static_cast<int>(s.integer());
                if (static_cast<int>(sigma.size()) != mat->n) throw std::runtime_error("sigma12.txt does not fit system.txt");
            }
            int ngen = 0, ngen_s = 0;
            auto ro = read_readout<Real>(region_name(R), ngen, false);
            auto ro_s = read_readout<Real>(sigma12_readout(R), ngen_s, true);
            if (ngen != ngen_s) throw std::runtime_error("the read-offs of M and sigma12 M do not fit");
            e.reset(new Engine<Real>(*mat, cmaps.at(R), dir + "/boundary.txt", ngen, ro, ro_s));
        }
        return *e;
    }

    // the values of M and sigma12 M at (X, Y) itself: deformation chosen per point (second path: another one), retries
    // on failure
    template <class Real>
    std::vector<std::complex<Real>> values_at(Region R, double X, double Y, bool second, double tol) {
        const Profile& pr = profiles.at(R);
        std::vector<std::array<double, 2>> tries = {choose_deformation(pr, X, Y, second)};
        tries.insert(tries.end(), RETRY.begin(), RETRY.end());
        for (const auto& d : tries) {
            try {
                Engine<Real>& e = engine<Real>(R);
                std::vector<std::complex<Real>> T = e.generators(X, Y, d, tol);
                std::vector<std::complex<Real>> Ts(T.begin() + e.ngen, T.end());
                T.resize(e.ngen);
                Real x, y;
                int rsign = 1;
                const Real Xr(X), Yr(Y), Zr = Real(1) - Xr - Yr;
                switch (R) {
                    case Region::decay: x = Xr; y = Yr; break;
                    case Region::z: x = -Xr / Zr; y = -Yr / Zr; rsign = -1; break;
                    case Region::x: x = Real(1) / Zr; y = -Xr / Zr; rsign = -1; break;
                    case Region::y: x = -Xr / Zr; y = Real(1) / Zr; rsign = -1; break;
                }
                auto v = amp.assemble<Real>(x, y, rsign, T);
                auto vs = amp.assemble<Real>(x, Real(1) - x - y, rsign, Ts);   // sigma12 M: (x, y) -> (x, z)
                v.insert(v.end(), vs.begin(), vs.end());
                return v;
            } catch (const diffeqs_exception&) {        // the solver's stopping criteria
            } catch (const boost::numeric::odeint::odeint_error&) {   // e.g. no new step size found
            }
        }
        std::ostringstream m;
        auto ph = physical_point(R, X, Y);
        m.precision(17);
        m << "no integration path works at (x, y) = (" << ph[0] << ", " << ph[1] << ") in region " << region_name(R);
        throw std::runtime_error(m.str());
    }

    // the values at (X, Y): near a spurious curve, the symmetric average across it (3 levels in double, 5 in quad)
    template <class Real>
    std::vector<std::complex<Real>> values(Region R, double X, double Y, bool second, double hfac, double tol) {
        const bool quad = !std::is_same<Real, double>::value;
        Near nr = spurious_near(profiles.at(R), X, Y, quad);
        if (nr.dist < 0) return values_at<Real>(R, X, Y, second, tol);
        double rho = std::min({X, Y, 1 - X - Y});
        double h = quad ? hfac * H_QUAD * rho : hfac * std::min(H_DOUBLE, rho / 100);
        int L = quad ? 5 : 3;
        std::vector<std::complex<Real>> sum(all.size(), std::complex<Real>(Real(0)));
        for (int k = 1; k <= L; ++k) {
            auto a = values_at<Real>(R, X + k * h * nr.n[0], Y + k * h * nr.n[1], second, tol);
            auto b = values_at<Real>(R, X - k * h * nr.n[0], Y - k * h * nr.n[1], second, tol);
            Real w = stencil_weight<Real>(L, k) / Real(2);
            for (size_t v = 0; v < sum.size(); ++v) sum[v] += (a[v] + b[v]) * w;
        }
        return sum;
    }

    // the values of M, then those of sigma12 M, at (X, Y) of region R
    template <class Real>
    std::vector<std::complex<Real>> point(Region R, double X, double Y, bool second, double hfac, double tol) {
        return values<Real>(R, X, Y, second, hfac, tol);
    }
};

Evaluator::Evaluator(const std::string& data_dir) : impl(new Impl(data_dir)) {}
Evaluator::~Evaluator() = default;
const std::vector<Value>& Evaluator::values() const { return impl->all; }

Result Evaluator::evaluate(double x, double y, const Options& o) const {
    auto t0 = std::chrono::steady_clock::now();
    auto rp = region_point(x, y);
    const Region R = rp.first;
    const double X = rp.second[0], Y = rp.second[1];
    Impl& I = *impl;
    const size_t nv = I.all.size();
    Result res;
    res.region = R; res.x = x; res.y = y;
    res.value.assign(nv, 0); res.value_lo.assign(nv, 0);
    res.error.assign(nv, std::numeric_limits<double>::quiet_NaN());

    auto run_double = [&](bool with_error) {
        auto v1 = I.point<double>(R, X, Y, false, 1.0, o.tol_double);
        res.value = v1; res.value_lo.assign(nv, 0); res.quad = false;
        if (with_error) {
            auto v2 = I.point<double>(R, X, Y, true, 0.5, o.tol_double);
            for (size_t i = 0; i < nv; ++i) res.error[i] = ERR_FACTOR * std::abs(v1[i] - v2[i]);
        }
    };
    auto run_quad = [&](bool with_error) {
        auto q1 = I.point<dd>(R, X, Y, false, 1.0, o.tol_quad);
        for (size_t i = 0; i < nv; ++i) {
            res.value[i] = cd(q1[i].real().hi, q1[i].imag().hi);
            res.value_lo[i] = cd(q1[i].real().lo, q1[i].imag().lo);
        }
        res.quad = true;
        if (with_error) {
            auto q2 = I.point<dd>(R, X, Y, true, 0.5, o.tol_quad);
            for (size_t i = 0; i < nv; ++i) res.error[i] = ERR_FACTOR * absd(q1[i] - q2[i]);
        }
    };

    if (o.precision == Precision::dbl) run_double(o.errors != Errors::none);
    else if (o.precision == Precision::quad) run_quad(o.errors != Errors::none);
    else if (!in_switch_region(R, X, Y, o)) run_double(o.errors == Errors::all);
    else {
        run_double(true);                           // two paths: the estimate decides
        double vmax = 0;
        bool bad = false;
        for (size_t i = 0; i < nv; ++i) vmax = std::max(vmax, std::abs(res.value[i]));
        for (size_t i = 0; i < nv; ++i)
            if (res.error[i] > o.crit_largest * vmax || res.error[i] > o.crit_self * std::abs(res.value[i])) bad = true;
        if (bad) {
            res.error.assign(nv, std::numeric_limits<double>::quiet_NaN());
            run_quad(o.errors == Errors::all);
        }
    }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

std::vector<ColourSum> Evaluator::colour_sum(const Result& r, double Nc, double Nf, double NfZ) const {
    std::vector<ColourSum> out;
    std::map<std::pair<std::string, int>, size_t> at;
    const auto& vals = values();
    for (size_t i = 0; i < vals.size(); ++i) {
        auto key = std::make_pair(vals[i].coefficient, vals[i].loop);
        if (!at.count(key)) { at[key] = out.size(); out.push_back({key.first, key.second, 0, 0}); }
        int c = 0, f = 0, z = 0;
        std::sscanf(vals[i].layer.c_str(), "Nc^%d*Nf^%d*NfZ^%d", &c, &f, &z);
        double w = std::pow(Nc, c) * std::pow(Nf, f) * std::pow(NfZ, z);
        ColourSum& s = out[at[key]];
        s.value += w * (r.value[i] + r.value_lo[i]);
        s.error += std::fabs(w) * r.error[i];
    }
    return out;
}

}  // namespace zgqq
