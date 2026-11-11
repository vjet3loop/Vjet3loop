// dd.hpp -- double-double arithmetic (about 32 significant digits), after Hida, Li, Bailey (QD).
// Error-free transformations: compile WITHOUT -ffast-math / -fassociative-math, and with -ffp-contract=off.
// Used as the Real type of diffeqs<Real> for quad precision (portable, unlike the x86-only QD library).
#ifndef DD_HPP
#define DD_HPP

#include <cmath>
#include <complex>
#include <istream>
#include <ostream>
#include <string>
#include <type_traits>

struct dd {
    double hi, lo;
    dd() : hi(0), lo(0) {}
    template <class T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type = 0>
    dd(T v) : hi(static_cast<double>(v)), lo(0) {}
    dd(double h, double l) : hi(h), lo(l) {}
};

inline dd quick_two_sum(double a, double b) { double s = a + b; return dd(s, b - (s - a)); }
inline dd two_sum(double a, double b) { double s = a + b, bb = s - a; return dd(s, (a - (s - bb)) + (b - bb)); }
inline dd two_prod(double a, double b) { double p = a * b; return dd(p, std::fma(a, b, -p)); }

inline dd operator+(const dd& a, const dd& b) {
    dd s = two_sum(a.hi, b.hi), t = two_sum(a.lo, b.lo);
    s.lo += t.hi; s = quick_two_sum(s.hi, s.lo); s.lo += t.lo;
    return quick_two_sum(s.hi, s.lo);
}
inline dd operator-(const dd& a) { return dd(-a.hi, -a.lo); }
inline dd operator-(const dd& a, const dd& b) { return a + (-b); }
inline dd operator*(const dd& a, const dd& b) {
    dd p = two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return quick_two_sum(p.hi, p.lo);
}
inline dd operator/(const dd& a, const dd& b) {
    double q1 = a.hi / b.hi; dd r = a - b * dd(q1);
    double q2 = r.hi / b.hi; r = r - b * dd(q2);
    double q3 = r.hi / b.hi;
    return quick_two_sum(q1, q2) + dd(q3);
}
inline dd& operator+=(dd& a, const dd& b) { return a = a + b; }
inline dd& operator-=(dd& a, const dd& b) { return a = a - b; }
inline dd& operator*=(dd& a, const dd& b) { return a = a * b; }
inline dd& operator/=(dd& a, const dd& b) { return a = a / b; }

inline bool operator<(const dd& a, const dd& b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
inline bool operator>(const dd& a, const dd& b) { return b < a; }
inline bool operator<=(const dd& a, const dd& b) { return !(b < a); }
inline bool operator>=(const dd& a, const dd& b) { return !(a < b); }
inline bool operator==(const dd& a, const dd& b) { return a.hi == b.hi && a.lo == b.lo; }
inline bool operator!=(const dd& a, const dd& b) { return !(a == b); }

inline double to_double(const dd& a) { return a.hi + a.lo; }
inline double to_double(double a) { return a; }

// the functions the Bulirsch-Stoer stepper and std::complex need (found by argument-dependent lookup)
inline dd abs(const dd& a) { return a.hi < 0 ? -a : a; }
inline dd sqrt(const dd& a) {
    if (a.hi <= 0) return dd(0);
    double x = 1.0 / std::sqrt(a.hi), ax = a.hi * x;
    return two_sum(ax, (a - two_prod(ax, ax)).hi * (x * 0.5));
}
inline dd pow(const dd& a, const dd& b) { return dd(std::pow(to_double(a), to_double(b))); }   // step control only
inline dd log10(const dd& a) { return dd(std::log10(to_double(a))); }                          // step control only
inline dd hypot(const dd& a, const dd& b) {
    dd x = abs(a), y = abs(b);
    if (x < y) std::swap(x, y);
    if (x.hi == 0) return dd(0);
    dd q = y / x;
    return x * sqrt(dd(1) + q * q);
}
// std::complex<dd>: plain arithmetic (libc++'s generic operators call std::isnan etc. and do not compile for dd).
// Non-template overloads, found by argument-dependent lookup, are preferred to the generic templates.
using cdd = std::complex<dd>;
inline cdd operator+(const cdd& a, const cdd& b) { return cdd(a.real() + b.real(), a.imag() + b.imag()); }
inline cdd operator-(const cdd& a, const cdd& b) { return cdd(a.real() - b.real(), a.imag() - b.imag()); }
inline cdd operator-(const cdd& a) { return cdd(-a.real(), -a.imag()); }
inline cdd operator*(const cdd& a, const cdd& b) {
    return cdd(a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real());
}
inline cdd operator/(const cdd& a, const cdd& b) {
    dd d = b.real() * b.real() + b.imag() * b.imag();
    return cdd((a.real() * b.real() + a.imag() * b.imag()) / d, (a.imag() * b.real() - a.real() * b.imag()) / d);
}
inline cdd operator+(const cdd& a, const dd& b) { return cdd(a.real() + b, a.imag()); }
inline cdd operator+(const dd& a, const cdd& b) { return cdd(a + b.real(), b.imag()); }
inline cdd operator-(const cdd& a, const dd& b) { return cdd(a.real() - b, a.imag()); }
inline cdd operator-(const dd& a, const cdd& b) { return cdd(a - b.real(), -b.imag()); }
inline cdd operator*(const cdd& a, const dd& b) { return cdd(a.real() * b, a.imag() * b); }
inline cdd operator*(const dd& a, const cdd& b) { return cdd(a * b.real(), a * b.imag()); }
inline cdd operator/(const cdd& a, const dd& b) { return cdd(a.real() / b, a.imag() / b); }
inline cdd operator/(const dd& a, const cdd& b) { return cdd(a) / b; }
inline dd abs(const cdd& z) { return hypot(z.real(), z.imag()); }

// decimal input and output (about 32 digits)
inline dd dd_pow10(int e) {
    dd r(1), b(10);
    for (int n = e < 0 ? -e : e; n; n >>= 1, b = b * b) if (n & 1) r = r * b;
    return e < 0 ? dd(1) / r : r;
}
inline dd dd_from_string(const std::string& s) {
    size_t i = 0;
    bool neg = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) neg = s[i++] == '-';
    dd v(0);
    int exp10 = 0;
    bool frac = false;
    for (; i < s.size(); ++i) {
        char c = s[i];
        if (c == '.') { frac = true; continue; }
        if (c < '0' || c > '9') break;
        v = v * dd(10) + dd(c - '0');
        if (frac) --exp10;
    }
    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) exp10 += std::stoi(s.substr(i + 1));
    v = exp10 >= 0 ? v * dd_pow10(exp10) : v / dd_pow10(-exp10);
    return neg ? -v : v;
}
inline std::istream& operator>>(std::istream& in, dd& a) {
    std::string s;
    if (in >> s) a = dd_from_string(s);
    return in;
}
inline std::string dd_to_string(dd a, int digits = 32) {
    if (a.hi == 0) return "0";
    std::string sign = a.hi < 0 ? "-" : "";
    a = abs(a);
    int e = static_cast<int>(std::floor(std::log10(a.hi)));
    dd m = e >= 0 ? a / dd_pow10(e) : a * dd_pow10(-e);
    if (m >= dd(10)) { m = m / dd(10); ++e; }
    if (m < dd(1)) { m = m * dd(10); --e; }
    std::string d;
    for (int k = 0; k < digits; ++k) {
        int q = static_cast<int>(std::floor(m.hi));
        if (q > 9) q = 9;
        if (q < 0) q = 0;
        d += static_cast<char>('0' + q);
        m = (m - dd(q)) * dd(10);
    }
    return sign + d.substr(0, 1) + "." + d.substr(1) + "e" + std::to_string(e);
}
inline std::ostream& operator<<(std::ostream& out, const dd& a) { return out << dd_to_string(a); }

#endif
