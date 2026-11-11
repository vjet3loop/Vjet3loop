/******************************************************************************
 *                                                                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 *                                                                            *
 ******************************************************************************/

// test.cpp -- make test: zgqq against reference values (tests/reference.txt)
//   1. decay: three loops, layer Nc^3 (double), against the exact planar values: error / largest < 1e-8;
//   2. tree level: omega = (1, 0, 0) for M and sigma12 M at every production point;
//   3. production regions: the recommended mode against quad at tolerance 1e-20 (including points where double
//      fails): error / largest value < 1e-6 and every value to better than 1e-2 of itself;
//   4. production regions, one point per region, quad: against quad at tolerance 1e-20 and against zgqq_alpha's quad
//      alpha and gamma (Schouten relations, all 84 three-loop values): error / largest value < 1e-10.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "dd.hpp"
#include "zgqq.hpp"

using namespace zgqq;
using cd = std::complex<double>;

int main() {
    std::ifstream in("reference.txt");
    if (!in) { std::cerr << "cannot open reference.txt\n"; return 1; }
    Evaluator ev("../data");
    const auto& vals = ev.values();
    const size_t nv = vals.size();
    double worst_planar = 0, worst_tree = 0, worst_a = 0, worst_b = 0, worst_quad = 0, worst_alpha = 0;
    int n_planar = 0, n_tree = 0, n_rec = 0, n_rec_quad = 0, n_quad = 0, n_point = 0;
    Result q;                                                       // the last point evaluated in quad
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line);
        std::string kind;
        double x, y;
        is >> kind >> x >> y;
        if (kind == "planar") {
            Options o;
            o.precision = Precision::dbl;
            o.errors = Errors::none;
            Result r = ev.evaluate(x, y, o);
            std::vector<cd> ref, got;
            for (int i = 0; i < 6; ++i) {
                std::getline(in, line);
                std::istringstream iv(line);
                std::string co;
                double re, im;
                iv >> co >> re >> im;
                for (size_t v = 0; v < nv; ++v)
                    if (vals[v].coefficient == co && vals[v].loop == 3 && vals[v].layer == "Nc^3*Nf^0*NfZ^0") got.push_back(r.value[v]);
                ref.push_back(cd(re, im));
            }
            double big = 0;
            for (const cd& v : ref) big = std::max(big, std::abs(v));
            for (size_t i = 0; i < ref.size(); ++i)
                worst_planar = std::max(worst_planar, i < got.size() ? std::abs(got[i] - ref[i]) / big : INFINITY);
            ++n_planar;
        } else if (kind == "point") {
            std::vector<dd> re(nv), im(nv);
            for (size_t v = 0; v < nv; ++v) {
                std::getline(in, line);
                std::istringstream iv(line);
                std::string a, b;
                iv >> a >> b;
                re[v] = dd_from_string(a); im[v] = dd_from_string(b);
            }
            auto err = [&](const Result& r, size_t v) {                 // |value - reference|
                return std::hypot(to_double(dd(r.value[v].real()) + dd(r.value_lo[v].real()) - re[v]),
                                  to_double(dd(r.value[v].imag()) + dd(r.value_lo[v].imag()) - im[v]));
            };
            double big = 0;
            for (size_t v = 0; v < nv; ++v) big = std::max(big, std::hypot(to_double(re[v]), to_double(im[v])));
            Result r = ev.evaluate(x, y);                               // recommended
            for (size_t v = 0; v < nv; ++v) {
                double e = err(r, v), a = std::hypot(to_double(re[v]), to_double(im[v]));
                worst_a = std::max(worst_a, e / big);
                if (a > 0) worst_b = std::max(worst_b, e / a);
                if (vals[v].loop == 0) {
                    bool one = vals[v].coefficient == "omega1" || vals[v].coefficient == "sigma12_omega1";
                    worst_tree = std::max(worst_tree, std::abs(r.value[v] + r.value_lo[v] - (one ? cd(1) : cd(0))));
                }
            }
            ++n_tree; ++n_rec; n_rec_quad += r.quad;
            if (n_point++ % 4 == 0) {                                   // the first (random) point of each region in quad
                Options o;
                o.precision = Precision::quad;
                o.errors = Errors::none;
                q = ev.evaluate(x, y, o);
                for (size_t v = 0; v < nv; ++v) worst_quad = std::max(worst_quad, err(q, v) / big);
                ++n_quad;
            } else {
                q = Result();
            }
        } else if (kind == "zgqq_alpha") {
            int n = 0;
            is >> n;
            std::map<std::pair<std::string, std::string>, cd> ref;
            for (int i = 0; i < n; ++i) {
                std::getline(in, line);
                std::istringstream iv(line);
                std::string co, layer;
                double a, b;
                iv >> co >> layer >> a >> b;
                ref[{co, layer}] = cd(a, b);
            }
            if (q.value.empty()) continue;                              // compared only at the points done in quad
            double big = 0, e = 0;
            for (const auto& kv : ref) big = std::max(big, std::abs(kv.second));
            int found = 0;
            for (size_t v = 0; v < nv; ++v) {
                auto it = ref.find({vals[v].coefficient, vals[v].layer});
                if (vals[v].loop != 3 || it == ref.end()) continue;
                e = std::max(e, std::abs(q.value[v] + q.value_lo[v] - it->second) / big);
                ++found;
            }
            worst_alpha = std::max(worst_alpha, found == n ? e : INFINITY);
        }
    }
    bool t1 = n_planar > 0 && worst_planar < 1e-8, t2 = n_tree > 0 && worst_tree < 1e-12,
         t3 = n_rec > 0 && worst_a < 1e-6 && worst_b < 1e-2, t4 = n_quad > 0 && worst_quad < 1e-10 && worst_alpha < 1e-10;
    std::printf("1. decay, three loops, layer Nc^3 vs exact planar values, %d points: worst error / largest %.1e  %s\n",
                n_planar, worst_planar, t1 ? "ok" : "FAILED");
    std::printf("2. tree level omega = (1, 0, 0), %d points: worst difference %.1e  %s\n", n_tree, worst_tree, t2 ? "ok" : "FAILED");
    std::printf("3. regions x, y, z, recommended mode vs quad 1e-20, %d points (%d evaluated in quad): worst error / largest "
                "value %.1e, worst value relative to itself %.1e  %s\n", n_rec, n_rec_quad, worst_a, worst_b, t3 ? "ok" : "FAILED");
    std::printf("4. regions x, y, z, quad, %d points: vs quad 1e-20 %.1e, vs zgqq_alpha (Schouten) %.1e  %s\n", n_quad, worst_quad,
                worst_alpha, t4 ? "ok" : "FAILED");
    bool ok = t1 && t2 && t3 && t4;
    std::printf("%s\n", ok ? "PASSED" : "FAILED");
    return ok ? 0 : 1;
}
