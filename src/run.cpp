/******************************************************************************
 *                                                                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 *                                                                            *
 ******************************************************************************/

// run.cpp -- the command-line driver: reads an input card, evaluates its points, writes the 144 values per point
// into the output file (default zgqq_output.txt); messages on the terminal (stderr).
//
//   zgqq input.txt [output.txt]
//
// input card (keywords in any order, '#' comments; see examples/input.txt):
//   data       data                 folder of the input files
//   precision  recommended          recommended, double or quad
//   errors     auto                 auto, all or none
//   tol        1e-10 1e-16          integration tolerances in double and quad
//   switch     3.2 43.8657          switch region: |eta*| > 3.2 or sqrt(s/q^2) > 43.8657 (4 TeV for q = M_Z)
//   criteria   1e-6 1e-2            quad if error / largest value > 1e-6 or error / value > 1e-2
//   Nc         3                    optional, all three together: also each coefficient and loop summed over the
//   Nf         5                      colour layers Nc^a*Nf^b*NfZ^c
//   NfZ        1
//   points                          then one point per line: x y  (x = s12/q^2, y = s13/q^2, z = 1 - x - y)
//
// The card is read strictly: an unknown keyword or value, a missing, extra or non-numeric argument, a keyword given
// twice, a malformed point line or no points stop the run (exit code 1) with the line.  Points in the DIS region (two
// invariants positive, one negative) or on a region boundary are rejected with a warning; the others are evaluated.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "dd.hpp"
#include "zgqq.hpp"

using namespace zgqq;

namespace {

const char* NAME = "zgqq";
const char* DEFAULT_OUTPUT = "zgqq_output.txt";
std::string card;                                   // the input card, for the messages

[[noreturn]] void fail(int line, const std::string& msg) {
    std::cerr << NAME << ": error: " << card << (line > 0 ? ", line " + std::to_string(line) : "") << ": " << msg << "\n";
    std::exit(1);
}

// a whole token as a finite number
bool number(const std::string& s, double& v) {
    char* e;
    v = std::strtod(s.c_str(), &e);
    return e != s.c_str() && *e == '\0' && std::isfinite(v);
}

// why a point cannot be evaluated, or nullptr (x + y + z = 1: all positive is the decay, one positive the production)
const char* rejected(double x, double y) {
    double z = 1 - x - y;
    int pos = (x > 0) + (y > 0) + (z > 0), neg = (x < 0) + (y < 0) + (z < 0);
    if (pos == 2 && neg == 1) return "DIS region not supported";
    if (pos + neg < 3) return "on a region boundary (an invariant is 0)";
    return nullptr;
}

struct Point { double x, y; int line; const char* why; };

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) { std::cerr << "usage: " << NAME << " input.txt [output.txt]\n"; return 1; }
    card = argv[1];
    std::ifstream in(card);
    if (!in) { std::cerr << NAME << ": error: cannot open " << card << "\n"; return 1; }
    std::fprintf(stderr, "%s: reading the input card %s\n", NAME, card.c_str());

    // ---- the input card
    std::string data = "data";
    Options o;
    int ncolour = 0;    // how many of Nc, Nf, NfZ are given
    double Nc = 0, Nf = 0, NfZ = 0;
    std::vector<Point> pts;
    std::set<std::string> seen;
    bool reading_points = false;
    std::string text;
    for (int ln = 1; std::getline(in, text); ++ln) {
        std::istringstream is(text.substr(0, text.find('#')));
        std::vector<std::string> t;
        for (std::string w; is >> w;) t.push_back(w);
        if (t.empty()) continue;
        if (reading_points) {
            double a, b;
            if (t.size() != 2 || !number(t[0], a) || !number(t[1], b))
                fail(ln, "expected a point, two numbers x y (keywords go before 'points')");
            pts.push_back({a, b, ln, rejected(a, b)});
            continue;
        }
        const std::string key = t[0];
        if (!seen.insert(key).second) fail(ln, "'" + key + "' given twice");
        auto args = [&](size_t n, const std::string& what) { if (t.size() != n + 1) fail(ln, "'" + key + "' needs " + what); };
        auto num = [&](size_t i, bool positive) {
            double v;
            if (!number(t[i], v) || (positive && !(v > 0)))
                fail(ln, "'" + key + "': '" + t[i] + "' is not a " + (positive ? "positive " : "") + "number");
            return v;
        };
        if (key == "data") { args(1, "one folder"); data = t[1]; }
        else if (key == "precision") {
            args(1, "one value: recommended, double or quad");
            if (t[1] == "recommended") o.precision = Precision::recommended;
            else if (t[1] == "double") o.precision = Precision::dbl;
            else if (t[1] == "quad") o.precision = Precision::quad;
            else fail(ln, "unknown precision '" + t[1] + "' (recommended, double or quad)");
        } else if (key == "errors") {
            args(1, "one value: auto, all or none");
            if (t[1] == "auto") o.errors = Errors::automatic;
            else if (t[1] == "all") o.errors = Errors::all;
            else if (t[1] == "none") o.errors = Errors::none;
            else fail(ln, "unknown errors '" + t[1] + "' (auto, all or none)");
        } else if (key == "tol") { args(2, "two positive numbers (double, quad)"); o.tol_double = num(1, true); o.tol_quad = num(2, true); }
        else if (key == "switch") { args(2, "two positive numbers"); o.switch_eta = num(1, true); o.switch_sqrts = num(2, true); }
        else if (key == "criteria") { args(2, "two positive numbers"); o.crit_largest = num(1, true); o.crit_self = num(2, true); }
        else if (key == "Nc") { args(1, "one number"); Nc = num(1, false); ++ncolour; }
        else if (key == "Nf") { args(1, "one number"); Nf = num(1, false); ++ncolour; }
        else if (key == "NfZ") { args(1, "one number"); NfZ = num(1, false); ++ncolour; }
        else if (key == "points") { args(0, "nothing on its line (the points follow, one per line)"); reading_points = true; }
        else fail(ln, "unknown keyword '" + key + "'");
    }
    if (ncolour != 0 && ncolour != 3) fail(0, "give all three of Nc, Nf, NfZ (or none)");
    if (!reading_points) fail(0, "no 'points' section");
    if (pts.empty()) fail(0, "no points after 'points'");
    int nrej = 0;
    for (size_t k = 0; k < pts.size(); ++k)
        if (pts[k].why) {
            ++nrej;
            std::fprintf(stderr, "%s: warning: point %zu (line %d), (x, y) = (%.17g, %.17g): %s; not evaluated\n", NAME, k + 1,
                         pts[k].line, pts[k].x, pts[k].y, pts[k].why);
        }
    const int neval = static_cast<int>(pts.size()) - nrej;
    if (neval == 0) fail(0, "no point to evaluate");

    // ---- the evaluation
    const std::string outname = argc > 2 ? argv[2] : DEFAULT_OUTPUT;
    std::ofstream out(outname);
    if (!out) { std::cerr << NAME << ": error: cannot write " << outname << "\n"; return 1; }
    const char* prec = o.precision == Precision::dbl ? "double" : o.precision == Precision::quad ? "quad" : "recommended";
    std::fprintf(stderr, "%s: %s: %zu points (%d rejected), precision %s\n", NAME, card.c_str(), pts.size(), nrej, prec);
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    try {
        std::fprintf(stderr, "%s: reading the input files in %s\n", NAME, data.c_str());
        Evaluator ev(data);
        const auto& vals = ev.values();
        std::fprintf(stderr, "%s: evaluating %d points, results into %s\n", NAME, neval, outname.c_str());
        out << "# zgqq: " << pts.size() << " points (" << nrej << " rejected), precision " << prec << "\n";
        out << "# per point: a header, then lines: coefficient loop layer Re Im error\n";
        int done = 0, nquad = 0;
        const int step = neval >= 10 ? neval / 10 : 1;
        for (size_t k = 0; k < pts.size(); ++k) {
            const Point& p = pts[k];
            char head[512];
            if (p.why) {
                std::snprintf(head, sizeof head, "# point %zu: rejected (%s), (x, y) = (%.17g, %.17g)", k + 1, p.why, p.x, p.y);
                out << head << "\n";
                continue;
            }
            Result r = ev.evaluate(p.x, p.y, o);
            std::snprintf(head, sizeof head, "# point %zu: region %s, %s, (x, y) = (%.17g, %.17g); %s; %.3f s", k + 1,
                          region_name(r.region), process_name(r.region), r.x, r.y, r.quad ? "quad" : "double", r.seconds);
            out << head << "\n";
            for (size_t v = 0; v < vals.size(); ++v) {
                char buf[256];
                if (r.quad) {
                    out << vals[v].coefficient << " " << vals[v].loop << " " << vals[v].layer << " "
                        << dd_to_string(dd(r.value[v].real()) + dd(r.value_lo[v].real())) << " "
                        << dd_to_string(dd(r.value[v].imag()) + dd(r.value_lo[v].imag()));
                    std::snprintf(buf, sizeof buf, " %.3e\n", r.error[v]);
                } else {
                    std::snprintf(buf, sizeof buf, "%s %d %s %.17g %.17g %.3e\n", vals[v].coefficient.c_str(), vals[v].loop,
                                  vals[v].layer.c_str(), r.value[v].real(), r.value[v].imag(), r.error[v]);
                }
                out << buf;
            }
            if (ncolour == 3)
                for (const auto& s : ev.colour_sum(r, Nc, Nf, NfZ)) {
                    char buf[256];
                    std::snprintf(buf, sizeof buf, "%s %d sum(Nc=%g;Nf=%g;NfZ=%g) %.17g %.17g %.3e\n", s.coefficient.c_str(),
                                  s.loop, Nc, Nf, NfZ, s.value.real(), s.value.imag(), s.error);
                    out << buf;
                }
            out.flush();
            ++done; nquad += r.quad;
            if (done % step == 0 || done == neval)
                std::fprintf(stderr, "%s: %d of %d points done (%.0f%%), %.1f s\n", NAME, done, neval, 100.0 * done / neval, elapsed());
        }
        std::fprintf(stderr, "%s: evaluation complete: %d points in %.1f s (%d in quad); results in %s\n", NAME, neval,
                     elapsed(), nquad, outname.c_str());
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s: error: %s (after %.1f s; the points done so far are in %s)\n", NAME, e.what(), elapsed(),
                     outname.c_str());
        return 2;
    }
    return 0;
}
