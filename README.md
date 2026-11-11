Authors: Thomas Gehrmann, Petr Jakubčík, Andreas von Manteuffel, Cesare Carlo Mella, Lorenzo Tancredi


# zgqq: full-colour finite remainders through 3 loops in QCD for vector boson plus three partons system V q̄ q g

A C++ library that evaluates the **finite remainders, 0 to 3 loops, of the helicity amplitude for
V(p₄) q̄(p₁) q(p₂) g(p₃) in full colour**, as the scalar coefficients ω₁, ω₂, ω₃ of

    M_{L+L} = [13] / (√2 p₄⁴ y z) · ( ω₁ ⟨25⟩²[23][56] + ω₂p₄² ⟨25⟩[36] + ω₃ ⟨12⟩[13]⟨53⟩[36] ),

in the decay region and in the three production regions. Every point gives the ω of M = M_{L+L} and of σ₁₂M (next
section).

The dimensionless invariants x, y, z are defined as

    x = s₁₂/p₄²,   y = s₁₃/p₄²,   z = s₂₃/p₄²,   with s_ij = (p_i + p_j)² and p₄² = (p₁ + p₂ + p₃)²,

so that x + y + z = 1. In the decay all three lie in (0, 1); in a production region one of them is larger than 1 and
the other two are negative (section Choosing points).

The differential equations are solved with the method presented in arXiv:2606.30354.

## Definitions

- **Finite remainders.** Throughout, M stands for the finite remainder of the amplitude: the UV-renormalised amplitude
  with its IR poles removed by the SCET Z factor,

      M = Z⁻¹ M_ren,

  in the MS-bar scheme at the scale μ² = p₄².
- **Expansion in α_s/2π.** With the MS-bar coupling α_s = α_s(μ²) at μ² = p₄²,

      M = Σ_L (α_s/2π)^L M^(L),   L = 0, 1, 2, 3,

  and the program gives the ω of each M^(L): the values of loop L are the coefficient of (α_s/2π)^L.
- **Normalisation.** At tree level ω = (1, 0, 0), for M and for σ₁₂M.
- **Colour.** At each loop, every ω is a sum of colour layers,

      ω^(L) = Σ Nc^a Nf^b NfZ^c ω^(L)_abc,

  with 1, 2, 7 and 14 layers at 0, 1, 2 and 3 loops. Nc is the number of colours and Nf the number of light quark
  flavours; NfZ multiplies the contributions in which V couples to a closed quark loop (see arXiv:2306.10170). The
  program gives each layer separately; with `Nc`, `Nf`, `NfZ` in the input card it also gives their sums.

## Helicities

All helicity amplitudes follow from M = M_{L+L} (quark L, gluon +, e⁺_L e⁻_R):

| (q, g) | e⁺_L e⁻_R | e⁺_R e⁻_L |
|---|---|---|
| L, + | M | σ₅₆M |
| L, − | −P σ₁₂σ₅₆M | −P σ₁₂M |
| R, − | P σ₅₆M | P M |
| R, + | −σ₁₂M | −σ₁₂σ₅₆M |

σ_ij exchanges p_i and p_j, P is parity (⟨ij⟩ ↔ [ij]). P and σ₅₆ leave the ω unchanged (they act on the spinor
structures only). σ₁₂ exchanges y and z: the ω of σ₁₂M are the ω evaluated at (x, z, y). The program gives both for
every point: the lines `omega1 …` are M, the lines `sigma12_omega1 …` are σ₁₂M.

**Production regions.** Each region is the same amplitude crossed into a scattering channel: region `x`
q q̄ → V g, region `y` q g → V q, region `z` q̄ g → V q̄, where the region label refers to the invariant larger
than 1, the partonic centre-of-mass energy ŝ. The jet pseudorapidity η* is that of the outgoing parton in the partonic
centre-of-mass frame, measured along the incoming quark (along the incoming antiquark for q̄ g → V q̄):
η* = ½ ln(z/y) in region `x`, ½ ln(z/x) in region `y`, ½ ln(y/x) in region `z`.

## Build and test

Requirements: a C++17 compiler and the Boost headers (nothing to link). Download the package and build, in its folder,
the library `libzgqq.a` and the program `zgqq`:

```sh
git clone https://github.com/cesarecarlo-mella/Vjet3loop.git
cd Vjet3loop
make
```

Optionally, test against reference values (about a minute):

```sh
make test
```

The default Boost location is `/opt/homebrew/include` (macOS with Homebrew). On Linux the Boost headers are typically
in `/usr/include`, where the compiler finds them by itself: there `make` works as it is. Elsewhere, give the location
(and, if needed, the compiler) to every `make`, for example:

```sh
make BOOST=/path/to/boost/include CXX=g++
```

## Run

Run `zgqq` in the folder of the package, with an input card (it reads the input files in `data/` there; from
another folder, give their path with the keyword `data`). To run the example card `examples/input.txt`:

```sh
./zgqq examples/input.txt
```

The results are written into `zgqq_output.txt` in the current folder. To write them
into another file, give its name:

```sh
./zgqq examples/input.txt output.txt
```

Everything is set in the input card, a text file; `examples/input.txt` is a commented example to copy and edit. Without
a file, the card can come from the shell:

```sh
printf 'points\n-0.5 -0.25\n' | ./zgqq /dev/stdin
```

### The input card

One keyword per line, each at most once; `#` starts a comment. Only `points` is needed; every other keyword has a
default. After `points` come only the points, one per line.

The card is checked strictly. An unknown keyword or value, a missing, extra or non-numeric argument, a keyword given
twice, a point line that is not exactly two numbers, or no points stop the run before anything is computed, with the
line of the card, for example

```
zgqq: error: card.txt, line 3: unknown precision 'quda' (recommended, double or quad)
```

| keyword | values | default | meaning |
|---|---|---|---|
| `data` | a folder | `data` | the input files (relative to where `zgqq` is run) |
| `precision` | `recommended`, `double`, `quad` | `recommended` | see Precision |
| `errors` | `auto`, `all`, `none` | `auto` | `auto`: an estimate where the procedure makes one (two paths); `all`: for every point (a second path where needed); `none`: one path, no estimate |
| `tol` | two numbers | `1e-10 1e-16` | integration tolerance in double and in quad |
| `switch` | a, b | `3.2 43.8657` | switch region of the recommended mode: \|η*\| > a or √(ŝ/p₄²) > b (43.8657 = 4 TeV / M_Z) |
| `criteria` | a, b | `1e-6 1e-2` | in the switch region, quad if error / largest value > a, or error / value > b for one value |
| `Nc`, `Nf`, `NfZ` | a number each | off | if all three are given: also each coefficient at each loop summed over its colour layers |
| `points` | | | from the next line on: one point per line, two numbers |

For example, quad for every point at a tighter tolerance, with error estimates and colour sums:

```sh
printf 'precision quad\nerrors all\ntol 1e-10 1e-20\nNc 3\nNf 5\nNfZ 1\npoints\n-0.6 1.9\n' | ./zgqq /dev/stdin
```

### Choosing points

A point is two numbers, **x y**, with x = s₁₂/p₄², y = s₁₃/p₄² (p₁ = q̄, p₂ = q, p₃ = g) and z = s₂₃/p₄² = 1 − x − y is computed internally.
The region, and with it the process, follows from them:

| (x, y, z) | region | process |
|---|---|---|
| all three in (0, 1) | `decay` | V → q̄ q g |
| x > 1, y < 0, z < 0 | `x` | q q̄ → V g |
| y > 1, x < 0, z < 0 | `y` | q g → V q |
| z > 1, x < 0, y < 0 | `z` | q̄ g → V q̄ |
| two positive, one negative | DIS | not supported |

A point in the DIS region, or on a region boundary (an invariant equal to 0), is rejected with a warning and not
evaluated;


### Output

An example. For the input card

```
Nc 3
Nf 5
NfZ 1
points
-0.5 -0.25
```

the output file contains (some lines left out, marked `...`):

```
# zgqq: 1 points (0 rejected), precision recommended
# per point: a header, then lines: coefficient loop layer Re Im error
# point 1: region z, qbar g -> V qbar, (x, y) = (-0.5, -0.25); double; 0.321 s
omega1 0 Nc^0*Nf^0*NfZ^0 1.0000000000000002 0 nan
omega1 1 Nc^1*Nf^0*NfZ^0 2.2951821682496996 -2.005681419440005 nan
omega1 1 Nc^-1*Nf^0*NfZ^0 -0.83837492562408311 -4.1761595267982248 nan
...
sigma12_omega1 0 Nc^0*Nf^0*NfZ^0 1.0000000000000002 0 nan
sigma12_omega1 1 Nc^1*Nf^0*NfZ^0 0.95438200136637308 -1.4087788152579404 nan
...
sigma12_omega3 3 Nc^-1*Nf^2*NfZ^0 0.85590832084261137 -0.47484365081563229 nan
omega1 0 sum(Nc=3;Nf=5;NfZ=1) 1.0000000000000002 0 nan
...
sigma12_omega3 3 sum(Nc=3;Nf=5;NfZ=1) -3.2627426118424019 199.80667052739417 nan
```

How to read it:

- `# point 1`: the point lies in region `z`, the process q̄ g → V q̄ (the other regions: `decay`, V → q̄ q g; `x`,
  q q̄ → V g; `y`, q g → V q); it was computed in double (`quad` for a point computed in quad), and it took 0.321 s.
- Then one line per value: **coefficient, loop, colour layer, real part, imaginary part, error**. First ω₁ of M at
  loops 0, 1, 2, 3 (1, 2, 7 and 14 colour layers), then ω₂ and ω₃ the same way: 72 lines. Then the same 72 for σ₁₂M,
  as `sigma12_omega1` … `sigma12_omega3`.
- The layer `Nc^a*Nf^b*NfZ^c` is the part multiplying Nc^a Nf^b NfZ^c.
- `nan` in the last column: no error estimate was made for this point (see Precision).
- The last 24 lines appear only when `Nc`, `Nf`, `NfZ` are in the card: each coefficient at each loop summed over its
  colour layers with these values.
- With several points, the header and its lines repeat for every point, numbered as in the card.
- A rejected point has only its header, for example
  `# point 2: rejected (DIS region not supported), (x, y) = (0.5, 0.69999999999999996)`.
- A point computed in quad has 32 digits (its colour sums 17).

To see only the headers and the colour sums: `grep -E '^#|sum' zgqq_output.txt`.

## Precision

The keyword `precision` chooses (all with the default tolerances `tol 1e-10 1e-16`: 1e-10 in double, 1e-16 in quad):

| `precision` | what is computed |
|---|---|
| `recommended` (default) | double (tolerance 1e-10); quad (tolerance 1e-16) only for the points where double is not accurate enough |
| `double` | double (tolerance 1e-10) for every point, with an error estimate |
| `quad` | quad (tolerance 1e-16) for every point, with an error estimate |

M and σ₁₂M of a point are always computed in the same precision.

### The recommended mode

In the production regions double loses digits only at large rapidity or high energy, the **switch region**: a point
with |η*| > 3.2 or √(ŝ/p₄²) > 43.8657 (that is √ŝ > 4 TeV for p₄^2 = M_Z^2). Decay points are never in it.

- **Outside the switch region** the point is computed once in double. No error estimate is made: the `error` column is
  `nan`.
- **Inside the switch region** the point is computed in double with an error estimate. If one of the values has
  - error > 10⁻⁶ relative to the the largest value, or
  - error > 10⁻² relative to the value itself,

  the point is computed again in quad, and the header says `quad`.

The thresholds can be changed with `switch` and `criteria` (keyword table above); 

### The error estimate

The `error` column is an absolute error estimate for each value: 10 × the difference between two evaluations of the
point on two different integration paths. 
`nan` means no estimate was made (the
point was computed on one path only).

| `errors` | `double`, `quad` | `recommended` |
|---|---|---|
| `auto` (default) | estimate for every point | estimate only for the double points in the switch region |
| `all` | estimate for every point | estimate for every point |
| `none` | no estimate (one path, twice as fast) | as `auto` |


### Points near spurious singularities

Some letters vanish on curves inside the regions where the amplitude is finite but its separate pieces are singular:
the iterated integrals have logarithmic branch points there, and on two of the curves,

    x = (1 − y)²   and   x = (1 − z)²,

the rational functions have poles. 

For a point closer to one of these curves than the threshold in the table below, the program evaluates the point's
neighbours P ± k h n on both sides, on the line through the point perpendicular to the curve (n the unit normal,
k = 1, 2, 3 in double, 1 to 5 in quad, h between 10⁻³ and 10⁻² of the point's distance to the edges of the region,
far outside the threshold below), and averages them with symmetric weights,

    f(P) = (15 f̄₁ − 6 f̄₂ + f̄₃) / 10,   f̄ₖ = [f(P + k h n) + f(P − k h n)] / 2   (double),

which reproduces the value at the point up to O(h⁶) (O(h¹⁰) in quad), since the amplitude is smooth across the curve.

| | threshold: distance to the curve, relative to the distance to the edges of the region | width in η* (region `x`) | evaluations per point |
|---|---|---|---|
| double | 10⁻⁵ | about ±5·10⁻⁶ | 6 |
| quad | 10⁻⁴ | about ±5·10⁻⁵ | 10 |

For the curves with logarithmic singularities only, the thresholds are 10⁻⁸ (double) and 10⁻¹² (quad).


### Tolerance

`tol` sets the integration tolerance in double and in quad (default `1e-10 1e-16`).

- **Double**: leave the default. A smaller tolerance does not give more digits, only a slower run.
- **Quad**: a smaller tolerance gives more digits and a slower run.

### Time

Per point, M and σ₁₂M together (they come from one integration), on one core:

| | laptop (Apple M3 Max) | Linux (AMD EPYC 7763 64-Core Processor) |
|---|---|---|
| a point in double | 0.2 s | 0.3 s |
| a point in quad | 5 s | 13 s |

On 10 000 random points in each production region, with 15 GeV < p_T < 2 TeV, |η*| ≤ 4.5, √ŝ ≤ 13.6 TeV and
p₄² = M_Z² (recommended mode, on the Linux machine), the median is 0.3 s per point and the mean 2.8 s: the 20% of the
points that go to quad take most of the time.

Each run of `zgqq` also reads the input files once at the start: put many points in one input card rather
than running `zgqq` once per point.

## Input files (`data/`)

Shipped with the package and read at the start of every run; nothing has to be generated.

| file | content |
|---|---|
| `system.txt` | the word system: 20 letter matrices (decay) |
| `letter_maps.txt` | the letter maps of the regions x, y, z |
| `boundary.txt` | the start values of the integration |
| `readout_<region>.txt` | read-off of the generators from the iterated integrals, per region, in double-double; `readout_decay_sigma12.txt`: the decay read-off used for σ₁₂M |
| `sigma12.txt` | the image of each iterated integral under σ₁₂ (y ↔ z): σ₁₂M comes from the same integration as M |
| `amplitude.txt` | the rational functions and the rows of the 72 values |

## Tests

`make test` (`tests/test.cpp`, reference values in `tests/reference.txt`): the decay Nc³ layer at three loops against
the exact planar values; tree level; the recommended mode and quad against quad at tolerance 1e-20 in the three
regions;

## License

GPL-3.0-or-later (`COPYING`), as the solver it builds on.
