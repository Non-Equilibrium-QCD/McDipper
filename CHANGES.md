# McDIPPER 1.2.1 (w_weight): changes relative to TheGarzo/McDipper-1.2.1 (c90a3d3)

These are changes I made while using McDIPPER for my own runs. They are offered
for your consideration; adopt whatever you find useful. Each entry says what
the code did before, what it does now, and why.

- **1. Bug fixes** change results and are worth a look in any case.
- **2. Improvements** are optional: new options, output format, robustness and
  speed.
- **3. Minor** are clean-ups with no effect on results.

**Tables built with the original version must be rebuilt** — items 1.1 and 1.2
change the tabulated gluon energy.

---

## 1. Bug fixes

### 1.1 Adjoint momentum-space dipole: spurious cutoff removed

- **Before.** `get_Q2` tabulates an effective saturation scale by searching
  r ∈ [rmin, rmax] for the radius where the dipole drops to e^{-1/4}. At the two
  largest x the dipole never reaches that value inside the window, so the entry
  stayed uninitialised (in practice zero; 161 fundamental and 70 adjoint grid
  points). `AdjointDipole_k` then took `log(Q2_1/Q2_2)` on exactly those points,
  got NaN, and the `tmp!=tmp` guard silently turned it into zero. On top of
  that, the branch x > xscaling and T < Tscaling was hardcoded to zero, which
  looks like a patch applied where the NaN showed up.
- **Now.** Q² is defined by D = exp(−G(x, r̄) T) = e^{-1/4}, Q² = 1/r̄². When
  r̄ lies beyond rmax, the existing large-r extrapolation of `Dipole::G`
  (constant G'', i.e. G = G(x, rmax)·(r/rmax)²) makes the condition solvable
  in closed form: Q² = 4 G(x, rmax) T / rmax² (× CA/CF for the adjoint). No new
  approximation is introduced; it matches the searched value to five digits
  where both exist and is linear in T. The hardcoded zero branch is gone, so the
  adjoint dipole has the same branch structure as the fundamental one, and the
  two scalings (x → xscaling, then T → Tscaling) are applied one after the
  other instead of one blocking the other.
- **Effect.** The jump of the gluon energy across T = Tscaling at forward
  rapidity drops from 9× (η=5), 6.5·10⁷ (η=9) and 1.3·10¹² (η=10) to 1.5–1.75,
  close to the 1.44 at mid-rapidity where the bug never triggered. At forward
  rapidity the table no longer decreases with increasing thickness. Nothing at
  T ≥ Tscaling changes.

### 1.2 Gluon energy integral: full azimuth, converged, no silent zeros

- **Before.** The integral runs over (q, k, φ), φ being the angle between q and
  k. The propagator 1/p² with p = q + k peaks at p → 0, i.e. at φ → π, so φ was
  cut to |φ| < 0.99π. In addition:
  - the φ Jacobian was hardcoded as 2π while the range was 1.98π (1% off);
  - Cuhre ran with maxEval = 50000, which never reached the requested
    EpsRel = 1e-4, and the fail flag was ignored. Compared with a converged
    reference, 214 of 1584 points were off by more than 1%, worst 22.8% at
    η=10, small T — the dilute forward region;
  - p² was computed as q² + k² + 2qk cos φ, which cancels catastrophically near
    φ = π and can go negative. The resulting NaNs were replaced by zero without
    any message.
- **Now.** The integral is still over (q, k, φ), but split into two pieces:
  the bulk φ ∈ [0, 0.99π] and the tail φ ∈ [0.99π, π], each integrated
  separately (the integrand is even in φ, so the result is 2× the half range).
  The peak at φ → π is not a divergence — the measure cancels the 1/p and p_reg
  caps it — but it is only ~p_reg/q wide. Over the full circle Cuhre never
  samples it; in its own interval of width 0.01π the first subdivision finds
  it. On top of that:
  - p² is computed in a numerically stable form;
  - the φ Jacobian equals the integrated range;
  - maxEval = 5·10⁶. Cuhre stops as soon as it reaches the tolerance, so this
    costs time only where the integrand is hard;
  - failed integrations and NaNs are counted and reported for every table.
- **Why not (q, p).** We first changed variables to (q, p), where the 1/p is
  cancelled by the measure and no cut is needed. That broke the symmetry
  e(y, T1, T2) = e(−y, T2, T1): small k becomes a diagonal ridge p ≈ q that
  Cuhre's axis-aligned boxes miss, and 16% of mirror pairs differed by more
  than 1% (worst 98%). In (q, k, φ) both dipoles sit on their own axis and the
  symmetry is exact at the integrand level. The (q, p) result was still useful
  as an independent cross-check: where it is reliable, the final integral
  agrees with it to 0.06–0.4%.
- **Effect.** The gluon energy is larger than with the truncated φ range: by
  ~12% at mid-rapidity and up to ~60% at η = 10.

### 1.3 Isospin and random nuclear orientation

- **Before.** Protons were assigned with `n <= Z`, so every nucleus had one
  proton too many (O: 9p + 7n). The polar Euler angle of the random rotation
  was sampled flat, giving ⟨u_z²⟩ = 0.25 or 0.5 depending on the axis instead
  of 1/3.
- **Now.** `n < Z`. The polar angle is sampled as asin(2u − 1), i.e. uniform in
  cos θ, so orientations are isotropic.
- **Effect.** Event-level output: charge densities, and, through the
  orientation, eccentricities. The bias in the orientation matters for
  deformed and clustered nuclei (O, Ne).

### 1.4 Hotspots: parton cross section computed instead of read off a fit

- **Before.** σ_parton came from `FitInverse_parton_sigma`, a function fitted
  offline to TRENTo, valid only for the Nq, Bq and BG it was made at (Nq = 3,
  BG = 0.16). For any other Nq, σ_NN was left at 0, so no nucleon pair
  interacted.
- **Now.** The program computes σ_parton at start-up for the run's own Nq, Bq,
  BG and √s: it samples 20000 hotspot configuration pairs (fixed seed), builds
  the overlap T̂(b), and bisects until ∫d²b [1 − ⟨exp(−σ_parton T̂)⟩] = σ_inel.
  No fit, about 4 s. σ_inel is reproduced to ~0.2%, better than eight of the
  ten TRENTo fits in the header (0.1–3.5%). Nq = 1 is unchanged.

---

## 2. Improvements

### 2.1 Smooth interpolation of the tables

- **Before.** The tables were interpolated linearly or bilinearly in thickness
  and rapidity: continuous, but with a derivative that steps at every grid
  node, by ~13% across a node for the dipole. Hydro is driven by gradients,
  ∂e/∂r = ∂e/∂T1 · ∂T1/∂r + (1 ↔ 2), so these steps put a staircase into the
  initial pressure gradients.
- **Now.** Catmull-Rom/bicubic interpolation (C¹) in:
  - thickness for the dipole;
  - (T1, T2) and η for the gluon energy table;
  - (η, T) for the quark tables.
  No measurable cost in time or memory. No negative values appear, so cubic
  overshoot is not an issue.

### 2.2 New hotspot mode: `Hotspots_mode`

A nucleon with hotspots is Nq Gaussians of width² Bq. Their centres are drawn
from a Gaussian of variance Br and then recentred on the nucleon position.
Recentring shrinks the realised spread of the centres to Br·(1 − 1/Nq). Br is
not configured directly; it is derived from BG, and the mode decides how:

- `Nucleon_Width` (default, the original behaviour, as in TRENTo). BG is the
  width² of the whole nucleon:
  Bq + Br(1 − 1/Nq) = BG  →  Br = (BG − Bq)/(1 − 1/Nq).
  The nucleon keeps the same size for any Nq and Bq.
- `Hotspot_Spread` (new). BG is the realised spread of the hotspot centres:
  Br(1 − 1/Nq) = BG  →  Br = BG/(1 − 1/Nq).
  The nucleon width² then becomes Bq + BG, so hotspots make the nucleon larger
  instead of redistributing a fixed size.

The two modes describe the same geometry when BG(Spread) = BG(Width) − Bq. For
example, Hotspot_Spread at BG = 0.116 and Nucleon_Width at BG = 0.156 with
Bq = 0.04 both give Br = 0.174. The switch only enters `Config::get_Br`,
σ_parton (1.4) is solved for whichever geometry results, and configs without
the key behave as before.

### 2.3 One HDF5 file per event

- **Before.** Every event wrote a set of ASCII `.dat` files (Global, Nuclei,
  EnergyMoments, charges, …) with six significant digits. A central PbPb event
  was 1.28 GB. A separate job array then converted them to HDF5 with
  `np.loadtxt`, keyed on hardcoded row counts that break silently when the
  grid changes.
- **Now.** `Output: H5: True` writes a single `Event_<id>.h5`:
  - `charges`: shape (NETA·NX·NY, 5) with columns e_g, e_q, n_u, n_d, n_s.
    Rows are ordered η, then x, then y, all ascending, with one chunk per η
    plane. The coordinates are not stored; the grid is rebuilt from the
    config with three linspaces.
  - `Nucleus1` and `Nucleus2`: shape (A, 6) with columns
    i, x, y, z, participant, nucleon type. This replaces `Nuclei`, which
    stacked both nuclei and carried a column for which nucleus a row belongs
    to.
  - `EnergyMoments`, `NuMoments` and `NdMoments` datasets, with the global
    quantities as attributes (`global/...`), plus `weight`.

  The file is compressed with shuffle + deflate, and float32 is more precise
  than the old six-digit text. A central PbPb event is ~210 MB with the
  coordinate columns; dropping them gives a further factor ~2.8 (70 → 25 MB
  on a fine-grid event).
  The `Format` tokens keep their meaning: they now select datasets instead of
  files. With `H5: False`, the old `.dat` output is unchanged.
- **Check.** Against the `.dat` output of the same event, every dataset and
  attribute agrees to 5e-6, the ASCII resolution.

### 2.4 The config a run writes out is valid YAML

- **Before.** `Config::dump` produced a file that neither pyyaml nor McDIPPER
  could read back:
  - `Version:` had no space after the colon, and several keys were wrongly
    indented;
  - booleans were written as 0/1, which the parser rejects;
  - `Weights`, `run_name`, `H5` and the whole `Thickness` section were missing.
- **Now.** The dump is a complete, valid config. A round trip — run, read the
  dumped config, rerun — reproduces the event bit for bit. This is also what
  makes it safe to drop the coordinate columns in 2.3.

### 2.5 Edge cases that used to fail silently

- **Quark moments at small thickness.** The integration was skipped
  for T < `gen_pars::T_tolerance`, which is dT/2 of a grid hardcoded in
  `params_gen.h` (TMax = 30, NT = 301), not of the grid the run configures. On
  the default grid only T = 0 is skipped, as intended. With a finer T grid,
  whole rows with real contributions were set to zero; for example, with
  dT = 0.05 the T = 0.05 row was exactly zero. Now the integration is skipped
  only at T = 0. The default grid is bit-identical.
- **Out-of-range lookups.** A charge lookup with T or η outside the
  table returned zero silently, which blanked out the densest cells of an
  event. It now stops with an error. The η = 0 grid point is located by
  rounding instead of truncation.
- **Impact parameter range.** bMin and bMax were uninitialised, so a config
  that selects sampled b without a `Range` (e.g. `Value` followed by
  `Sampling`) read them from uninitialised memory. They now default to 0.

### 2.6 Tables: reuse, splitting and checks

- **Reuse when the grid covers the run.** Previously the run's T and
  η grids had to equal the table's exactly, so changing TMax, NT, ETA_RANGE or
  NETA forced a rebuild of tens of CPU hours. Now the table only has to cover
  the run. NETA and the η nodes are taken from the table, and the η interval
  is found by searching the tabulated nodes, so the run and table grids may
  differ.
- **Build one rapidity slice per job.** `Output: TableOnly`
  writes the tables and exits, and `merge_table.py` joins the slices. A single
  slice used to give dETA = 1/0 and turn the whole table into NaN.
- **Half the rapidities.** Since e(y, T1, T2) = e(−y, T2, T1) holds
  exactly at the integrand level (1.2), only y ≥ 0 is computed and the rest is
  mirrored: 61 slices instead of 121. The merged table is exactly symmetric.
  Checked against a directly computed y = −5 slice: agreement to 1.3e-4.
- **OpenMP over a rapidity plane.** The
  (T1, T2) points of one rapidity plane are computed in parallel. To make that
  safe, the GSL interpolation accelerators were removed (shared mutable state,
  no measurable speed gain), the Cuhre tallies were made thread-safe, and
  Cuba's own process forking was turned off (it was 30% slower even serially).
  Output is bit-identical to the serial build. One slice takes 2h59 on 32
  threads and 1h32 on 64, about 99% of linear scaling.
- **Resume an interrupted build.** With `Output: Reuse` and
  `ReuseTable: <file>`, a build fills each point it finds in an earlier,
  possibly partial table and integrates only the rest. Points are matched by
  coordinates, so the old table may also be on a different grid. This was
  added after two slices hit the 24 h wall clock at 93% and 95%.
- **Checks on read-back.** A missing,
  truncated or over-long table file now stops with a message naming the file
  and the expected size; previously these cases were, respectively, a
  segfault, a vague error, and silently reading the wrong rows. Grids with
  fewer than 4 nodes (too few for bicubic) are rejected up front. The merge
  stops on any non-finite value and reports the worst violation of the
  mirror symmetry.

---

## 3. Minor

- **Memory and resource handling.** Wrong `delete`
  forms, splines and dipole slices that were never freed, `ChargeMaker` and
  `Config_weights` not released in their destructors, and the integration
  workspace leaking on two exit paths of `BesselIntegrate`. ~3 MB read buffers
  were moved off the stack, and `jT` is clamped so `jT+1` stays inside the
  table. No numerical change.
- **Config parsing.** Comments are stripped. A line without
  `:` no longer resets the section and silently drops every following field.
  A line that fails to parse is now an error.
- **η bin lookup.** Guarded against float rounding at grid nodes.
- **Diagnostic dumps.** The charge test dumps at η = ±3 and ±7 are
  commented out. They made |η| ≥ 7 a hidden requirement of every table.
- **Build.** The Makefile tracks header dependencies, so editing
  `params_*.h` or `routines.h` actually triggers a rebuild.
- **Verbose output.** The per-point table print only appears with
  `Verbose`.
