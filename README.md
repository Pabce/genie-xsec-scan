# Standalone GENIE Inclusive Cross-Section Scanner

`xsec_scan` evaluates GENIE cross-section models directly (at the inclusive level) and writes CSV curves.
It does not generate events or histogram an event sample.

The scanner is a standalone add-on: it is versioned separately from GENIE and
links against an already-built GENIE installation selected with `$GENIE`.
See [PROVENANCE.md](PROVENANCE.md) for the extraction point and reproducibility
notes.

## Build

From this repository:

```bash
GENIE=/path/to/GENIE ./build.sh
```

`$GENIE` must point to the root of an already-built GENIE installation. The
wrapper requires `root-config` and `$GENIE/bin/genie-config`. If Pythia6 is not
already discoverable, set `PYTHIA6` to its installation prefix or
`PYTHIA6_LIBDIR` to its library directory. For unusual installs, add flags with:

```bash
GENIE=/path/to/GENIE \
  GENIE_EXTRA_CFLAGS="..." GENIE_EXTRA_LDFLAGS="..." \
  ./build.sh
```

The executable is written to:

```bash
./out/xsec_scan
```

### Source layout

The implementation lives under `src/`:

- `main.cxx` is the executable entry point.
- `xsec_scan.hpp` exposes the small command-line runner interface.
- `xsec_scan.cxx` is the unity translation unit.
- `detail/options.inc` owns CLI and batch configuration parsing.
- `detail/runtime_points.inc` owns runtime paths and scan-point construction.
- `detail/kinematics.inc` owns kinematic solving and transformations.
- `detail/qel_folding.inc` and `detail/initial_state_folding.inc` own the two
  deterministic fold implementations.
- `detail/evaluation.inc` owns tune resolution, caching, and grouped component
  evaluation.
- `detail/output_run.inc` owns CSV metadata/output and process workers.

The unity arrangement keeps scanner-only implementation symbols internal while
making the formerly 4,900-line source navigable by responsibility.

## Batch mode

`--batch-config` evaluates several named kinematic scans in one process. This
reuses tune contexts for repeated beam energies and shares cached nuclear-state
quadrature and samples across all energies and angles:

```bash
./out/xsec_scan \
  --batch-config config/c12_fig6_gem21.batch.ini
```

The configuration is dependency-free INI syntax. `[global]` uses CLI option
names without the leading `--`; boolean flags accept `true` or `false`.
Repeated `[scan NAME]` sections accept an optional target PDG plus `fixed` and
`range` entries:

```ini
[global]
mode = tune
tune = GEM21_11a_00_000
event-generator-list = EM
probe = 11
target = 1000060120
observable = d2
diff = Eprime,costheta_l
components = true
fold = auto
jobs = 1
output = out/batch/example.csv

[scan E056_th36]
target = 1000060120
fixed = E=0.56
fixed = costheta_l=0.809016994375
range = Eprime:0.555:0.155:92

[scan E056_th60]
fixed = E=0.56
fixed = costheta_l=0.5
range = Eprime:0.555:0.105:96
```

Multiple `range` entries form a Cartesian product, matching repeated `--scan`
arguments. A scan-local `target` overrides the global target and causes the
scanner to resolve/cache the appropriate GENIE context for that target and
beam energy. Global `fixed` entries apply to every scan; scan-local values can
override them. CLI options after `--batch-config` override scalar global
settings, which is useful for `--jobs`, fold resolution, and output paths.
Full-line comments begin with `#` or `;`.

Batch CSVs add `batch` and `target` columns containing the section name and
resolved target PDG. `row_id` remains globally increasing across sections, and
rows retain configuration-file order.
The bundled three-panel figure-6 example produces 252 points and 14,112 rows:
9.43 s with `--jobs 1` and 2.96 s with `--jobs 10` on the benchmark M1 Pro.
The serial and parallel outputs are identical.

The mixed-target example evaluates identical 0.56 GeV, 60-degree kinematics
for C12 and Fe56:

```bash
./out/xsec_scan \
  --batch-config config/c12_fe56_e056_th60_gem21.batch.ini
```

It produces 10,752 valid rows in 8.26 s serial or 3.86 s with ten workers on
the benchmark machine. Each target slice exactly matches its corresponding
standalone invocation.

At runtime the executable prepends `config` to `GXMLPATH`
when it can infer the path from its own location. This is harmless for stock
GENIE releases where the electromagnetic Q2 floor is compiled into GENIE.

The scanner's `--em-q2-min` option is opt-in. If the linked GENIE exposes a
configurable `EM-Q2-min`, the CSV metadata records
`# em_q2_min_status,applied`. If the linked GENIE is a base/hard-coded build,
the scanner keeps running with GENIE's built-in threshold and records
`# em_q2_min_status,unavailable`; add `--strict` to make that a fatal error.
For student plug-and-play use with base GENIE, leave this option unset. Only
pass it when you know the linked GENIE build supports it:

```bash
--em-q2-min 0.02
```

Do not put `$GENIE/config` ahead of tune directories in `GXMLPATH` for tune
scans. GENIE already adds the default config path after the active tune path.
If `$GENIE/config` is injected too early, tune-specific files such as
`G21_11a/CommonParam.xml` can be shadowed by the base XML.

## Tune Mode

Tune mode uses `GEVGDriver` and sums the cross-section algorithms selected by
the requested GENIE tune and event-generator list.

```bash
./out/xsec_scan \
  --mode tune \
  --tune G18_01a_02_11a \
  --probe 11 \
  --target 1000010010 \
  --event-generator-list EM \
  --observable d2 \
  --diff Eprime,costheta_l \
  --fixed E=2.222 \
  --scan Eprime:0.6:2.0:40 \
  --scan costheta_l:0.75:0.98:30 \
  --xsec-unit nb \
  --output out/em_eprime_costh.csv
```

If no event-generator list is supplied, electron and positron probes default to
`EM`; all other probes default to `Default`.

Some GENIE models use `XSec()` as a differential shape while `Integral()`
carries the tuned absolute normalization. To reproduce event-generator
normalization without random event histogramming, pass `--shape-norm` directly
to the executable:

```bash
--shape-norm auto
```

The scanner deterministically integrates each component's raw
`d2sigma/dEprime/dcostheta_l` shape over the outgoing-lepton plane and compares
it to `Integral()`. `auto` applies the factor only when the mismatch is larger
than `--shape-norm-auto-threshold` (default factor 2). Adjust the midpoint grid
with:

```bash
--shape-norm-ne 120 --shape-norm-ncosth 120
```

### Conditional normalization for G18 event folding

`--fold auto` now defaults to `--fold-normalization event` for the supported
Rosenbluth/uncorrelated LocalFGM QE and EM empirical-MEC paths. At each nuclear
state it divides the differential shape by that state's shape integral, then
uses GENIE's channel `Integral()` for the absolute rate. QE additionally accounts
for the Pauli veto and the generator's retry distribution. This changes both
the shape and rate relative to the historical raw fold; a single constant
rescaling of an old curve cannot reproduce it. Tensor models retain their
existing path; this option does not assert event equivalence for RES or DIS.

For speed, the smooth QE mean density (state integral divided by phase-space
width) and the MEC state integral are interpolated in validated cached tables.
The exact native QE width is applied after its phase-space veto, including
near empty-state thresholds. QE Pauli acceptance is computed once per
nucleus/channel/energy using eight shifted deterministic integration sequences
with analytic azimuth averaging. Tables and acceptance are prepared before
forking scan workers and reused across output points. Progress and convergence
estimates are logged. `--fold-norm-rel-tol 0.0002` controls the normalization
estimate (default 0.02%); it is separate from `--qel-rel-tol`. Neither estimate
is a rigorous error bound. Failure to meet the requested tolerance is an error.
The MEC cache can refine to 513 nodes near low-energy thresholds (needed at
120 MeV); it stops early once the same interpolation tolerance is satisfied.

Empirical MEC uses the moving cluster's conditional W,Q2 shape, the generator's
sampling limits, and its separate QE-derived channel rate. The stationary
`--shape-norm` correction is bypassed for this path to avoid applying the rate
twice. The other models' shape-normalization behavior is unchanged.

**Required GENIE repairs:** the empirical-MEC branch must be reachable in
`KPhaseSpace::Q2Lim`, and `MECGenerator::SelectEmpiricalKinematics` must use a
valid rejection envelope over its physical low-energy domain. The old coarse
envelope grid can miss the domain completely or underestimate the maximum.
For EM `EmpiricalMECPXSec2015`, the repaired generator uses an analytic upper
bound and fails if a sampled value violates it. Scanner-only changes cannot
repair previously generated biased MEC events. Record the exact GENIE source
and library hashes with comparisons; see `PROVENANCE.md`.

`--fold-normalization raw` retains the historical unnormalized folding measure
for diagnostics. Legacy QE backends require this explicit choice; the event
mode supports `auto`, `adaptive-theta`, and `native-q2-reference` only.

### Rosenbluth QE Fold

Rosenbluth QE models expose `dσ/dQ2`, so a requested `d2` curve in
`Eprime,costheta_l` needs an explicit deterministic fold. The scanner keeps
that behavior opt-in:

```bash
./out/xsec_scan \
  --mode tune \
  --tune G18_10a_02_11a \
  --event-generator-list EM \
  --probe 11 \
  --target 1000060120 \
  --observable d2 \
  --diff Eprime,costheta_l \
  --fixed E=0.56 \
  --fixed costheta_l=0.5 \
  --scan Eprime:0.55:0.105:96 \
  --xsec-unit nb \
  --components \
  --qel-bin-fold \
  --qel-fold-density exact-theta \
  --fold-normalization raw \
  --qel-fold-nr 80 \
  --qel-fold-np 192 \
  --qel-fold-nphi-p 24 \
  --output out/c12_g18_qel_fold.csv
```

**Legacy warning:** for G18 Ar40/Ca40/Ca48 production, use the controlled
`adaptive-theta` backend documented below. The legacy fixed-angle Jacobian has
a normalization bias which grid refinement does not fix.

`--qel-fold-density exact-theta --fold-normalization raw` is an explicit legacy option. It returns
`d2σ/dEprime/dcostheta_l` at the requested lab angle using deterministic
quadrature over the configured nuclear state. For comparison with inclusive
electron-scattering figures in `d2σ/dΩdE`, divide by `2*pi`.

### RES, MEC, and DIS Initial-State Folds

Direct `XSec()` evaluation does not run the event-generation modules that add
initial-state motion. For event-generator-equivalent inclusive lepton curves,
use the model-aware mode:

```bash
--fold auto
```

Auto mode follows the standard GENIE generator ordering for the configured
cross-section model:

- Rosenbluth QE is folded because its generator applies nuclear motion before
  sampling lepton kinematics.
- SuSAv2/hadron-tensor QE is evaluated directly because its specialized
  generator samples lepton kinematics from the nuclear tensor.
- RES and DIS are folded because `FermiMover` precedes their kinematics
  generators.
- Empirical MEC is folded because it samples the dinucleon momentum before
  lepton kinematics. SuSAv2 and Nieves/Valencia MEC are evaluated directly
  because they sample lepton kinematics from the nuclear tensor first.

Auto mode also applies GENIE's event-chain threshold and W/Q2 phase-space
vetoes. This removes cross section from folded nuclear states that the event
generator would reject, including unphysical low-energy-transfer RES tails.
It requires tune mode, where the resolved event-generator module chain is
available for inspection, and cannot be combined with manual fold overrides.

For scan-heavy folded curves, use process-isolated parallel workers:

```bash
--jobs auto
```

Workers are forked only after GENIE has built the tune and resolved the
interaction list. This avoids ROOT/GENIE initialization races, shares the
read-only initialized state through copy-on-write, and merges rows back in the
original scan order. The paper-reproduction helper enables `--jobs auto` by
default; use `--jobs 1` for serial diagnostics.

The folded path also caches the nuclear radius/momentum quadrature by target,
hit nucleon, nuclear model, and grid settings; uses the exact Lorentz-invariant
Jacobian instead of six finite-difference solves per state; reuses identical EM
sea-quark/antiquark DIS folds; and skips the exactly zero RS proton amplitudes.
Within each scan point, RES or DIS components sharing a hit nucleon now prepare
the bound nuclear state, boosted lepton kinematics, event phase-space checks,
and Jacobian once per fold sample before evaluating the channel-specific model
terms. The cached nuclear samples do not depend on beam energy or angle, so the
same machinery is reusable by multi-kinematic batch scans.
For the 96-point GEM21 C12 panel at 0.56 GeV and 60 degrees (512 fold states),
the original 338.66 s run dropped to 5.94 s with `--jobs 1` (57.0x) and 3.27 s
with ten workers (103.6x) on an Apple M1 Pro. All 5,376 output rows are
byte-for-byte identical to the preceding optimized implementation and retain
the original statuses and physics values.

For diagnostic/manual studies, select processes explicitly:

```bash
--initial-state-fold RES,MEC,DIS \
--initial-state-fold-samples 512 \
--initial-state-fold-nr 12 \
--initial-state-fold-np 64
```

RES and DIS reproduce the configured nuclear-model sampling plus the
`FermiMover/Default` binding prescription. MEC instead reproduces
`MECGenerator::GenerateFermiMomentum`: it samples two constituent nucleons at
radius zero, sums their momenta, and keeps the di-nucleon cluster on shell.
The fold is deterministic (Halton points), opt-in in `xsec_scan`, and currently
supports `d2` in `Eprime/omega,costheta_l`. Explicit process selection is a
manual override and can be physically inappropriate—for example, SuSAv2 MEC
has no dinucleon code until after lepton kinematics are selected.

The manually selected cross-section fold normally retains model values outside GENIE's
event-generator phase space. To reproduce the event chain's rejection behavior,
including the RES `WLim()` and `Q2Lim_W()` checks applied after Fermi motion,
enable:

```bash
--initial-state-fold-event-phase-space
```

This is deliberately not the default: it changes the deterministic model fold
into an event-sample-equivalent density. The CSV metadata records whether the
veto is on. With `--initial-state-fold-debug`, diagnostics distinguish states
rejected below threshold from states outside the allowed W/Q2 phase space.

## Direct Algorithm Mode

Algorithm mode evaluates one `XSecAlgorithmI` while still using GENIE
interaction-list generators to build the interaction objects.

```bash
./out/xsec_scan \
  --mode alg \
  --tune G18_01a_02_11a \
  --xsec-alg genie::KNOTunedQPMDISPXSec/Default \
  --process DIS-EM \
  --probe 11 \
  --target 1000010010 \
  --observable d2 \
  --diff W,Q2 \
  --fixed E=2.222 \
  --scan W:1.1:2.5:40 \
  --scan Q2:0.05:2.0:40 \
  --xsec-unit nb \
  --output out/dis_w_q2.csv
```

For inclusive algorithms used with resonance-style interaction lists, avoid
overcounting resonance tags with:

```bash
--interaction-sum hit-state
```

The default is `--interaction-sum generated`, which keeps every generated
interaction channel.

## CSV Points

Instead of a Cartesian grid, pass a CSV with a header using the same variable
names accepted by the CLI:

```csv
E,Eprime,costheta_l
2.222,1.500,0.95
2.222,1.250,0.90
```

Then run:

```bash
./out/xsec_scan \
  --mode tune \
  --tune G18_01a_02_11a \
  --probe 11 \
  --target 1000010010 \
  --observable d2 \
  --diff Eprime,costheta_l \
  --points points.csv \
  --xsec-unit nb \
  --output out/points.csv
```

Repeated `--fixed VAR=VALUE` entries fill values that are missing from the CSV.

## Projections

Point densities are the default. Lower-dimensional projections are explicit
midpoint integrations over the intersection of the requested range with
physical lepton phase space. Kinematically forbidden midpoints contribute
zero; malformed or underconstrained inputs remain errors. Integration ranges
must be increasing, and an integrated variable must not also be fixed, scanned,
or supplied in a points CSV row.

For example, project `d2sigma/dW dQ2` onto `dσ/dQ2`:

```bash
./out/xsec_scan \
  --mode alg \
  --tune G18_01a_02_11a \
  --xsec-alg genie::KNOTunedQPMDISPXSec/Default \
  --process DIS-EM \
  --probe 11 \
  --target 1000010010 \
  --observable d1 \
  --diff Q2 \
  --fixed E=2.222 \
  --scan Q2:0.05:2.0:50 \
  --integrate-over W:1.1:2.5:100 \
  --xsec-unit nb \
  --output out/dis_dq2_projected.csv
```

For a fully integrated cross section over a finite rectangle:

```bash
--observable total --integrate-over W:1.1:2.5:100 --integrate-over Q2:0.05:2.0:100
```

Without `--integrate-over`, `--observable total` calls the GENIE model
`Integral()` for each interaction. That is the most literal GENIE total, but it
can be much slower than explicit finite-range midpoint integration for some
models.

## Variables

Accepted names and aliases:

- `E`
- `Eprime`, `Ep`, `E'`
- `Q2`
- `W`
- `omega`, `q0`
- `q3`, `absq`, `|q|`
- `xB`, `x`
- `y`
- `costheta_l`, `costhl`, `ctl`
- `theta_l_deg`
- `Tl`

`E` is the incoming beam energy and is treated as a conditioning variable, not
as a final-state differential variable.

## Output

The output is long-form CSV. Metadata lines start with `#`. Data rows have:

```text
row_id,component,observable,value,unit,status,message,E,Eprime,Q2,W,omega,q3,xB,y,costheta_l,theta_l_deg,Tl
```

By default each point has one `component=total` row. Add `--components` to also
write one row per GENIE interaction/model contribution. A total row is `ok`
only when every required component was evaluated successfully; otherwise its
value is zero and the first component error is reported.

Rows outside the physical phase space are written with `status=invalid` and
`value=0`. A projection is invalid only when its entire integration grid has no
physical midpoint. Because an integrated variable has no unique point value,
it and kinematic quantities derived from it are left blank in projected rows.
Add `--strict` to make the first invalid row abort the job.

## Native Phase Space

By default the scanner chooses a native GENIE phase space per interaction:

- QEL/MEC: `Tl,costheta_l` for two-dimensional densities, `Q2` for one-dimensional densities.
- RES/DIS: `W,Q2` for two-dimensional densities, `W` for one-dimensional densities.

Override this with:

```bash
--native-phase-space auto|WQ2|Tlctl|xy|Q2|W
```

Arbitrary two-variable output pairs are transformed from the native phase space
with a finite-difference Jacobian.

## C12 Inclusive Electron Reproduction

The helper renders the three C12 kinematic settings from paper figure 6 as a
six-panel PNG/PDF:

```bash
./reproduce_c12_ee.py \
  --panel-set fig6 \
  --tunes G21_11a_00_000,G18_10a_02_11a \
  --fold auto \
  --jobs auto \
  --stem c12_ee_repro \
  --outdir out/c12_ee_repro
```

Use `--panel-set fig7` or `--panel-set fig8` for the paper's six remaining C12
energy/angle settings. Give each run a distinct `--outdir` and `--stem`.

It converts the scanner output from `nb/GeV/dcostheta_l` to
`microbarn/sr/GeV` using `1/(2*pi*1000)` and groups component rows into QE,
MEC, RES, DIS, and total curves. The helper defaults to `--fold auto`.
The C12 helper now defaults to the controlled `auto` QE density. Use
`--qel-fold-density exact-theta --fold-normalization raw` only to reproduce legacy calculations; their
Jacobian bias is not fixed by increasing the grid. Historical benchmark
comparisons should also retain the old non-QE folding settings.

### Controlled G18 QE integration

For the uncorrelated `LocalFGM` configuration used by `G18_10a_02_11a`
on C12, Ar40, Ca40 and Ca48, select:

```bash
--fold auto --qel-fold-density adaptive-theta \
--qel-rel-tol 0.001 --qel-abs-tol 1e-12 --qel-max-eval 2000000
```

`adaptive-theta` is a specialized, controlled alternative to the legacy
`exact-theta` midpoint calculation. It preserves GENIE's native piecewise
momentum histogram, caches its radial-shell probabilities, marginalizes the
radius with the local Pauli cut, splits momentum support boundaries and uses
embedded Gauss–Kronrod integration for momentum and azimuth. The density uses
the invariant two-body measure; this also corrects a normalization problem in
the legacy fixed-angle Jacobian. Existing `exact-theta` results should not be
assumed accurate merely because a denser grid looks smoother.

This backend explicitly rejects unsupported nuclear models, densities with multiple interior extrema,
correlated tails, momentum-dependent removal energy and manual
nuclear-state overrides. It does not silently substitute a different physics
model. The default `--qel-fold-density auto` selects this controlled backend and fails
clearly for unsupported configurations. Legacy reproduction requires explicitly
selecting `exact-theta --fold-normalization raw`; there is no silent fallback to its biased Jacobian. `qel-fold-method`, grid sizes
and lattice sample counts do not control this backend.

The absolute tolerance is in GENIE's internal cross-section units, before
conversion to the requested output units. Component CSV messages record the
estimated quadrature error, evaluation count and requested relative tolerance.
These estimates do not include model uncertainty or provide a substitute for
independent convergence checks. A native empty-histogram radial tail is logged;
this specialization rejects tail probability above 1e-5. Budget exhaustion produces `status=nonconverged`; unsupported configurations
produce `status=invalid`. Both cause a process failure with
`--strict`. Progress is logged every 100,000 kernel/phase-space calls. For large
campaigns, use small scan chunks with atomic, hash-verified receipts so completed
chunks can be reused after interruption.

`--qel-fold-density native-q2-reference` is a validation-only finite-bin
estimator. It samples continuous momenta within the native histogram bins and
uses the old generator's independent Q2-to-lab-lepton mapping. Use
`--qel-fold-samples N --qel-reference-seed S` for independently shifted Halton
replicates, and explicitly set both `--qel-bin-width` values. Compare to the
adaptive density *integrated over the same bin*, not to its centre value.
Several seeds are necessary to estimate its numerical uncertainty. Its radial
CDF uses 8192 intervals; it shares the native distribution adapter but not the
adaptive Jacobian or radial Pauli marginal. The older `generator-q2` modes now
include the previously missing struck-proton/neutron count.

Tests:

```bash
c++ -std=c++17 -O2 tests/qe_quadrature_test.cxx -o /tmp/qe_quadrature_test
/tmp/qe_quadrature_test
python3 tests/check_qe_adaptive.py --scanner "$PWD/out/xsec_scan" \
  --output /tmp/qe-adaptive-check-fresh
```

For resumable production, `scripts/run_qe_campaign.py` accepts a JSON list of
`E_MeV`, `theta_deg`, and `omega_max_MeV` settings and scans all three Ar/Ca nuclei
in 25-point chunks. It hashes the executable, settings, GENIE libraries and XML,
locks its output directory, atomically records completed chunks, and reuses
only matching receipts. It preserves unreceipted partial files for inspection
and stops scheduling new work on failure. For example:

```bash
python3 scripts/run_qe_campaign.py --scanner "$PWD/out/xsec_scan" \
  --settings kinematics.json --output out/ar-ca-adaptive-v2 --jobs 4
```
