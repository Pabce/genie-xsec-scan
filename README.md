# Standalone GENIE Cross-Section Scanner

`xsec_scan` evaluates GENIE cross-section models directly and writes CSV curves.
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
  --qel-fold-nr 80 \
  --qel-fold-np 192 \
  --qel-fold-nphi-p 24 \
  --output out/c12_g18_qel_fold.csv
```

`--qel-fold-density exact-theta` is the default. It returns
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
  --stem c12_ee_repro \
  --outdir out/c12_ee_repro
```

Use `--panel-set fig7` or `--panel-set fig8` for the paper's six remaining C12
energy/angle settings. Give each run a distinct `--outdir` and `--stem`.

It converts the scanner output from `nb/GeV/dcostheta_l` to
`microbarn/sr/GeV` using `1/(2*pi*1000)` and groups component rows into QE,
MEC, RES, DIS, and total curves. The helper defaults to `--fold auto`.
When auto mode selects a Rosenbluth QE fold, it uses
`--qel-fold-density exact-theta` by default: the scanner evaluates the fixed
lab-angle QE energy-conservation condition deterministically instead of
histogramming generated events. Increase `--qel-fold-nr`, `--qel-fold-np`, and
`--qel-fold-nphi-p` for smoother Rosenbluth curves.
