# Provenance

This repository was extracted from the locally ignored
`validation/xsec_scan/` directory of a GENIE Generator checkout.

- GENIE repository: `https://github.com/Pabce/GENIE-Generator.git`
- GENIE branch at extraction: `maid-dev`
- GENIE commit at extraction: `c373c89ec1ed44ffa833170807617f776f493bab`
- Extraction date: 2026-08-13

The GENIE checkout contained uncommitted MAID-related changes at extraction
time. Those changes are not part of this repository. Results can depend on the
GENIE build, tune, configuration XML, and data files selected at runtime, so
record the corresponding GENIE commit and working-tree state with published
scanner output.

## Event-fold normalization repair (2026-09-30)

The conditional QE / EM empirical-MEC fold needs corresponding GENIE source
repairs in `src/Framework/Interaction/KPhaseSpace.cxx` and
`src/Physics/Multinucleon/EventGen/MECGenerator.cxx`. The latter replaces the
low-energy EM rejection grid with a physical domain and analytic upper bound;
the weak sampler is unchanged. Keep both scanner and native-generator source
patches and runtime hashes with validation results. Old empirical-MEC event
samples cannot validate the repaired sampler.

In the local development checkout, the `KPhaseSpaceCuts.h/.cxx` source pair
was restored byte-for-byte from the authoritative LXPLUS checkout because
existing shared libraries already referenced it. Rebuilding Framework/Interaction
without those sources left unresolved runtime calls. This is a build consistency
repair, distinct from the two MEC physics fixes. Unrelated MAID working-tree
changes have been preserved. The corresponding MEC/Q2-limit edits were synced
to LXPLUS while retaining its newer phase-space API; its binaries were not rebuilt.

The corresponding local GENIE repairs and restored phase-space helper are
versioned in `Pabce/GENIE-Generator`, commit
`475092fa6` on branch `maid-dev`. The scanner's MEC cache also includes the
513-node refinement used for the 120-MeV production spectra; this raises only
the refinement limit and retains the existing accuracy tolerance and early stop.
