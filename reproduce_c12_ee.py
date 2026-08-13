#!/usr/bin/env python3
"""Reproduce the C12 inclusive e,e' panel layout with xsec_scan.

The scanner writes d2sigma/dEprime/dcos(theta_l) in nb/GeV.  For an
azimuthally symmetric inclusive electron cross section this script converts to
microbarn/sr/GeV with

    d2sigma/dEprime/dOmega = d2sigma/dEprime/dcos(theta_l) / (2*pi)

and 1000 nb = 1 microbarn.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import pathlib
import re
import subprocess
import sys
from collections import defaultdict


SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
SCANNER = SCRIPT_DIR / "out" / "xsec_scan"
DEFAULT_OUTDIR = SCRIPT_DIR / "out" / "c12_ee_repro"
TARGET_C12 = "1000060120"
PROBE_ELECTRON = "11"
CONVERT_NB_DCOSTH_TO_UB_SR = 1.0 / (2.0 * math.pi * 1000.0)

PANEL_SETS = {
    "fig6": [
        {
            "tag": "E024_th60",
            "E": 0.24,
            "theta": 60.0,
            "omega_min": 0.005,
            "omega_max": 0.115,
            "n": 64,
            "label": "12C, 0.24 GeV, theta = 60 deg",
        },
        {
            "tag": "E056_th36",
            "E": 0.56,
            "theta": 36.0,
            "omega_min": 0.005,
            "omega_max": 0.405,
            "n": 92,
            "label": "12C, 0.56 GeV, theta = 36 deg",
        },
        {
            "tag": "E056_th60",
            "E": 0.56,
            "theta": 60.0,
            "omega_min": 0.005,
            "omega_max": 0.455,
            "n": 96,
            "label": "12C, 0.56 GeV, theta = 60 deg",
        },
    ],
    "fig7": [
        {
            "tag": "E0961_th37p5",
            "E": 0.961,
            "theta": 37.5,
            "omega_min": 0.08,
            "omega_max": 0.70,
            "n": 125,
            "label": "12C, 0.961 GeV, theta = 37.5 deg",
        },
        {
            "tag": "E1299_th37p5",
            "E": 1.299,
            "theta": 37.5,
            "omega_min": 0.12,
            "omega_max": 0.85,
            "n": 147,
            "label": "12C, 1.299 GeV, theta = 37.5 deg",
        },
        {
            "tag": "E2222_th15p54",
            "E": 2.222,
            "theta": 15.54,
            "omega_min": 0.005,
            "omega_max": 0.995,
            "n": 199,
            "label": "12C, 2.222 GeV, theta = 15.54 deg",
        },
    ],
    "fig8": [
        {
            "tag": "E1501_th37p5",
            "E": 1.501,
            "theta": 37.5,
            "omega_min": 0.25,
            "omega_max": 0.95,
            "n": 141,
            "label": "12C, 1.501 GeV, theta = 37.5 deg",
        },
        {
            "tag": "E3595_th16",
            "E": 3.595,
            "theta": 16.0,
            "omega_min": 0.25,
            "omega_max": 1.65,
            "n": 141,
            "label": "12C, 3.595 GeV, theta = 16 deg",
        },
        {
            "tag": "E3595_th20",
            "E": 3.595,
            "theta": 20.0,
            "omega_min": 0.35,
            "omega_max": 0.90,
            "n": 111,
            "label": "12C, 3.595 GeV, theta = 20 deg",
        },
    ],
}

COMPONENTS = ("QE", "MEC", "RES", "DIS")
COLORS = {
    "total": 1,
    "QE": 4,
    "MEC": 2,
    "RES": 3,
    "DIS": 800,
}


def compact_tune(tune: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", tune)


def runtime_env(preserve_gxmlpath: bool) -> dict[str, str]:
    env = os.environ.copy()
    if not env.get("GENIE"):
        raise SystemExit(
            "GENIE is not set; point it at the root of an already-built "
            "GENIE installation"
        )
    env["PATH"] = f"{env['GENIE']}/bin:{env.get('PATH', '')}"
    if not preserve_gxmlpath:
        env.pop("GXMLPATH", None)

    root_lib = ""
    try:
        root_lib = subprocess.check_output(
            ["root-config", "--libdir"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        pass

    lib_parts = [
        f"{env['GENIE']}/lib",
    ]
    pythia6_libdir = env.get("PYTHIA6_LIBDIR", "")
    if not pythia6_libdir and env.get("PYTHIA6"):
        pythia6_root = pathlib.Path(env["PYTHIA6"])
        pythia6_libdir = str(pythia6_root / "lib" if (pythia6_root / "lib").is_dir() else pythia6_root)
    if pythia6_libdir:
        lib_parts.append(pythia6_libdir)
    if root_lib:
        lib_parts.append(root_lib)
    for name in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
        old = env.get(name, "")
        env[name] = ":".join([p for p in lib_parts + [old] if p])
    return env


def run_scan(tune: str, panel: dict, outdir: pathlib.Path, args: argparse.Namespace) -> pathlib.Path:
    out_csv = outdir / f"{compact_tune(tune)}_{panel['tag']}.csv"
    log_file = outdir / f"{compact_tune(tune)}_{panel['tag']}.log"
    if args.reuse and out_csv.exists():
        return out_csv

    e0 = panel["E"]
    costh = math.cos(math.radians(panel["theta"]))
    ep_hi = e0 - panel["omega_min"]
    ep_lo = e0 - panel["omega_max"]
    if ep_lo <= 0.001:
        raise RuntimeError(f"panel {panel['tag']} has nonphysical Eprime minimum")

    cmd = [
        str(SCANNER),
        "--mode",
        "tune",
        "--tune",
        tune,
        "--event-generator-list",
        "EM",
        "--probe",
        PROBE_ELECTRON,
        "--target",
        TARGET_C12,
        "--observable",
        "d2",
        "--diff",
        "Eprime,costheta_l",
        "--fixed",
        f"E={e0:.12g}",
        "--fixed",
        f"costheta_l={costh:.12g}",
        "--scan",
        f"Eprime:{ep_hi:.12g}:{ep_lo:.12g}:{int(panel['n'] * args.nscale)}",
        "--xsec-unit",
        "nb",
        "--components",
        "--jobs",
        str(args.jobs),
        "--output",
        str(out_csv),
    ]
    if args.fold == "auto":
        cmd += ["--fold", "auto"]
    if args.qel_bin_fold or args.fold == "auto":
        cmd += [
            "--qel-bin-width",
            f"Eprime={args.qel_bin_width_energy:.12g}",
            "--qel-bin-width",
            f"costheta_l={args.qel_bin_width_costh:.12g}",
            "--qel-fold-method",
            args.qel_fold_method,
            "--qel-fold-density",
            args.qel_fold_density,
            "--qel-fold-samples",
            str(args.qel_fold_samples),
            "--qel-fold-nr",
            str(args.qel_fold_nr),
            "--qel-fold-np",
            str(args.qel_fold_np),
            "--qel-fold-ncosth-p",
            str(args.qel_fold_ncosth_p),
            "--qel-fold-nphi-p",
            str(args.qel_fold_nphi_p),
            "--qel-fold-ncos0",
            str(args.qel_fold_ncos0),
            "--qel-fold-nphi0",
            str(args.qel_fold_nphi0),
        ]
        if args.qel_bin_fold:
            cmd += ["--qel-bin-fold"]
        if args.no_qel_fold_scan_cache:
            cmd += ["--no-qel-fold-scan-cache"]
        if args.qel_fold_kf is not None:
            cmd += ["--qel-fold-kf", f"{args.qel_fold_kf:.12g}"]
        if args.qel_fold_removal_energy is not None:
            cmd += [
                "--qel-fold-removal-energy",
                f"{args.qel_fold_removal_energy:.12g}",
            ]
    if args.initial_state_fold.lower() != "off":
        cmd += [
            "--initial-state-fold",
            args.initial_state_fold,
            "--initial-state-fold-samples",
            str(args.initial_state_fold_samples),
            "--initial-state-fold-nr",
            str(args.initial_state_fold_nr),
            "--initial-state-fold-np",
            str(args.initial_state_fold_np),
        ]
        if args.initial_state_fold_debug:
            cmd += ["--initial-state-fold-debug"]
        if args.initial_state_fold_event_phase_space:
            cmd += ["--initial-state-fold-event-phase-space"]
    if args.shape_norm != "off":
        cmd += [
            "--shape-norm",
            args.shape_norm,
            "--shape-norm-ne",
            str(args.shape_norm_ne),
            "--shape-norm-ncosth",
            str(args.shape_norm_ncosth),
            "--shape-norm-auto-threshold",
            f"{args.shape_norm_auto_threshold:.12g}",
        ]
    env = runtime_env(args.preserve_gxmlpath)
    with log_file.open("w") as log:
        log.write(" ".join(cmd) + "\n\n")
        subprocess.run(cmd, cwd=SCRIPT_DIR, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    return out_csv


def category(component: str) -> str | None:
    if component == "total":
        return "total"
    if component.startswith("QES"):
        return "QE"
    if component.startswith("MEC"):
        return "MEC"
    if component.startswith("RES"):
        return "RES"
    if component.startswith("DIS"):
        return "DIS"
    return None


def read_series(csv_path: pathlib.Path) -> tuple[dict[str, list[tuple[float, float]]], dict[str, int]]:
    grouped = defaultdict(lambda: defaultdict(float))
    invalid = defaultdict(int)

    with csv_path.open() as f:
        rows = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in rows:
            cat = category(row["component"])
            if cat is None:
                continue
            omega = float(row["omega"])
            status = row["status"]
            value = float(row["value"]) * CONVERT_NB_DCOSTH_TO_UB_SR
            if status != "ok":
                invalid[cat] += 1
            grouped[cat][omega] += value

    series = {}
    for cat, points in grouped.items():
        series[cat] = sorted(points.items())
    return series, dict(invalid)


def write_summary(
    outdir: pathlib.Path,
    scan_table: dict,
    invalid_table: dict,
    tunes: list[str],
    panels: list[dict],
    args: argparse.Namespace,
) -> None:
    summary = outdir / "summary.txt"
    with summary.open("w") as f:
        f.write(f"C12 inclusive e,e' analytic scan ({args.panel_set})\n")
        f.write(f"scanner: {SCANNER}\n")
        f.write(f"target: {TARGET_C12}\n")
        f.write(f"probe: {PROBE_ELECTRON}\n")
        f.write(f"scanner jobs: {args.jobs}\n")
        f.write("unit conversion: ub/sr/GeV = nb/GeV/dcostheta / (2*pi*1000)\n\n")
        f.write("panel omega ranges:\n")
        for panel in panels:
            f.write(
                f"  {panel['tag']}: {panel['omega_min']} to "
                f"{panel['omega_max']} GeV, n={int(panel['n'] * args.nscale)}\n"
            )
        f.write("\n")
        if args.qel_bin_fold or args.fold == "auto":
            f.write(
                f"Rosenbluth QE fold: mode={args.fold}, "
                f"dE={args.qel_bin_width_energy}, "
                f"dcostheta={args.qel_bin_width_costh}, "
                f"method={args.qel_fold_method}, "
                f"density={args.qel_fold_density}, "
                f"samples={args.qel_fold_samples}, "
                f"nr={args.qel_fold_nr}, "
                f"np={args.qel_fold_np}, "
                f"ncosth_p={args.qel_fold_ncosth_p}, "
                f"nphi_p={args.qel_fold_nphi_p}, "
                f"ncos0={args.qel_fold_ncos0}, "
                f"nphi0={args.qel_fold_nphi0}, "
                f"scan_cache={'off' if args.no_qel_fold_scan_cache else 'on'}, "
                f"kF={args.qel_fold_kf if args.qel_fold_kf is not None else 'auto'}, "
                "removalE="
                f"{args.qel_fold_removal_energy if args.qel_fold_removal_energy is not None else 'auto'}\n\n"
            )
        if args.shape_norm != "off":
            f.write(
                f"Shape normalization: {args.shape_norm}, "
                f"grid={args.shape_norm_ne}x{args.shape_norm_ncosth}, "
                f"auto-threshold={args.shape_norm_auto_threshold}\n\n"
            )
        if args.fold == "auto" or args.initial_state_fold.lower() != "off":
            f.write(
                "Initial-state fold: "
                f"{args.fold if args.fold == 'auto' else args.initial_state_fold}, "
                f"samples={args.initial_state_fold_samples}, "
                f"nr={args.initial_state_fold_nr}, "
                f"np={args.initial_state_fold_np}, "
                "event-phase-space="
                f"{'on' if args.fold == 'auto' or args.initial_state_fold_event_phase_space else 'off'}\n\n"
            )
        for tune in tunes:
            f.write(f"{tune}\n")
            for panel in panels:
                csv_path = scan_table[(tune, panel["tag"])]
                f.write(f"  {panel['tag']}: {csv_path.name}\n")
                invalid = invalid_table.get((tune, panel["tag"]), {})
                if invalid:
                    f.write(f"    invalid component rows: {invalid}\n")


def plot(
    scan_table: dict,
    tunes: list[str],
    panels: list[dict],
    outdir: pathlib.Path,
    stem: str,
    args: argparse.Namespace,
) -> tuple[pathlib.Path, pathlib.Path]:
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFont(132, "XYZ")
    ROOT.gStyle.SetLabelFont(132, "XYZ")

    canvas = ROOT.TCanvas("c12_ee_repro", "c12_ee_repro", 1120, 930)
    canvas.Divide(2, 3, 0.003, 0.003)
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(132)

    all_series = {}
    invalid_table = {}
    ymax_by_row = []
    for row, panel in enumerate(panels):
        ymax = 0.0
        for tune in tunes:
            series, invalid = read_series(scan_table[(tune, panel["tag"])])
            all_series[(tune, panel["tag"])] = series
            invalid_table[(tune, panel["tag"])] = invalid
            for cat in ("total",) + COMPONENTS:
                for _, y in series.get(cat, []):
                    ymax = max(ymax, y)
        ymax_by_row.append(max(1.0, ymax * 1.18))

    graph_refs = []
    for row, panel in enumerate(panels):
        for col, tune in enumerate(tunes):
            pad = canvas.cd(row * 2 + col + 1)
            pad.SetTicks(1, 1)
            pad.SetLeftMargin(0.12 if col == 0 else 0.03)
            pad.SetRightMargin(0.03)
            pad.SetTopMargin(0.12 if row == 0 else 0.04)
            pad.SetBottomMargin(0.16)

            xmin = panel["omega_min"]
            xmax = panel["omega_max"]
            frame = pad.DrawFrame(xmin, 0.0, xmax, ymax_by_row[row])
            frame.GetXaxis().SetTitle("Energy Transfer [GeV]")
            frame.GetYaxis().SetTitle("d^{2}#sigma/d#Omega dE [#mub/sr/GeV]" if col == 0 else "")
            frame.GetXaxis().CenterTitle(True)
            frame.GetXaxis().SetTitleSize(0.052)
            frame.GetYaxis().SetTitleSize(0.075)
            frame.GetXaxis().SetLabelSize(0.052)
            frame.GetYaxis().SetLabelSize(0.065 if col == 0 else 0.0)
            frame.GetYaxis().SetTitleOffset(0.78)
            frame.GetXaxis().SetTitleOffset(1.05)
            frame.GetXaxis().SetNdivisions(505)
            frame.GetYaxis().SetNdivisions(505)

            if row == 0:
                latex.SetTextAlign(22)
                latex.SetTextSize(0.095)
                latex.DrawLatex(0.5, 0.965, tune.replace("_00_000", ""))

            latex.SetTextAlign(12)
            latex.SetTextSize(0.06)
            row_label_y = 0.84 if row == 0 else 0.92
            latex.DrawLatex(0.16 if col == 0 else 0.06, row_label_y, panel["label"])

            series = all_series[(tune, panel["tag"])]
            for cat in ("QE", "MEC", "RES", "DIS", "total"):
                pts = series.get(cat, [])
                if not pts:
                    continue
                graph = ROOT.TGraph(len(pts))
                for i, (x, y) in enumerate(pts):
                    graph.SetPoint(i, x, y)
                graph.SetLineColor(COLORS[cat])
                graph.SetLineWidth(3 if cat == "total" else 2)
                graph.Draw("L SAME")
                graph_refs.append(graph)

            legend_row = 2
            if row == legend_row and col == 1:
                legend_box = (
                    (0.06, 0.52, 0.32, 0.90)
                    if args.panel_set == "fig8"
                    else (0.58, 0.56, 0.96, 0.92)
                )
                legend = ROOT.TLegend(*legend_box)
                legend.SetBorderSize(0)
                legend.SetFillColor(0)
                legend.SetFillStyle(1001)
                legend.SetTextFont(132)
                legend.SetTextSize(0.055)
                for cat in ("total", "QE", "MEC", "RES", "DIS"):
                    dummy = ROOT.TGraph()
                    dummy.SetLineColor(COLORS[cat])
                    dummy.SetLineWidth(3 if cat == "total" else 2)
                    legend.AddEntry(dummy, cat, "l")
                    graph_refs.append(dummy)
                legend.Draw()
                graph_refs.append(legend)

    canvas.Update()
    png = outdir / f"{stem}.png"
    pdf = outdir / f"{stem}.pdf"
    canvas.SaveAs(str(png))
    canvas.SaveAs(str(pdf))
    write_summary(outdir, scan_table, invalid_table, tunes, panels, args)
    return png, pdf


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panel-set",
        choices=tuple(PANEL_SETS),
        default="fig6",
        help="Paper figure kinematics to scan (default: fig6).",
    )
    parser.add_argument(
        "--tunes",
        default="G21_11a_00_000,G18_10a_02_11a",
        help="Comma-separated pair of GENIE tune IDs for the two plot columns.",
    )
    parser.add_argument("--outdir", type=pathlib.Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--stem", default="c12_ee_repro")
    parser.add_argument(
        "--jobs",
        default="auto",
        help="Forked scanner workers per panel (default: auto).",
    )
    parser.add_argument("--reuse", action="store_true", help="Reuse existing panel CSVs.")
    parser.add_argument("--plot-only", action="store_true", help="Do not run xsec_scan.")
    parser.add_argument(
        "--nscale",
        type=float,
        default=1.0,
        help="Scale the hard-coded number of points per panel.",
    )
    parser.add_argument(
        "--preserve-gxmlpath",
        action="store_true",
        help="Keep the caller's GXMLPATH instead of clearing it for tune-safe default lookup.",
    )
    parser.add_argument(
        "--fold",
        choices=("auto", "manual", "off"),
        default="auto",
        help=(
            "Select generator-path-aware folding (default: auto). Auto folds "
            "Rosenbluth QE, Empirical MEC, RES, and DIS only when their "
            "standard event-generator path samples after initial-state motion."
        ),
    )
    parser.add_argument(
        "--qel-bin-fold",
        action="store_true",
        help="Enable the scanner's opt-in Rosenbluth QE fold.",
    )
    parser.add_argument(
        "--qel-bin-width-energy",
        type=float,
        default=0.025,
        help="Finite Eprime/omega bin width in GeV for non-default fold densities.",
    )
    parser.add_argument(
        "--qel-bin-width-costh",
        type=float,
        default=0.08,
        help="Finite costheta_l bin width for non-default fold densities.",
    )
    parser.add_argument(
        "--qel-fold-method",
        choices=("lattice", "grid"),
        default="lattice",
        help="Deterministic integration method for --qel-bin-fold.",
    )
    parser.add_argument(
        "--qel-fold-density",
        choices=("generator-q2", "exact-theta", "q2-jacobian", "qel-delta"),
        default="exact-theta",
        help="Rosenbluth density used inside --qel-bin-fold.",
    )
    parser.add_argument(
        "--qel-fold-samples",
        type=int,
        default=200000,
        help="Advanced Rosenbluth fold sample count for non-default densities.",
    )
    parser.add_argument(
        "--no-qel-fold-scan-cache",
        action="store_true",
        help="Disable the scanner's one-pass scan fill for supported advanced densities.",
    )
    parser.add_argument(
        "--qel-fold-nr",
        type=int,
        default=12,
        help="Radius grid/CDF resolution for --qel-bin-fold.",
    )
    parser.add_argument(
        "--qel-fold-np",
        type=int,
        default=64,
        help="Momentum grid/CDF resolution for --qel-bin-fold.",
    )
    parser.add_argument(
        "--qel-fold-ncosth-p",
        type=int,
        default=16,
        help="Advanced nucleon-direction costheta samples.",
    )
    parser.add_argument(
        "--qel-fold-nphi-p",
        type=int,
        default=8,
        help="Nucleon-direction phi samples for exact-theta.",
    )
    parser.add_argument(
        "--qel-fold-ncos0",
        type=int,
        default=24,
        help="Advanced Rosenbluth fold sample count for non-default densities.",
    )
    parser.add_argument(
        "--qel-fold-nphi0",
        type=int,
        default=8,
        help="Advanced Rosenbluth fold azimuth count for non-default densities.",
    )
    parser.add_argument(
        "--qel-fold-kf",
        type=float,
        default=None,
        help="Override Fermi momentum in GeV for --qel-bin-fold.",
    )
    parser.add_argument(
        "--qel-fold-removal-energy",
        type=float,
        default=None,
        help="Override removal energy in GeV for --qel-bin-fold.",
    )
    parser.add_argument(
        "--shape-norm",
        choices=("off", "auto", "all", "empirical-mec"),
        default="empirical-mec",
        help="Forward scanner shape normalization mode.",
    )
    parser.add_argument(
        "--initial-state-fold",
        default="off",
        help="Manual processes to fold deterministically, or 'off' (default: off).",
    )
    parser.add_argument(
        "--initial-state-fold-samples",
        type=int,
        default=512,
        help="Low-discrepancy initial-state samples per differential point.",
    )
    parser.add_argument(
        "--initial-state-fold-nr",
        type=int,
        default=12,
        help="Radius CDF resolution for the RES/DIS fold.",
    )
    parser.add_argument(
        "--initial-state-fold-np",
        type=int,
        default=64,
        help="Momentum CDF resolution for RES/DIS and each MEC nucleon.",
    )
    parser.add_argument(
        "--initial-state-fold-event-phase-space",
        action="store_true",
        help=(
            "Apply GENIE's event-chain threshold and W/Q2 phase-space vetoes "
            "inside the deterministic fold."
        ),
    )
    parser.add_argument(
        "--initial-state-fold-debug",
        action="store_true",
        help="Forward per-component accepted-state diagnostics.",
    )
    parser.add_argument(
        "--shape-norm-ne",
        type=int,
        default=120,
        help="Eprime bins for --shape-norm.",
    )
    parser.add_argument(
        "--shape-norm-ncosth",
        type=int,
        default=120,
        help="costheta_l bins for --shape-norm.",
    )
    parser.add_argument(
        "--shape-norm-auto-threshold",
        type=float,
        default=2.0,
        help="Apply --shape-norm auto when Integral()/shape differs by this factor.",
    )
    parser.add_argument(
        "--empirical-mec-shape-norm",
        action="store_true",
        help="Compatibility alias for --shape-norm empirical-mec.",
    )
    parser.add_argument(
        "--empirical-mec-shape-norm-ne",
        type=int,
        default=None,
        help="Compatibility alias for --shape-norm-ne.",
    )
    parser.add_argument(
        "--empirical-mec-shape-norm-ncosth",
        type=int,
        default=None,
        help="Compatibility alias for --shape-norm-ncosth.",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if args.fold == "auto" and args.qel_bin_fold:
        raise SystemExit(
            "--qel-bin-fold is a manual override and cannot be combined with --fold auto"
        )
    if args.fold == "auto" and args.initial_state_fold.lower() not in ("off", "none"):
        raise SystemExit(
            "an explicit --initial-state-fold process list cannot be combined with --fold auto"
        )
    if args.empirical_mec_shape_norm:
        args.shape_norm = "empirical-mec"
    if args.empirical_mec_shape_norm_ne is not None:
        args.shape_norm_ne = args.empirical_mec_shape_norm_ne
    if args.empirical_mec_shape_norm_ncosth is not None:
        args.shape_norm_ncosth = args.empirical_mec_shape_norm_ncosth

    tunes = [x.strip() for x in args.tunes.split(",") if x.strip()]
    if len(tunes) != 2:
        raise SystemExit("--tunes must contain exactly two comma-separated tune IDs")
    if not SCANNER.exists():
        raise SystemExit(f"scanner is missing: {SCANNER}; run ./build.sh first")

    panels = PANEL_SETS[args.panel_set]
    args.outdir.mkdir(parents=True, exist_ok=True)
    scan_table = {}
    for tune in tunes:
        for panel in panels:
            if args.plot_only:
                csv_path = args.outdir / f"{compact_tune(tune)}_{panel['tag']}.csv"
                if not csv_path.exists():
                    raise SystemExit(f"missing cached CSV for --plot-only: {csv_path}")
            else:
                print(f"[scan] {tune} {panel['tag']}", flush=True)
                csv_path = run_scan(tune, panel, args.outdir, args)
            scan_table[(tune, panel["tag"])] = csv_path

    png, pdf = plot(scan_table, tunes, panels, args.outdir, args.stem, args)
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
