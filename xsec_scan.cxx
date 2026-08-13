// Standalone GENIE cross-section scanner.
//
// This file intentionally lives outside src/ so it can be copied next to a
// stock GENIE installation and compiled against that installation.

#include <TLorentzVector.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TVector3.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/ProcessInfo.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Interaction/XclsTag.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/RunOpt.h"
#if __has_include("Framework/Interaction/KPhaseSpaceCuts.h")
#include "Framework/Interaction/KPhaseSpaceCuts.h"
#define XSEC_SCAN_HAS_KPHASESPACECUTS 1
#else
#define XSEC_SCAN_HAS_KPHASESPACECUTS 0
#endif
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kDegToRad = kPi / 180.0;
constexpr double kDefaultRelTol = 1e-4;

std::string g_em_q2_min_status = "not-requested";
std::string g_em_q2_min_message;

class OutOfPhaseSpaceError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

enum class Var {
  E,
  Eprime,
  Q2,
  W,
  Omega,
  Q3,
  XB,
  Y,
  CosThetaL,
  ThetaLDeg,
  Tl,
};

enum class Observable {
  Total,
  D1,
  D2,
};

enum class Mode {
  Tune,
  Alg,
};

enum class NativePS {
  Auto,
  WQ2,
  TlCtl,
  XY,
  Q2,
  W,
};

enum class ShapeNormMode {
  Off,
  Auto,
  All,
  EmpiricalMEC,
};

struct RangeSpec {
  Var var;
  double min = 0.0;
  double max = 0.0;
  int n = 0;
};

struct Point {
  int row_id = 0;
  std::map<Var, double> values;
};

struct KinePoint {
  double E = 0.0;
  double Eprime = 0.0;
  double Q2 = 0.0;
  double W = 0.0;
  double omega = 0.0;
  double q3 = 0.0;
  double xB = 0.0;
  double y = 0.0;
  double costheta_l = 0.0;
  double theta_l_deg = 0.0;
  double Tl = 0.0;
};

struct ComponentValue {
  std::string component;
  double value = 0.0;
  std::string status = "ok";
  std::string message;
  std::map<Var, double> kin_values;
};

struct Options {
  Mode mode = Mode::Tune;
  std::string tune = "Default";
  std::string event_generator_list;
  std::string xsec_alg;
  std::string process = "DIS-EM";
  std::string intlist_alg;
  std::string intlist_config;
  std::string interaction_sum = "generated";
  int probe = genie::kPdgElectron;
  int target = genie::kPdgTgtFreeP;
  Observable observable = Observable::Total;
  std::vector<Var> diff_vars;
  std::vector<RangeSpec> scans;
  std::vector<RangeSpec> integrate_over;
  std::map<Var, double> fixed;
  std::string points_path;
  std::string output_path = "out/xsec_scan.csv";
  std::string xsec_unit = "cm2";
  NativePS native_ps = NativePS::Auto;
  double jac_step = 1e-5;
  double em_q2_min = 0.02;
  bool em_q2_min_set = false;
  bool qel_bin_fold = false;
  double qel_bin_width_energy = 0.005;
  double qel_bin_width_costh = 0.01;
  std::string qel_fold_method = "lattice";
  std::string qel_fold_density = "exact-theta";
  int qel_fold_samples = 200000;
  int qel_fold_nr = 12;
  int qel_fold_np = 64;
  int qel_fold_ncosth_p = 8;
  int qel_fold_nphi_p = 4;
  int qel_fold_ncos0 = 24;
  int qel_fold_nphi0 = 8;
  double qel_fold_pmax = -1.0;
  double qel_fold_kf = -1.0;
  double qel_fold_removal_energy = -1.0;
  bool qel_fold_scan_cache = true;
  bool qel_fold_debug = false;
  std::set<std::string> initial_state_fold;
  int initial_state_fold_samples = 512;
  int initial_state_fold_nr = 12;
  int initial_state_fold_np = 64;
  bool initial_state_fold_event_phase_space = false;
  bool initial_state_fold_debug = false;
  ShapeNormMode shape_norm = ShapeNormMode::Off;
  int shape_norm_ne = 120;
  int shape_norm_ncosth = 120;
  double shape_norm_auto_threshold = 2.0;
  bool components = false;
  bool strict = false;
  std::string message_thresholds = "config/Messenger_whisper.xml";
};

std::string trim(const std::string &s)
{
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = s.find_last_not_of(" \t\r\n");
  return s.substr(first, last - first + 1);
}

std::string lower(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

std::vector<std::string> split(const std::string &s, char delim)
{
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) out.push_back(trim(item));
  return out;
}

std::string csv_escape(const std::string &s)
{
  if (s.find_first_of(",\"\n\r") == std::string::npos) return s;
  std::string out = "\"";
  for (char c : s) {
    if (c == '"') out += "\"\"";
    else out += c;
  }
  out += "\"";
  return out;
}

std::string canonical_var_name(Var v)
{
  switch (v) {
    case Var::E: return "E";
    case Var::Eprime: return "Eprime";
    case Var::Q2: return "Q2";
    case Var::W: return "W";
    case Var::Omega: return "omega";
    case Var::Q3: return "q3";
    case Var::XB: return "xB";
    case Var::Y: return "y";
    case Var::CosThetaL: return "costheta_l";
    case Var::ThetaLDeg: return "theta_l_deg";
    case Var::Tl: return "Tl";
  }
  return "unknown";
}

std::string var_unit(Var v)
{
  switch (v) {
    case Var::E:
    case Var::Eprime:
    case Var::W:
    case Var::Omega:
    case Var::Q3:
    case Var::Tl:
      return "GeV";
    case Var::Q2:
      return "GeV2";
    case Var::ThetaLDeg:
      return "deg";
    case Var::XB:
    case Var::Y:
    case Var::CosThetaL:
      return "1";
  }
  return "1";
}

Var parse_var(std::string name)
{
  name = lower(trim(name));
  name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
  if (name == "e" || name == "ein" || name == "enu" || name == "ebeam") {
    return Var::E;
  }
  if (name == "eprime" || name == "ep" || name == "e'" || name == "eout" ||
      name == "efinal") {
    return Var::Eprime;
  }
  if (name == "q2" || name == "q^2") return Var::Q2;
  if (name == "w") return Var::W;
  if (name == "omega" || name == "q0" || name == "nu" || name == "v") {
    return Var::Omega;
  }
  if (name == "q3" || name == "absq" || name == "|q|" || name == "q") {
    return Var::Q3;
  }
  if (name == "xb" || name == "x_b" || name == "x" || name == "bjorkenx") {
    return Var::XB;
  }
  if (name == "y") return Var::Y;
  if (name == "costheta_l" || name == "costhetal" || name == "costhl" ||
      name == "ctl" || name == "cos_theta_l") {
    return Var::CosThetaL;
  }
  if (name == "theta_l_deg" || name == "thetadeg" || name == "theta_l" ||
      name == "theta") {
    return Var::ThetaLDeg;
  }
  if (name == "tl" || name == "t_l" || name == "leptonkinetic") {
    return Var::Tl;
  }
  throw std::runtime_error("Unknown variable name: " + name);
}

std::string observable_name(Observable obs)
{
  if (obs == Observable::Total) return "total";
  if (obs == Observable::D1) return "d1";
  return "d2";
}

std::string native_name(NativePS ps)
{
  switch (ps) {
    case NativePS::Auto: return "auto";
    case NativePS::WQ2: return "WQ2";
    case NativePS::TlCtl: return "Tlctl";
    case NativePS::XY: return "xy";
    case NativePS::Q2: return "Q2";
    case NativePS::W: return "W";
  }
  return "auto";
}

std::string shape_norm_name(ShapeNormMode mode)
{
  switch (mode) {
    case ShapeNormMode::Off: return "off";
    case ShapeNormMode::Auto: return "auto";
    case ShapeNormMode::All: return "all";
    case ShapeNormMode::EmpiricalMEC: return "empirical-mec";
  }
  return "off";
}

std::string diff_label(const std::vector<Var> &vars)
{
  std::string out;
  for (size_t i = 0; i < vars.size(); ++i) {
    if (i) out += ",";
    out += canonical_var_name(vars[i]);
  }
  return out;
}

std::string range_label(const std::vector<RangeSpec> &ranges)
{
  std::string out;
  for (size_t i = 0; i < ranges.size(); ++i) {
    if (i) out += ";";
    out += canonical_var_name(ranges[i].var) + ":" +
           std::to_string(ranges[i].min) + ":" +
           std::to_string(ranges[i].max) + ":" +
           std::to_string(ranges[i].n);
  }
  return out;
}

void print_usage(std::ostream &os)
{
  os
      << "GENIE standalone cross-section scanner\n\n"
      << "Examples:\n"
      << "  xsec_scan --mode tune --tune G18_01a_02_11a --probe 11 --target 1000010010 \\\n"
      << "    --observable d2 --diff W,Q2 --fixed E=2 --scan W:1.1:2.0:10 --scan Q2:0.1:1.0:10\n\n"
      << "  xsec_scan --mode alg --xsec-alg genie::KNOTunedQPMDISPXSec/Default \\\n"
      << "    --process DIS-EM --observable d2 --diff Eprime,costheta_l \\\n"
      << "    --fixed E=2 --scan Eprime:0.5:1.8:20 --scan costheta_l:0.5:0.95:10\n\n"
      << "Options:\n"
      << "  --mode tune|alg\n"
      << "  --tune NAME\n"
      << "  --event-generator-list NAME\n"
      << "  --xsec-alg genie::Class/Config\n"
      << "  --process QEL-EM|RES-EM|DIS-EM|MEC-EM|... for algorithm mode\n"
      << "  --probe PDG --target PDG\n"
      << "  --observable total|d1|d2\n"
      << "  --diff VAR or --diff VAR1,VAR2\n"
      << "  --fixed VAR=VALUE, repeatable\n"
      << "  --scan VAR:MIN:MAX:N, repeatable\n"
      << "  --integrate-over VAR:MIN:MAX:N, repeatable; requires MIN < MAX\n"
      << "  --points input.csv\n"
      << "  --output output.csv\n"
      << "  --xsec-unit cm2|nb|ub|pb\n"
      << "  --native-phase-space auto|WQ2|Tlctl|xy|Q2|W\n"
      << "  --em-q2-min VALUE opt-in EM phase-space override, if supported by GENIE\n"
      << "  --qel-bin-fold opt-in Rosenbluth QE fold for d2 Eprime/omega,costheta_l\n"
      << "  --qel-bin-width Eprime|omega|costheta_l=WIDTH\n"
      << "  --qel-fold-method lattice|grid\n"
      << "  --qel-fold-density exact-theta (default; advanced alternatives accepted)\n"
      << "  --qel-fold-samples N for the lattice method\n"
      << "  --qel-fold-nr N --qel-fold-np N --qel-fold-ncosth-p N --qel-fold-nphi-p N\n"
      << "  --qel-fold-ncos0 N --qel-fold-nphi0 N\n"
      << "  --qel-fold-pmax VALUE --qel-fold-kf VALUE --qel-fold-removal-energy VALUE\n"
      << "  --no-qel-fold-scan-cache disable one-pass scan filling for supported QE scans\n"
      << "  --initial-state-fold RES,MEC,DIS deterministic event-chain initial-state fold\n"
      << "  --initial-state-fold-samples N --initial-state-fold-nr N --initial-state-fold-np N\n"
      << "  --initial-state-fold-event-phase-space apply GENIE event-chain phase-space vetoes\n"
      << "  --initial-state-fold-debug print accepted-state diagnostics\n"
      << "  --shape-norm off|auto|all|empirical-mec rescale d2 shapes to Integral()\n"
      << "  --shape-norm-ne N --shape-norm-ncosth N --shape-norm-auto-threshold FACTOR\n"
      << "  --empirical-mec-shape-norm compatibility alias for --shape-norm empirical-mec\n"
      << "  --components\n"
      << "  --strict\n";
}

RangeSpec parse_range(const std::string &s)
{
  const auto fields = split(s, ':');
  if (fields.size() != 4) {
    throw std::runtime_error("Range must have form VAR:MIN:MAX:N: " + s);
  }
  RangeSpec r;
  r.var = parse_var(fields[0]);
  r.min = std::stod(fields[1]);
  r.max = std::stod(fields[2]);
  r.n = std::stoi(fields[3]);
  if (!std::isfinite(r.min) || !std::isfinite(r.max)) {
    throw std::runtime_error("Range endpoints must be finite: " + s);
  }
  if (r.n <= 0) throw std::runtime_error("Range N must be positive: " + s);
  return r;
}

std::pair<Var, double> parse_fixed(const std::string &s)
{
  const auto pos = s.find('=');
  if (pos == std::string::npos) {
    throw std::runtime_error("Fixed value must have form VAR=VALUE: " + s);
  }
  const double value = std::stod(s.substr(pos + 1));
  if (!std::isfinite(value)) {
    throw std::runtime_error("Fixed value must be finite: " + s);
  }
  return {parse_var(s.substr(0, pos)), value};
}

Mode parse_mode(const std::string &s)
{
  const std::string v = lower(s);
  if (v == "tune") return Mode::Tune;
  if (v == "alg" || v == "algorithm") return Mode::Alg;
  throw std::runtime_error("Unknown mode: " + s);
}

Observable parse_observable(const std::string &s)
{
  const std::string v = lower(s);
  if (v == "total" || v == "sigma") return Observable::Total;
  if (v == "d1" || v == "dsigma") return Observable::D1;
  if (v == "d2" || v == "d2sigma") return Observable::D2;
  throw std::runtime_error("Unknown observable: " + s);
}

NativePS parse_native_ps(const std::string &s)
{
  const std::string v = lower(s);
  if (v == "auto") return NativePS::Auto;
  if (v == "wq2" || v == "w,q2") return NativePS::WQ2;
  if (v == "tlctl" || v == "tl,costheta_l" || v == "tl,costhl") return NativePS::TlCtl;
  if (v == "xy" || v == "x,y") return NativePS::XY;
  if (v == "q2") return NativePS::Q2;
  if (v == "w") return NativePS::W;
  throw std::runtime_error("Unknown native phase space: " + s);
}

ShapeNormMode parse_shape_norm_mode(const std::string &s)
{
  const std::string v = lower(s);
  if (v == "off" || v == "none" || v == "false" || v == "0") {
    return ShapeNormMode::Off;
  }
  if (v == "auto") return ShapeNormMode::Auto;
  if (v == "all" || v == "on" || v == "true" || v == "1") {
    return ShapeNormMode::All;
  }
  if (v == "empirical-mec" || v == "empirical_mec" ||
      v == "empiricalmec" || v == "mec") {
    return ShapeNormMode::EmpiricalMEC;
  }
  throw std::runtime_error("Unknown shape normalization mode: " + s);
}

Options parse_options(int argc, char **argv)
{
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto need_value = [&](const std::string &name) -> std::string {
      if (i + 1 >= argc) throw std::runtime_error("Missing value for " + name);
      return argv[++i];
    };
    if (arg == "--help" || arg == "-h") {
      print_usage(std::cout);
      std::exit(0);
    } else if (arg == "--mode") {
      opt.mode = parse_mode(need_value(arg));
    } else if (arg == "--tune") {
      opt.tune = need_value(arg);
    } else if (arg == "--event-generator-list") {
      opt.event_generator_list = need_value(arg);
    } else if (arg == "--xsec-alg") {
      opt.xsec_alg = need_value(arg);
    } else if (arg == "--process") {
      opt.process = need_value(arg);
    } else if (arg == "--interaction-list-generator") {
      const std::string spec = need_value(arg);
      const auto slash = spec.find('/');
      opt.intlist_alg = spec.substr(0, slash);
      opt.intlist_config = (slash == std::string::npos) ? "Default" : spec.substr(slash + 1);
    } else if (arg == "--interaction-sum") {
      opt.interaction_sum = lower(need_value(arg));
    } else if (arg == "--probe") {
      opt.probe = std::stoi(need_value(arg));
    } else if (arg == "--target") {
      opt.target = std::stoi(need_value(arg));
    } else if (arg == "--observable") {
      opt.observable = parse_observable(need_value(arg));
    } else if (arg == "--diff") {
      opt.diff_vars.clear();
      for (const std::string &v : split(need_value(arg), ',')) {
        opt.diff_vars.push_back(parse_var(v));
      }
    } else if (arg == "--fixed") {
      const auto fixed = parse_fixed(need_value(arg));
      opt.fixed[fixed.first] = fixed.second;
    } else if (arg == "--scan") {
      opt.scans.push_back(parse_range(need_value(arg)));
    } else if (arg == "--integrate-over") {
      opt.integrate_over.push_back(parse_range(need_value(arg)));
    } else if (arg == "--points") {
      opt.points_path = need_value(arg);
    } else if (arg == "--output") {
      opt.output_path = need_value(arg);
    } else if (arg == "--xsec-unit") {
      opt.xsec_unit = lower(need_value(arg));
    } else if (arg == "--native-phase-space") {
      opt.native_ps = parse_native_ps(need_value(arg));
    } else if (arg == "--jac-step") {
      opt.jac_step = std::stod(need_value(arg));
    } else if (arg == "--em-q2-min") {
      opt.em_q2_min = std::stod(need_value(arg));
      opt.em_q2_min_set = true;
    } else if (arg == "--qel-bin-fold") {
      opt.qel_bin_fold = true;
    } else if (arg == "--qel-bin-width") {
      const auto width = parse_fixed(need_value(arg));
      if (width.second <= 0.0) {
        throw std::runtime_error("--qel-bin-width must be positive");
      }
      if (width.first == Var::Eprime || width.first == Var::Omega) {
        opt.qel_bin_width_energy = width.second;
      } else if (width.first == Var::CosThetaL) {
        opt.qel_bin_width_costh = width.second;
      } else {
        throw std::runtime_error("--qel-bin-width supports Eprime, omega, or costheta_l");
      }
    } else if (arg == "--qel-fold-method") {
      opt.qel_fold_method = lower(need_value(arg));
    } else if (arg == "--qel-fold-density") {
      opt.qel_fold_density = lower(need_value(arg));
    } else if (arg == "--qel-fold-samples") {
      opt.qel_fold_samples = std::stoi(need_value(arg));
    } else if (arg == "--qel-fold-np") {
      opt.qel_fold_np = std::stoi(need_value(arg));
    } else if (arg == "--qel-fold-nr") {
      opt.qel_fold_nr = std::stoi(need_value(arg));
    } else if (arg == "--qel-fold-ncosth-p") {
      opt.qel_fold_ncosth_p = std::stoi(need_value(arg));
    } else if (arg == "--qel-fold-nphi-p") {
      opt.qel_fold_nphi_p = std::stoi(need_value(arg));
    } else if (arg == "--qel-fold-ncos0") {
      opt.qel_fold_ncos0 = std::stoi(need_value(arg));
    } else if (arg == "--qel-fold-nphi0") {
      opt.qel_fold_nphi0 = std::stoi(need_value(arg));
    } else if (arg == "--qel-fold-nangle") {
      opt.qel_fold_ncos0 = std::stoi(need_value(arg));
    } else if (arg == "--qel-fold-pmax") {
      opt.qel_fold_pmax = std::stod(need_value(arg));
    } else if (arg == "--qel-fold-kf") {
      opt.qel_fold_kf = std::stod(need_value(arg));
    } else if (arg == "--qel-fold-removal-energy") {
      opt.qel_fold_removal_energy = std::stod(need_value(arg));
    } else if (arg == "--qel-fold-scan-cache") {
      opt.qel_fold_scan_cache = true;
    } else if (arg == "--no-qel-fold-scan-cache") {
      opt.qel_fold_scan_cache = false;
    } else if (arg == "--qel-fold-debug") {
      opt.qel_fold_debug = true;
    } else if (arg == "--initial-state-fold") {
      for (std::string process : split(need_value(arg), ',')) {
        process = lower(process);
        if (process == "all") {
          opt.initial_state_fold.insert("res");
          opt.initial_state_fold.insert("mec");
          opt.initial_state_fold.insert("dis");
        } else {
          opt.initial_state_fold.insert(process);
        }
      }
    } else if (arg == "--initial-state-fold-samples") {
      opt.initial_state_fold_samples = std::stoi(need_value(arg));
    } else if (arg == "--initial-state-fold-nr") {
      opt.initial_state_fold_nr = std::stoi(need_value(arg));
    } else if (arg == "--initial-state-fold-np") {
      opt.initial_state_fold_np = std::stoi(need_value(arg));
    } else if (arg == "--initial-state-fold-event-phase-space") {
      opt.initial_state_fold_event_phase_space = true;
    } else if (arg == "--initial-state-fold-debug") {
      opt.initial_state_fold_debug = true;
    } else if (arg == "--shape-norm") {
      opt.shape_norm = parse_shape_norm_mode(need_value(arg));
    } else if (arg == "--shape-norm-ne") {
      opt.shape_norm_ne = std::stoi(need_value(arg));
    } else if (arg == "--shape-norm-ncosth") {
      opt.shape_norm_ncosth = std::stoi(need_value(arg));
    } else if (arg == "--shape-norm-auto-threshold") {
      opt.shape_norm_auto_threshold = std::stod(need_value(arg));
    } else if (arg == "--empirical-mec-shape-norm") {
      opt.shape_norm = ShapeNormMode::EmpiricalMEC;
    } else if (arg == "--empirical-mec-shape-norm-ne") {
      opt.shape_norm_ne = std::stoi(need_value(arg));
    } else if (arg == "--empirical-mec-shape-norm-ncosth") {
      opt.shape_norm_ncosth = std::stoi(need_value(arg));
    } else if (arg == "--message-thresholds") {
      opt.message_thresholds = need_value(arg);
    } else if (arg == "--components") {
      opt.components = true;
    } else if (arg == "--strict") {
      opt.strict = true;
    } else {
      throw std::runtime_error("Unknown option: " + arg);
    }
  }

  if (opt.event_generator_list.empty()) {
    opt.event_generator_list =
        (opt.probe == genie::kPdgElectron || opt.probe == genie::kPdgPositron)
            ? "EM"
            : "Default";
  }
  if (opt.observable == Observable::D1 && opt.diff_vars.size() != 1) {
    throw std::runtime_error("--observable d1 requires exactly one --diff variable");
  }
  if (opt.observable == Observable::D2 && opt.diff_vars.size() != 2) {
    throw std::runtime_error("--observable d2 requires exactly two --diff variables");
  }
  if (opt.observable == Observable::Total && !opt.diff_vars.empty()) {
    throw std::runtime_error("--observable total does not take --diff");
  }
  if (opt.observable == Observable::D2 && !opt.integrate_over.empty()) {
    throw std::runtime_error("--observable d2 does not take --integrate-over");
  }
  if (opt.observable == Observable::D1 && opt.integrate_over.size() > 1) {
    throw std::runtime_error("--observable d1 supports at most one --integrate-over range");
  }
  if (opt.observable == Observable::Total && opt.integrate_over.size() > 2) {
    throw std::runtime_error("--observable total supports at most two --integrate-over ranges");
  }
  std::set<Var> seen_diff;
  for (Var var : opt.diff_vars) {
    if (var == Var::E) {
      throw std::runtime_error("E is a conditioning variable and cannot be used with --diff");
    }
    if (!seen_diff.insert(var).second) {
      throw std::runtime_error("Duplicate --diff variable: " + canonical_var_name(var));
    }
  }
  std::set<Var> seen_scans;
  for (const RangeSpec &scan : opt.scans) {
    if (!seen_scans.insert(scan.var).second) {
      throw std::runtime_error("Duplicate --scan variable: " + canonical_var_name(scan.var));
    }
    if (opt.fixed.count(scan.var)) {
      throw std::runtime_error("Variable cannot be both --fixed and --scan: " +
                               canonical_var_name(scan.var));
    }
  }
  std::set<Var> seen_integrals;
  for (const RangeSpec &range : opt.integrate_over) {
    if (range.var == Var::E) {
      throw std::runtime_error("E is a conditioning variable and cannot be integrated over");
    }
    if (!(range.max > range.min)) {
      throw std::runtime_error("--integrate-over requires MIN < MAX for " +
                               canonical_var_name(range.var));
    }
    if (!seen_integrals.insert(range.var).second) {
      throw std::runtime_error("Duplicate --integrate-over variable: " +
                               canonical_var_name(range.var));
    }
    if (seen_diff.count(range.var)) {
      throw std::runtime_error("Variable cannot be both --diff and --integrate-over: " +
                               canonical_var_name(range.var));
    }
    if (seen_scans.count(range.var) || opt.fixed.count(range.var)) {
      throw std::runtime_error("Integrated variable must not also be supplied as a point value: " +
                               canonical_var_name(range.var));
    }
  }
  if (opt.mode == Mode::Alg && opt.xsec_alg.empty()) {
    throw std::runtime_error("--mode alg requires --xsec-alg");
  }
  if (opt.interaction_sum != "generated" && opt.interaction_sum != "hit-state") {
    throw std::runtime_error("--interaction-sum must be generated or hit-state");
  }
  if (!std::isfinite(opt.jac_step) || opt.jac_step <= 0.0) {
    throw std::runtime_error("--jac-step must be positive");
  }
  if (opt.em_q2_min_set &&
      (!std::isfinite(opt.em_q2_min) || opt.em_q2_min < 0.0)) {
    throw std::runtime_error("--em-q2-min must be finite and non-negative");
  }
  if (opt.qel_fold_np <= 0) {
    throw std::runtime_error("--qel-fold-np must be positive");
  }
  if (opt.qel_fold_method != "lattice" && opt.qel_fold_method != "grid") {
    throw std::runtime_error("--qel-fold-method must be lattice or grid");
  }
  if (opt.qel_fold_density != "generator-q2" &&
      opt.qel_fold_density != "exact-theta" &&
      opt.qel_fold_density != "q2-jacobian" &&
      opt.qel_fold_density != "qel-delta") {
    throw std::runtime_error(
        "--qel-fold-density must be generator-q2, exact-theta, q2-jacobian, or qel-delta");
  }
  if (opt.qel_fold_samples <= 0) {
    throw std::runtime_error("--qel-fold-samples must be positive");
  }
  if (opt.qel_fold_nr <= 0) {
    throw std::runtime_error("--qel-fold-nr must be positive");
  }
  if (opt.qel_fold_ncosth_p <= 0) {
    throw std::runtime_error("--qel-fold-ncosth-p must be positive");
  }
  if (opt.qel_fold_nphi_p <= 0) {
    throw std::runtime_error("--qel-fold-nphi-p must be positive");
  }
  if (opt.qel_fold_ncos0 <= 0) {
    throw std::runtime_error("--qel-fold-ncos0 must be positive");
  }
  if (opt.qel_fold_nphi0 <= 0) {
    throw std::runtime_error("--qel-fold-nphi0 must be positive");
  }
  if (opt.qel_fold_pmax == 0.0 || opt.qel_fold_pmax < -1.0) {
    throw std::runtime_error("--qel-fold-pmax must be positive or left unset");
  }
  if (!std::isfinite(opt.qel_fold_pmax) ||
      !std::isfinite(opt.qel_fold_kf) ||
      !std::isfinite(opt.qel_fold_removal_energy)) {
    throw std::runtime_error("QE fold momentum and energy settings must be finite");
  }
  for (const std::string &process : opt.initial_state_fold) {
    if (process != "res" && process != "mec" && process != "dis") {
      throw std::runtime_error(
          "--initial-state-fold accepts only RES, MEC, DIS, or all");
    }
  }
  if (opt.initial_state_fold_samples <= 0 ||
      opt.initial_state_fold_nr <= 0 || opt.initial_state_fold_np <= 0) {
    throw std::runtime_error("Initial-state fold grid sizes must be positive");
  }
  if (opt.shape_norm_ne <= 0 || opt.shape_norm_ncosth <= 0) {
    throw std::runtime_error("--shape-norm grid sizes must be positive");
  }
  if (!std::isfinite(opt.shape_norm_auto_threshold) ||
      opt.shape_norm_auto_threshold <= 1.0) {
    throw std::runtime_error("--shape-norm-auto-threshold must be greater than 1");
  }
  return opt;
}

void prepend_existing_path(std::vector<std::filesystem::path> &paths,
                           const std::filesystem::path &path)
{
  std::error_code ec;
  if (std::filesystem::exists(path, ec) && std::filesystem::is_directory(path, ec)) {
    paths.push_back(std::filesystem::weakly_canonical(path, ec));
  }
}

void configure_runtime_xml_path(const char *argv0)
{
  std::vector<std::filesystem::path> paths;
  std::error_code ec;
  const std::filesystem::path exe =
      std::filesystem::weakly_canonical(std::filesystem::path(argv0), ec);
  const std::filesystem::path exe_dir = exe.parent_path();

  prepend_existing_path(paths, exe_dir.parent_path() / "config");
  prepend_existing_path(paths, std::filesystem::current_path(ec) / "config");

  std::string combined;
  std::set<std::string> seen;
  for (const auto &path : paths) {
    const std::string s = path.string();
    if (s.empty() || seen.count(s)) continue;
    if (!combined.empty()) combined += ":";
    combined += s;
    seen.insert(s);
  }
  if (const char *old = std::getenv("GXMLPATH")) {
    std::set<std::string> default_paths;
    if (const char *genie = std::getenv("GENIE")) {
      const std::filesystem::path genie_config =
          std::filesystem::weakly_canonical(std::filesystem::path(genie) / "config", ec);
      default_paths.insert(genie_config.string());
    }

    std::string token;
    std::stringstream ss(old);
    while (std::getline(ss, token, ':')) {
      std::stringstream ss2(token);
      std::string subtoken;
      while (std::getline(ss2, subtoken, ';')) {
        std::stringstream ss3(subtoken);
        std::string path_token;
        while (std::getline(ss3, path_token, ',')) {
          if (path_token.empty()) continue;
          const std::filesystem::path old_path =
              std::filesystem::weakly_canonical(std::filesystem::path(path_token), ec);
          const std::string canonical = old_path.string();
          if (!canonical.empty() && default_paths.count(canonical)) continue;
          const std::string keep = canonical.empty() ? path_token : canonical;
          if (seen.count(keep)) continue;
          if (!combined.empty()) combined += ":";
          combined += keep;
          seen.insert(keep);
        }
      }
    }
  }
  if (!combined.empty()) setenv("GXMLPATH", combined.c_str(), 1);
}

bool set_registry_double(const std::string &file_id, const std::string &set_name,
                         const std::string &param, double value)
{
  genie::Registry *r =
      genie::AlgConfigPool::Instance()->CommonList(file_id, set_name);
  if (!r || !r->Exists(param)) return false;
  r->UnLock();
  r->Set(param, value);
  r->Lock();
  return true;
}

void apply_em_q2_min_override(const Options &opt)
{
  g_em_q2_min_status = "not-requested";
  g_em_q2_min_message.clear();
  if (!opt.em_q2_min_set) return;

  const bool kinematics_applied =
      set_registry_double("Param", "Kinematics", "EM-Q2-min", opt.em_q2_min);

#if XSEC_SCAN_HAS_KPHASESPACECUTS
  genie::KPhaseSpaceCuts::Instance()->SetQ2MinOverride(opt.em_q2_min);
  const bool phase_space_applied =
      set_registry_double("PhaseSpaceCuts", "Default", "EM-Q2-min", opt.em_q2_min);

  g_em_q2_min_status = "applied";
  std::string targets = "KPhaseSpaceCuts";
  if (phase_space_applied) targets += ";PhaseSpaceCuts/Default";
  if (kinematics_applied) targets += ";Param/Kinematics";
  g_em_q2_min_message = "set EM-Q2-min through " + targets;
  return;
#else
  if (kinematics_applied) {
    g_em_q2_min_status = "applied";
    g_em_q2_min_message = "set EM-Q2-min through Param/Kinematics";
    return;
  }

  g_em_q2_min_status = "unavailable";
  g_em_q2_min_message =
      "EM-Q2-min override requested, but this GENIE install was built without "
      "KPhaseSpaceCuts and has no Param/Kinematics EM-Q2-min entry; using the "
      "linked GENIE built-in/default phase-space threshold";
  std::cerr << "warning: " << g_em_q2_min_message << "\n";
  if (opt.strict) throw std::runtime_error(g_em_q2_min_message);
#endif
}

double xsec_unit_factor(const std::string &unit)
{
  if (unit == "cm2") return genie::units::cm2;
  if (unit == "nb" || unit == "nanobarn") return genie::units::nanobarn;
  if (unit == "ub" || unit == "microbarn") return genie::units::microbarn;
  if (unit == "pb" || unit == "picobarn") return genie::units::picobarn;
  throw std::runtime_error("Unknown cross-section unit: " + unit);
}

std::string value_unit_label(const Options &opt)
{
  std::string unit = opt.xsec_unit;
  if (opt.observable == Observable::D1) {
    unit += "/" + var_unit(opt.diff_vars.at(0));
  } else if (opt.observable == Observable::D2) {
    unit += "/" + var_unit(opt.diff_vars.at(0)) + "/" + var_unit(opt.diff_vars.at(1));
  }
  return unit;
}

bool has(const std::map<Var, double> &m, Var v)
{
  return m.find(v) != m.end();
}

double get(const std::map<Var, double> &m, Var v)
{
  const auto it = m.find(v);
  if (it == m.end()) throw std::runtime_error("Missing " + canonical_var_name(v));
  return it->second;
}

void add_scan_points_recursive(const std::vector<RangeSpec> &ranges, size_t idx,
                               std::map<Var, double> current,
                               std::vector<Point> &points)
{
  if (idx == ranges.size()) {
    Point p;
    p.row_id = static_cast<int>(points.size());
    p.values = std::move(current);
    points.push_back(std::move(p));
    return;
  }
  const RangeSpec &r = ranges[idx];
  for (int i = 0; i < r.n; ++i) {
    const double t = (r.n == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(r.n - 1);
    current[r.var] = r.min + t * (r.max - r.min);
    add_scan_points_recursive(ranges, idx + 1, current, points);
  }
}

std::vector<Point> read_points_csv(const std::string &path,
                                   const std::map<Var, double> &fixed)
{
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Could not open points CSV: " + path);

  std::string line;
  std::vector<Var> columns;
  std::set<Var> seen_columns;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;
    for (const std::string &name : split(line, ',')) {
      const Var var = parse_var(name);
      if (!seen_columns.insert(var).second) {
        throw std::runtime_error("Duplicate points CSV column: " +
                                 canonical_var_name(var));
      }
      columns.push_back(var);
    }
    break;
  }
  if (columns.empty()) throw std::runtime_error("No header found in points CSV: " + path);

  std::vector<Point> points;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;
    const auto fields = split(line, ',');
    if (fields.size() != columns.size()) {
      throw std::runtime_error("CSV row has wrong number of columns: " + line);
    }
    Point p;
    p.row_id = static_cast<int>(points.size());
    p.values = fixed;
    for (size_t i = 0; i < fields.size(); ++i) {
      const double value = std::stod(fields[i]);
      if (!std::isfinite(value)) {
        throw std::runtime_error("Points CSV values must be finite in row: " + line);
      }
      p.values[columns[i]] = value;
    }
    points.push_back(std::move(p));
  }
  if (points.empty()) {
    throw std::runtime_error("No data rows found in points CSV: " + path);
  }
  return points;
}

std::vector<Point> make_points(const Options &opt)
{
  if (!opt.points_path.empty() && !opt.scans.empty()) {
    throw std::runtime_error("Use either --points or --scan, not both");
  }
  if (!opt.points_path.empty()) return read_points_csv(opt.points_path, opt.fixed);

  std::vector<Point> points;
  add_scan_points_recursive(opt.scans, 0, opt.fixed, points);
  if (points.empty()) {
    Point p;
    p.row_id = 0;
    p.values = opt.fixed;
    points.push_back(std::move(p));
  }
  return points;
}

void validate_projection_point_values(const std::vector<Point> &points,
                                      const Options &opt)
{
  for (const Point &point : points) {
    for (const RangeSpec &range : opt.integrate_over) {
      if (has(point.values, range.var)) {
        throw std::runtime_error(
            "Integrated variable must not also be present in point row " +
            std::to_string(point.row_id) + ": " +
            canonical_var_name(range.var));
      }
    }
  }
}

double particle_mass(int pdg)
{
  TParticlePDG *p = genie::PDGLibrary::Instance()->Find(pdg, false);
  return p ? p->Mass() : 0.0;
}

double final_lepton_mass(const genie::Interaction &interaction)
{
  TParticlePDG *p = interaction.FSPrimLepton();
  if (p) return p->Mass();
  return particle_mass(interaction.InitState().ProbePdg());
}

double hit_mass(const genie::Interaction &interaction)
{
  if (interaction.InitState().Tgt().HitNucIsSet()) {
    return interaction.InitState().Tgt().HitNucMass();
  }
  return interaction.InitState().Tgt().Mass();
}

double sqr(double x) { return x * x; }

double clamp(double x, double lo, double hi)
{
  return std::max(lo, std::min(hi, x));
}

void check_close(const std::string &name, double supplied, double derived,
                 double reltol = kDefaultRelTol)
{
  const double scale = std::max(1.0, std::max(std::fabs(supplied), std::fabs(derived)));
  if (std::fabs(supplied - derived) > reltol * scale) {
    std::ostringstream ss;
    ss << "Inconsistent " << name << ": supplied " << supplied
       << ", derived " << derived;
    throw std::runtime_error(ss.str());
  }
}

std::optional<double> solve_ep_bisection(double lo, double hi,
                                         const std::function<double(double)> &f)
{
  const int scan_bins = 200;
  double a = lo;
  double fa = f(a);
  for (int i = 1; i <= scan_bins; ++i) {
    const double b = lo + (hi - lo) * static_cast<double>(i) / scan_bins;
    const double fb = f(b);
    if (!std::isfinite(fa)) {
      a = b;
      fa = fb;
      continue;
    }
    if (std::isfinite(fb) && fa * fb <= 0.0) {
      double left = a;
      double right = b;
      double fleft = fa;
      for (int iter = 0; iter < 100; ++iter) {
        const double mid = 0.5 * (left + right);
        const double fmid = f(mid);
        if (!std::isfinite(fmid)) break;
        if (std::fabs(fmid) < 1e-12 || std::fabs(right - left) < 1e-10) {
          return mid;
        }
        if (fleft * fmid <= 0.0) {
          right = mid;
        } else {
          left = mid;
          fleft = fmid;
        }
      }
      return 0.5 * (left + right);
    }
    a = b;
    fa = fb;
  }
  return std::nullopt;
}

KinePoint solve_kinematics(const std::map<Var, double> &values,
                           const genie::Interaction &interaction)
{
  if (!has(values, Var::E)) {
    throw std::runtime_error("Kinematics are underconstrained: E is required");
  }

  const double E = get(values, Var::E);
  const double mi = particle_mass(interaction.InitState().ProbePdg());
  const double mf = final_lepton_mass(interaction);
  const double M = hit_mass(interaction);
  if (E <= mi) throw OutOfPhaseSpaceError("Incoming E is below probe mass");

  const double pi = std::sqrt(std::max(0.0, E * E - mi * mi));
  std::vector<std::pair<std::string, double>> ep_candidates;
  auto add_ep = [&](const std::string &why, double ep) {
    if (std::isfinite(ep)) ep_candidates.push_back({why, ep});
  };

  if (has(values, Var::Eprime)) add_ep("Eprime", get(values, Var::Eprime));
  if (has(values, Var::Tl)) add_ep("Tl", get(values, Var::Tl) + mf);
  if (has(values, Var::Omega)) add_ep("omega", E - get(values, Var::Omega));
  if (has(values, Var::Y)) add_ep("y", E * (1.0 - get(values, Var::Y)));
  if (has(values, Var::W) && has(values, Var::Q2)) {
    const double omega = (sqr(get(values, Var::W)) + get(values, Var::Q2) - M * M) / (2.0 * M);
    add_ep("W,Q2", E - omega);
  }
  if (has(values, Var::XB) && has(values, Var::Q2)) {
    const double x = get(values, Var::XB);
    if (x != 0.0) add_ep("xB,Q2", E - get(values, Var::Q2) / (2.0 * M * x));
  }
  if (has(values, Var::Q3) && has(values, Var::Q2)) {
    const double rad = sqr(get(values, Var::Q3)) - get(values, Var::Q2);
    if (rad < 0.0) {
      throw OutOfPhaseSpaceError("Supplied q3 and Q2 have no physical energy transfer");
    }
    add_ep("q3,Q2", E - std::sqrt(rad));
  }
  if (has(values, Var::W) && has(values, Var::Q3)) {
    const double omega = -M + std::sqrt(std::max(0.0, sqr(get(values, Var::W)) + sqr(get(values, Var::Q3))));
    add_ep("W,q3", E - omega);
  }
  if (has(values, Var::XB) && has(values, Var::Q3)) {
    const double x = get(values, Var::XB);
    const double omega = -M * x + std::sqrt(std::max(0.0, sqr(M * x) + sqr(get(values, Var::Q3))));
    add_ep("xB,q3", E - omega);
  }

  auto q2_from_angle = [&](double ep, double ctl) {
    if (ep < mf) return std::numeric_limits<double>::quiet_NaN();
    const double pf = std::sqrt(std::max(0.0, ep * ep - mf * mf));
    return -mi * mi - mf * mf + 2.0 * E * ep - 2.0 * pi * pf * ctl;
  };

  bool attempted_ep_solution = false;
  if (has(values, Var::CosThetaL) || has(values, Var::ThetaLDeg)) {
    if (has(values, Var::ThetaLDeg)) {
      const double theta = get(values, Var::ThetaLDeg);
      if (theta < 0.0 || theta > 180.0) {
        throw OutOfPhaseSpaceError("Supplied theta_l_deg is outside [0,180]");
      }
    }
    const double ctl = has(values, Var::CosThetaL)
                           ? get(values, Var::CosThetaL)
                           : std::cos(get(values, Var::ThetaLDeg) * kDegToRad);
    if (ctl < -1.0 || ctl > 1.0) {
      throw OutOfPhaseSpaceError("Supplied costheta_l is outside [-1,1]");
    }
    if (has(values, Var::Q2)) {
      attempted_ep_solution = true;
      auto ep = solve_ep_bisection(mf, E, [&](double x) {
        return q2_from_angle(x, ctl) - get(values, Var::Q2);
      });
      if (ep) add_ep("costheta_l,Q2", *ep);
    }
    if (has(values, Var::W)) {
      attempted_ep_solution = true;
      auto ep = solve_ep_bisection(mf, E, [&](double x) {
        const double q2 = q2_from_angle(x, ctl);
        return M * M + 2.0 * M * (E - x) - q2 - sqr(get(values, Var::W));
      });
      if (ep) add_ep("costheta_l,W", *ep);
    }
    if (has(values, Var::XB)) {
      attempted_ep_solution = true;
      auto ep = solve_ep_bisection(mf, E, [&](double x) {
        const double q2 = q2_from_angle(x, ctl);
        return q2 - 2.0 * M * get(values, Var::XB) * (E - x);
      });
      if (ep) add_ep("costheta_l,xB", *ep);
    }
    if (has(values, Var::Q3)) {
      attempted_ep_solution = true;
      auto ep = solve_ep_bisection(mf, E, [&](double x) {
        const double q2 = q2_from_angle(x, ctl);
        return sqr(E - x) + q2 - sqr(get(values, Var::Q3));
      });
      if (ep) add_ep("costheta_l,q3", *ep);
    }
  }

  if (ep_candidates.empty()) {
    if (attempted_ep_solution) {
      throw OutOfPhaseSpaceError("Supplied kinematics have no physical Eprime solution");
    }
    throw std::runtime_error(
        "Kinematics are underconstrained: supply E plus enough of Eprime, omega, y, W, Q2, xB, q3, costheta_l");
  }

  double Eprime = ep_candidates.front().second;
  for (const auto &candidate : ep_candidates) {
    check_close("Eprime from " + candidate.first, candidate.second, Eprime);
  }
  if (Eprime < mf || Eprime > E + kDefaultRelTol * std::max(1.0, E)) {
    throw OutOfPhaseSpaceError("Derived Eprime is outside the physical range");
  }

  const double omega = E - Eprime;
  std::vector<std::pair<std::string, double>> q2_candidates;
  auto add_q2 = [&](const std::string &why, double q2) {
    if (std::isfinite(q2)) q2_candidates.push_back({why, q2});
  };
  if (has(values, Var::Q2)) add_q2("Q2", get(values, Var::Q2));
  if (has(values, Var::CosThetaL) || has(values, Var::ThetaLDeg)) {
    const double ctl = has(values, Var::CosThetaL)
                           ? get(values, Var::CosThetaL)
                           : std::cos(get(values, Var::ThetaLDeg) * kDegToRad);
    add_q2("costheta_l", q2_from_angle(Eprime, ctl));
  }
  if (has(values, Var::W)) add_q2("W", M * M + 2.0 * M * omega - sqr(get(values, Var::W)));
  if (has(values, Var::XB)) add_q2("xB", 2.0 * M * omega * get(values, Var::XB));
  if (has(values, Var::Q3)) add_q2("q3", sqr(get(values, Var::Q3)) - omega * omega);

  if (q2_candidates.empty()) {
    throw std::runtime_error("Kinematics are underconstrained: could not derive Q2");
  }
  double Q2 = q2_candidates.front().second;
  for (const auto &candidate : q2_candidates) {
    check_close("Q2 from " + candidate.first, candidate.second, Q2);
  }
  if (Q2 < -kDefaultRelTol) throw OutOfPhaseSpaceError("Derived Q2 is negative");
  Q2 = std::max(0.0, Q2);

  const double W2 = M * M + 2.0 * M * omega - Q2;
  if (W2 <= 0.0) throw OutOfPhaseSpaceError("Derived W2 is not positive");

  const double pf = std::sqrt(std::max(0.0, Eprime * Eprime - mf * mf));
  double ctl = 1.0;
  if (pi > 0.0 && pf > 0.0) {
    ctl = (-mi * mi - mf * mf + 2.0 * E * Eprime - Q2) / (2.0 * pi * pf);
  }
  if (ctl < -1.0 - 5e-4 || ctl > 1.0 + 5e-4) {
    throw OutOfPhaseSpaceError("Derived costheta_l is outside [-1,1]");
  }
  ctl = clamp(ctl, -1.0, 1.0);

  KinePoint kp;
  kp.E = E;
  kp.Eprime = Eprime;
  kp.Q2 = Q2;
  kp.W = std::sqrt(W2);
  kp.omega = omega;
  kp.q3 = std::sqrt(std::max(0.0, omega * omega + Q2));
  kp.xB = (omega != 0.0) ? Q2 / (2.0 * M * omega) : std::numeric_limits<double>::quiet_NaN();
  kp.y = omega / E;
  kp.costheta_l = ctl;
  kp.theta_l_deg = std::acos(ctl) / kDegToRad;
  kp.Tl = Eprime - mf;

  if (has(values, Var::Eprime)) check_close("Eprime", get(values, Var::Eprime), kp.Eprime);
  if (has(values, Var::Tl)) check_close("Tl", get(values, Var::Tl), kp.Tl);
  if (has(values, Var::Omega)) check_close("omega", get(values, Var::Omega), kp.omega);
  if (has(values, Var::Y)) check_close("y", get(values, Var::Y), kp.y);
  if (has(values, Var::Q2)) check_close("Q2", get(values, Var::Q2), kp.Q2);
  if (has(values, Var::W)) check_close("W", get(values, Var::W), kp.W);
  if (has(values, Var::Q3)) check_close("q3", get(values, Var::Q3), kp.q3);
  if (has(values, Var::XB)) check_close("xB", get(values, Var::XB), kp.xB);
  if (has(values, Var::CosThetaL)) check_close("costheta_l", get(values, Var::CosThetaL), kp.costheta_l);
  if (has(values, Var::ThetaLDeg)) check_close("theta_l_deg", get(values, Var::ThetaLDeg), kp.theta_l_deg);

  return kp;
}

double value_of(const KinePoint &kp, Var v)
{
  switch (v) {
    case Var::E: return kp.E;
    case Var::Eprime: return kp.Eprime;
    case Var::Q2: return kp.Q2;
    case Var::W: return kp.W;
    case Var::Omega: return kp.omega;
    case Var::Q3: return kp.q3;
    case Var::XB: return kp.xB;
    case Var::Y: return kp.y;
    case Var::CosThetaL: return kp.costheta_l;
    case Var::ThetaLDeg: return kp.theta_l_deg;
    case Var::Tl: return kp.Tl;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

std::map<Var, double> output_map_from_kine(const KinePoint &kp)
{
  std::map<Var, double> out;
  for (Var v : {Var::E, Var::Eprime, Var::Q2, Var::W, Var::Omega,
                Var::Q3, Var::XB, Var::Y, Var::CosThetaL,
                Var::ThetaLDeg, Var::Tl}) {
    out[v] = value_of(kp, v);
  }
  return out;
}

void set_probe_p4(genie::InitialState *init, double E)
{
  const double m = particle_mass(init->ProbePdg());
  const double pz = std::sqrt(std::max(0.0, E * E - m * m));
  TLorentzVector p4(0.0, 0.0, pz, E);
  init->SetProbeP4(p4);
}

void apply_kinematics(genie::Interaction *interaction, const KinePoint &kp)
{
  set_probe_p4(interaction->InitStatePtr(), kp.E);

  genie::Kinematics *k = interaction->KinePtr();
  k->ClearRunningValues();
  k->SetW(kp.W);
  k->SetQ2(kp.Q2);
  k->Setq2(-kp.Q2);
  k->Setx(kp.xB);
  k->Sety(kp.y);
  k->SetKV(genie::kKVTl, kp.Tl);
  k->SetKV(genie::kKVctl, kp.costheta_l);
  k->SetKV(genie::kKVv, kp.omega);
  k->SetKV(genie::kKVQ0, kp.omega);
  k->SetKV(genie::kKVQ3, kp.q3);

  const double mf = final_lepton_mass(*interaction);
  const double pf = std::sqrt(std::max(0.0, kp.Eprime * kp.Eprime - mf * mf));
  const double sintheta = std::sqrt(std::max(0.0, 1.0 - kp.costheta_l * kp.costheta_l));
  TLorentzVector p4l(pf * sintheta, 0.0, pf * kp.costheta_l, kp.Eprime);
  k->SetFSLeptonP4(p4l);

  const double mi = particle_mass(interaction->InitState().ProbePdg());
  const double pi = std::sqrt(std::max(0.0, kp.E * kp.E - mi * mi));
  TLorentzVector q4(-p4l.Px(), -p4l.Py(), pi - p4l.Pz(), kp.omega);
  TLorentzVector had4 = interaction->InitState().Tgt().HitNucP4() + q4;
  k->SetHadSystP4(had4);
}

genie::InitialState make_initial_state(const Options &opt, double E)
{
  genie::InitialState init(opt.target, opt.probe);
  set_probe_p4(&init, E);
  return init;
}

std::pair<std::string, std::string> split_alg_id(const std::string &spec,
                                                 const std::string &default_config)
{
  const auto slash = spec.find('/');
  if (slash == std::string::npos) return {spec, default_config};
  return {spec.substr(0, slash), spec.substr(slash + 1)};
}

const genie::XSecAlgorithmI *get_xsec_alg(const std::string &spec)
{
  const auto [name, config] = split_alg_id(spec, "Default");
  const genie::Algorithm *alg = genie::AlgFactory::Instance()->GetAlgorithm(name, config);
  const auto *model = dynamic_cast<const genie::XSecAlgorithmI *>(alg);
  if (!model) throw std::runtime_error("Not an XSecAlgorithmI: " + spec);
  return model;
}

std::pair<std::string, std::string> process_to_intlist(const std::string &process)
{
  std::string p = process;
  std::replace(p.begin(), p.end(), '_', '-');
  const auto dash = p.find('-');
  if (dash == std::string::npos) {
    throw std::runtime_error("Process must look like DIS-EM, RES-CC, QEL-NC, MEC-EM: " + process);
  }
  std::string sc = lower(p.substr(0, dash));
  std::string current = p.substr(dash + 1);
  std::transform(current.begin(), current.end(), current.begin(),
                 [](unsigned char c) { return std::toupper(c); });

  std::string alg;
  if (sc == "qel" || sc == "qe" || sc == "qes") alg = "genie::QELInteractionListGenerator";
  else if (sc == "res") alg = "genie::RESInteractionListGenerator";
  else if (sc == "dis") alg = "genie::DISInteractionListGenerator";
  else if (sc == "mec" || sc == "2p2h") alg = "genie::MECInteractionListGenerator";
  else throw std::runtime_error("Unsupported first-release process family: " + process);

  std::string config;
  if (current == "EM") config = "EM-Default";
  else if (current == "CC") config = "CC-Default";
  else if (current == "NC") config = "NC-Default";
  else throw std::runtime_error("Unsupported process current: " + process);
  return {alg, config};
}

std::unique_ptr<genie::InteractionList> make_algorithm_interactions(
    const Options &opt, const genie::InitialState &init)
{
  std::string alg_name = opt.intlist_alg;
  std::string alg_config = opt.intlist_config;
  if (alg_name.empty()) {
    const auto mapped = process_to_intlist(opt.process);
    alg_name = mapped.first;
    alg_config = mapped.second;
  }

  const genie::Algorithm *alg =
      genie::AlgFactory::Instance()->GetAlgorithm(alg_name, alg_config);
  const auto *generator =
      dynamic_cast<const genie::InteractionListGeneratorI *>(alg);
  if (!generator) {
    throw std::runtime_error("Not an InteractionListGeneratorI: " + alg_name + "/" + alg_config);
  }
  std::unique_ptr<genie::InteractionList> list(generator->CreateInteractionList(init));
  if (!list || list->empty()) {
    throw std::runtime_error("Interaction list is empty for " + alg_name + "/" + alg_config);
  }
  return list;
}

std::string hit_state_key(const genie::Interaction &i)
{
  std::ostringstream ss;
  ss << i.ProcInfo().ScatteringTypeId() << ":"
     << i.ProcInfo().InteractionTypeId() << ":"
     << i.InitState().Tgt().HitNucPdg() << ":"
     << i.InitState().Tgt().HitQrkPdg() << ":"
     << i.InitState().Tgt().HitSeaQrk() << ":"
     << i.ExclTag().FinalQuarkPdg();
  return ss.str();
}

std::string component_name(const genie::Interaction &i,
                           const genie::XSecAlgorithmI &model)
{
  std::ostringstream ss;
  ss << genie::ScatteringType::AsString(i.ProcInfo().ScatteringTypeId())
     << ":" << genie::InteractionType::AsString(i.ProcInfo().InteractionTypeId())
     << ":hit=" << i.InitState().Tgt().HitNucPdg()
     << ":model=" << model.Id().Key();
  if (i.ExclTag().KnownResonance()) ss << ":res=" << i.ExclTag().Resonance();
  if (i.InitState().Tgt().HitQrkIsSet()) {
    ss << ":q=" << i.InitState().Tgt().HitQrkPdg()
       << (i.InitState().Tgt().HitSeaQrk() ? "sea" : "val");
  }
  return ss.str();
}

NativePS auto_native_ps(const genie::Interaction &interaction,
                        const genie::XSecAlgorithmI &model, int dim)
{
  const genie::ProcessInfo &proc = interaction.ProcInfo();
  if (dim == 1) {
    if (proc.IsQuasiElastic() || proc.IsMEC() || proc.IsElectronScattering()) {
      return NativePS::Q2;
    }
    return NativePS::W;
  }
  if (proc.IsMEC()) return NativePS::TlCtl;
  if (proc.IsQuasiElastic()) {
    const std::string name = model.Id().Name();
    const std::string key = model.Id().Key();
    if (name.find("SuSAv2") != std::string::npos ||
        name.find("HadronTensor") != std::string::npos ||
        key.find("SuSAv2") != std::string::npos ||
        key.find("HadronTensor") != std::string::npos) {
      return NativePS::TlCtl;
    }
    return NativePS::Q2;
  }
  return NativePS::WQ2;
}

genie::KinePhaseSpace_t phase_space(NativePS native)
{
  switch (native) {
    case NativePS::WQ2: return genie::kPSWQ2fE;
    case NativePS::TlCtl: return genie::kPSTlctl;
    case NativePS::XY: return genie::kPSxyfE;
    case NativePS::Q2: return genie::kPSQ2fE;
    case NativePS::W: return genie::kPSWfE;
    case NativePS::Auto: break;
  }
  throw std::runtime_error("Internal error: unresolved native phase space");
}

bool is_rosenbluth_qel(const genie::Interaction &interaction,
                       const genie::XSecAlgorithmI &model)
{
  if (!interaction.ProcInfo().IsQuasiElastic()) return false;
  const std::string name = model.Id().Name();
  const std::string key = model.Id().Key();
  return name.find("RosenbluthPXSec") != std::string::npos ||
         key.find("RosenbluthPXSec") != std::string::npos;
}

bool is_empirical_mec(const genie::Interaction &interaction,
                      const genie::XSecAlgorithmI &model)
{
  if (!interaction.ProcInfo().IsMEC()) return false;
  const std::string name = model.Id().Name();
  const std::string key = model.Id().Key();
  return name.find("EmpiricalMECPXSec2015") != std::string::npos ||
         key.find("EmpiricalMECPXSec2015") != std::string::npos;
}

bool qel_bin_fold_supports_diff(const std::vector<Var> &diff_vars)
{
  if (diff_vars.size() != 2) return false;
  bool has_energy = false;
  bool has_costh = false;
  for (Var v : diff_vars) {
    if (v == Var::Eprime || v == Var::Omega) has_energy = true;
    if (v == Var::CosThetaL) has_costh = true;
  }
  return has_energy && has_costh;
}

const genie::NuclearModelI *configured_nuclear_model(
    const genie::Interaction &interaction)
{
  genie::Registry *global =
      genie::AlgConfigPool::Instance()->GlobalParameterList();
  if (!global) return nullptr;

  const int tgt_pdg = interaction.InitState().Tgt().Pdg();
  std::vector<std::string> keys;
  if (tgt_pdg != 0) keys.push_back("NuclearModel@Pdg=" + std::to_string(tgt_pdg));
  keys.push_back("NuclearModel");

  for (const std::string &key : keys) {
    if (!global->Exists(key)) continue;
    const RgAlg alg = global->GetAlg(key);
    const genie::Algorithm *a =
        genie::AlgFactory::Instance()->GetAlgorithm(alg.name, alg.config);
    const auto *model = dynamic_cast<const genie::NuclearModelI *>(a);
    if (model) return model;
  }
  return nullptr;
}

double fallback_fermi_momentum(const genie::Target &target)
{
  const int A = target.A();
  if (A <= 1) return 0.0;
  if (A <= 3) return 0.115;
  if (A <= 7) return 0.190;
  if (A <= 16) return 0.221;
  if (A <= 25) return 0.230;
  if (A <= 38) return 0.236;
  if (A <= 60) return 0.241;
  return 0.245;
}

double qel_fold_fermi_momentum(const genie::Interaction &interaction,
                               const Options &opt)
{
  if (opt.qel_fold_kf >= 0.0) return opt.qel_fold_kf;
  const genie::Target &target = interaction.InitState().Tgt();
  const int hit = target.HitNucPdg();
  const genie::NuclearModelI *nucl_model = configured_nuclear_model(interaction);
  if (nucl_model && target.IsNucleus() && target.HitNucIsSet()) {
    const double kf = nucl_model->FermiMomentum(target, hit);
    if (std::isfinite(kf) && kf > 0.0) return kf;
  }
  return fallback_fermi_momentum(target);
}

double qel_fold_removal_energy(const genie::Interaction &interaction,
                               const Options &opt)
{
  if (opt.qel_fold_removal_energy >= 0.0) {
    return opt.qel_fold_removal_energy;
  }
  const genie::Target &target = interaction.InitState().Tgt();
  if (!target.IsNucleus()) return 0.0;

  const std::string key = "RFG-NucRemovalE@Pdg=" + std::to_string(target.Pdg());
  genie::Registry *r =
      genie::AlgConfigPool::Instance()->CommonList("Param", "FermiGas");
  if (r && r->Exists(key)) {
    const double eb = r->GetDouble(key);
    if (std::isfinite(eb) && eb >= 0.0) return eb;
  }

  if (target.A() <= 4) return 0.015;
  if (target.A() <= 16) return 0.025;
  if (target.A() <= 40) return 0.030;
  if (target.A() <= 60) return 0.036;
  return 0.044;
}

struct WeightedValue {
  double value = 0.0;
  double weight = 0.0;
};

struct QELFoldConfig {
  const genie::NuclearModelI *nucl_model = nullptr;
  genie::QELEvGen_BindingMode_t binding_mode = genie::kUseNuclearModel;
  double min_angle_em = 0.0;
  bool fermi_keep_on_shell = false;
  bool fermi_momentum_dependent_ermv = false;
  double nucl_r0 = 1.4;
};

struct QELFoldDebugCounters {
  long states = 0;
  long cos0_cells = 0;
  long q2_cells = 0;
  long kin_ok = 0;
  long q2_ok = 0;
  long jac_ok = 0;
  long factor_ok = 0;
  long xsec_ok = 0;
  long pauli_ok = 0;
  long pauli_blocked = 0;
  long bin_hits = 0;
  double cos0_max_min = std::numeric_limits<double>::infinity();
  double cos0_max_max = -std::numeric_limits<double>::infinity();
  double hit_e_min = std::numeric_limits<double>::infinity();
  double hit_e_max = -std::numeric_limits<double>::infinity();
  double sqrt_s_min = std::numeric_limits<double>::infinity();
  double sqrt_s_max = -std::numeric_limits<double>::infinity();
};

TVector3 qel_com_beta_lab(const genie::InitialState &initial_state)
{
  std::unique_ptr<TLorentzVector> probe(initial_state.GetProbeP4(genie::kRfLab));
  const TLorentzVector &hit = initial_state.Tgt().HitNucP4();
  TLorentzVector total = *probe + hit;
  return total.BoostVector();
}

bool qel_set_evgen_kinematics(genie::Interaction *interaction,
                              double cos_theta0, double phi0,
                              double min_angle_em)
{
  cos_theta0 = clamp(cos_theta0, -1.0, 1.0);

  TParticlePDG *lep = interaction->FSPrimLepton();
  if (!lep) return false;
  const double lep_mass = lep->Mass();

  TParticlePDG *recoil =
      TDatabasePDG::Instance()->GetParticle(interaction->RecoilNucleonPdg());
  if (!recoil) return false;
  const double mNf = recoil->Mass();

  const double sqrt_s = interaction->InitState().CMEnergy();
  if (sqrt_s <= lep_mass + mNf) return false;

  const double out_lep_E =
      (sqrt_s * sqrt_s - mNf * mNf + lep_mass * lep_mass) / (2.0 * sqrt_s);
  const double p2 = out_lep_E * out_lep_E - lep_mass * lep_mass;
  if (p2 < 0.0) return false;
  const double out_p = std::sqrt(p2);

  const TVector3 beta = qel_com_beta_lab(interaction->InitState());
  TVector3 lepton3(0.0, 0.0, out_p);
  lepton3.SetTheta(std::acos(cos_theta0));
  lepton3.SetPhi(phi0);

  TVector3 zvec(0.0, 0.0, 1.0);
  TVector3 rot = (zvec.Cross(beta)).Unit();
  double angle = beta.Angle(zvec);
  if (beta.Perp() == 0.0 && beta.Z() < 0.0) {
    rot = TVector3(0.0, 1.0, 0.0);
    angle = kPi;
  }
  if (rot.Mag() >= genie::controls::kASmallNum) lepton3.Rotate(angle, rot);

  TLorentzVector lepton(lepton3, out_lep_E);
  TLorentzVector out_nucleon(-lepton.Px(), -lepton.Py(), -lepton.Pz(),
                             std::sqrt(out_p * out_p + mNf * mNf));
  lepton.Boost(beta);
  out_nucleon.Boost(beta);

  if (interaction->ProcInfo().IsEM() &&
      180.0 * lepton.Theta() / kPi < min_angle_em) {
    return false;
  }

  std::unique_ptr<TLorentzVector> probe(
      interaction->InitState().GetProbeP4(genie::kRfLab));
  TLorentzVector q4 = *probe - lepton;
  const double Q2 = -q4.Mag2();
  if (!std::isfinite(Q2) || Q2 < 0.0) return false;

  genie::Kinematics *k = interaction->KinePtr();
  k->SetFSLeptonP4(lepton);
  k->SetHadSystP4(out_nucleon);
  k->SetQ2(Q2);
  k->Setq2(-Q2);
  return true;
}

double wrap_angle(double x)
{
  while (x <= -kPi) x += 2.0 * kPi;
  while (x > kPi) x -= 2.0 * kPi;
  return x;
}

double scattering_phi_hit_rest(const genie::Interaction &interaction)
{
  std::unique_ptr<TLorentzVector> probe(
      interaction.InitState().GetProbeP4(genie::kRfLab));
  TLorentzVector k = *probe;
  TLorentzVector l = interaction.Kine().FSLeptonP4();
  const TLorentzVector &hit = interaction.InitState().Tgt().HitNucP4();
  const TVector3 boost_to_hit_rest = -hit.BoostVector();
  k.Boost(boost_to_hit_rest);
  l.Boost(boost_to_hit_rest);

  TVector3 z = k.Vect();
  if (z.Mag() <= 0.0) return 0.0;
  z = z.Unit();
  TVector3 ref(1.0, 0.0, 0.0);
  if (std::fabs(ref.Dot(z)) > 0.95) ref = TVector3(0.0, 1.0, 0.0);
  TVector3 x = ref - z * ref.Dot(z);
  if (x.Mag() <= 0.0) return 0.0;
  x = x.Unit();
  TVector3 y = z.Cross(x).Unit();
  const TVector3 lp = l.Vect();
  return std::atan2(lp.Dot(y), lp.Dot(x));
}

bool qel_evgen_map(const genie::Interaction &bound_interaction,
                   double cos_theta0, double phi0, double min_angle_em,
                   double &Q2, double &phi_hit_rest)
{
  genie::Interaction tmp(bound_interaction);
  if (!qel_set_evgen_kinematics(&tmp, cos_theta0, phi0, min_angle_em)) {
    return false;
  }
  Q2 = tmp.Kine().Q2();
  phi_hit_rest = scattering_phi_hit_rest(tmp);
  return std::isfinite(Q2) && std::isfinite(phi_hit_rest);
}

double qel_evgen_jacobian_q2phi(const genie::Interaction &bound_interaction,
                                double cos_theta0, double phi0,
                                double min_angle_em,
                                double cos_theta0_max)
{
  const double h_u = 1e-4 * std::max(1.0, std::fabs(cos_theta0));
  const double h_phi = 1e-4;

  auto eval = [&](double u, double p, double &Q2, double &phi_hr) {
    return qel_evgen_map(bound_interaction, clamp(u, -1.0, cos_theta0_max),
                         p, min_angle_em, Q2, phi_hr);
  };

  double Q2_up = 0.0, ph_up = 0.0, Q2_um = 0.0, ph_um = 0.0;
  double u_plus = std::min(cos_theta0_max, cos_theta0 + h_u);
  double u_minus = std::max(-1.0, cos_theta0 - h_u);
  if (u_plus == u_minus ||
      !eval(u_plus, phi0, Q2_up, ph_up) ||
      !eval(u_minus, phi0, Q2_um, ph_um)) {
    return 0.0;
  }
  const double du = u_plus - u_minus;
  const double dQ2_du = (Q2_up - Q2_um) / du;
  const double dphihr_du = wrap_angle(ph_up - ph_um) / du;

  double Q2_pp = 0.0, ph_pp = 0.0, Q2_pm = 0.0, ph_pm = 0.0;
  if (!eval(cos_theta0, phi0 + h_phi, Q2_pp, ph_pp) ||
      !eval(cos_theta0, phi0 - h_phi, Q2_pm, ph_pm)) {
    return 0.0;
  }
  const double dQ2_dphi = (Q2_pp - Q2_pm) / (2.0 * h_phi);
  const double dphihr_dphi = wrap_angle(ph_pp - ph_pm) / (2.0 * h_phi);

  return std::fabs(dQ2_du * dphihr_dphi - dQ2_dphi * dphihr_du);
}

QELFoldConfig qel_fold_config(const genie::Interaction &interaction)
{
  QELFoldConfig cfg;
  cfg.nucl_model = configured_nuclear_model(interaction);
  if (!cfg.nucl_model) {
    throw std::runtime_error("--qel-bin-fold could not find configured NuclearModel");
  }

  std::string gen_config = interaction.ProcInfo().IsEM() ? "EM-Default" : "Default";
  const genie::Algorithm *qel_gen =
      genie::AlgFactory::Instance()->GetAlgorithm("genie::QELEventGenerator",
                                                  gen_config);
  const genie::Registry &qel_gen_config = qel_gen->GetConfig();
  std::string binding_mode = "UseNuclearModel";
  if (qel_gen_config.Exists("HitNucleonBindingMode")) {
    binding_mode = qel_gen_config.GetString("HitNucleonBindingMode");
  }
  cfg.binding_mode = genie::utils::StringToQELBindingMode(binding_mode);
  if (qel_gen_config.Exists("SF-MinAngleEMscattering")) {
    cfg.min_angle_em = qel_gen_config.GetDouble("SF-MinAngleEMscattering");
  }

  const genie::Algorithm *fermi =
      genie::AlgFactory::Instance()->GetAlgorithm("genie::FermiMover",
                                                  "Default");
  const genie::Registry &fermi_config = fermi->GetConfig();
  if (fermi_config.Exists("KeepHitNuclOnMassShell")) {
    cfg.fermi_keep_on_shell =
        fermi_config.GetBool("KeepHitNuclOnMassShell");
  }
  bool lfg_mom_dep = false;
  if (fermi_config.Exists("LFG-MomentumDependentErmv")) {
    lfg_mom_dep = fermi_config.GetBool("LFG-MomentumDependentErmv");
  }
  cfg.fermi_momentum_dependent_ermv = lfg_mom_dep;
  if (fermi_config.Exists("MomentumDependentErmv")) {
    cfg.fermi_momentum_dependent_ermv =
        fermi_config.GetBool("MomentumDependentErmv");
  }

  genie::Registry *nucl =
      genie::AlgConfigPool::Instance()->CommonList("Param", "NUCL");
  if (nucl && nucl->Exists("NUCL-R0")) {
    const double r0 = nucl->GetDouble("NUCL-R0");
    if (std::isfinite(r0) && r0 > 0.0) cfg.nucl_r0 = r0;
  }
  return cfg;
}

std::vector<WeightedValue> qel_radius_grid(const genie::Interaction &interaction,
                                           const QELFoldConfig &cfg,
                                           const Options &opt)
{
  const genie::Target &target = interaction.InitState().Tgt();
  if (!target.IsNucleus() ||
      cfg.nucl_model->ModelType(target) != genie::kNucmLocalFermiGas) {
    return {{0.0, 1.0}};
  }

  const double rmax =
      3.0 * cfg.nucl_r0 * std::pow(static_cast<double>(target.A()), 1.0 / 3.0);
  const double dr = rmax / static_cast<double>(opt.qel_fold_nr);
  std::vector<WeightedValue> grid;
  double norm = 0.0;
  for (int i = 0; i < opt.qel_fold_nr; ++i) {
    const double r = (static_cast<double>(i) + 0.5) * dr;
    const double w = r * r * genie::utils::nuclear::Density(r, target.A()) * dr;
    if (std::isfinite(w) && w > 0.0) {
      grid.push_back({r, w});
      norm += w;
    }
  }
  if (norm <= 0.0) return {{0.0, 1.0}};
  for (auto &x : grid) x.weight /= norm;
  return grid;
}

double qel_auto_pmax(const genie::Interaction &interaction,
                     const QELFoldConfig &cfg, double radius,
                     const Options &opt)
{
  if (opt.qel_fold_pmax > 0.0) return opt.qel_fold_pmax;
  if (opt.qel_fold_kf > 0.0) return opt.qel_fold_kf;

  const genie::Target &target = interaction.InitState().Tgt();
  if (!target.IsNucleus()) return 0.0;

  double max_p = 0.0;
  const double search_max = 4.0;
  const int nscan = 160;
  for (int i = 0; i < nscan; ++i) {
    const double p = (static_cast<double>(i) + 0.5) * search_max / nscan;
    const double prob = cfg.nucl_model->Prob(p, -1.0, target, radius);
    if (std::isfinite(prob) && prob > 0.0) max_p = p;
  }
  if (max_p > 0.0) return std::min(search_max, max_p + search_max / nscan);
  return qel_fold_fermi_momentum(interaction, opt);
}

std::vector<WeightedValue> qel_momentum_grid(const genie::Interaction &interaction,
                                             const QELFoldConfig &cfg,
                                             double radius,
                                             const Options &opt)
{
  const genie::Target &target = interaction.InitState().Tgt();
  if (!target.IsNucleus()) return {{0.0, 1.0}};

  const double pmax = qel_auto_pmax(interaction, cfg, radius, opt);
  if (pmax <= 0.0) return {{0.0, 1.0}};

  std::vector<WeightedValue> grid;
  double norm = 0.0;
  for (int i = 0; i < opt.qel_fold_np; ++i) {
    const double p = (static_cast<double>(i) + 0.5) * pmax /
                     static_cast<double>(opt.qel_fold_np);
    double w = 0.0;
    if (opt.qel_fold_kf > 0.0) {
      w = (p <= opt.qel_fold_kf) ? p * p : 0.0;
    } else {
      w = cfg.nucl_model->Prob(p, -1.0, target, radius);
    }
    if (std::isfinite(w) && w > 0.0) {
      grid.push_back({p, w});
      norm += w;
    }
  }
  if (norm <= 0.0) return {{0.0, 1.0}};
  for (auto &x : grid) x.weight /= norm;
  return grid;
}

double radical_inverse(std::uint64_t n, int base)
{
  double inv_base = 1.0 / static_cast<double>(base);
  double f = inv_base;
  double x = 0.0;
  while (n > 0) {
    x += f * static_cast<double>(n % static_cast<std::uint64_t>(base));
    n /= static_cast<std::uint64_t>(base);
    f *= inv_base;
  }
  return x;
}

double qel_lattice_u(std::uint64_t sample, int dim)
{
  static constexpr int bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  if (dim < 0 || dim >= static_cast<int>(sizeof(bases) / sizeof(bases[0]))) {
    throw std::runtime_error("internal error: lattice dimension out of range");
  }
  return radical_inverse(sample + 1, bases[dim]);
}

std::size_t select_weighted_index(const std::vector<WeightedValue> &grid,
                                  double u)
{
  if (grid.empty()) {
    throw std::runtime_error("internal error: empty weighted grid");
  }
  u = clamp(u, 0.0, std::nextafter(1.0, 0.0));
  double acc = 0.0;
  for (std::size_t i = 0; i < grid.size(); ++i) {
    acc += grid[i].weight;
    if (u < acc) return i;
  }
  return grid.size() - 1;
}

bool qel_bind_old_fermi_mover(genie::Interaction *interaction,
                              const QELFoldConfig &cfg,
                              const TVector3 &p3,
                              double removal_energy)
{
  genie::Target *target = interaction->InitStatePtr()->TgtPtr();
  if (!target || !target->HitNucIsSet()) return false;

  TParticlePDG *hit_pdg =
      TDatabasePDG::Instance()->GetParticle(target->HitNucPdg());
  if (!hit_pdg) return false;
  const double mN = hit_pdg->Mass();
  const double p2 = p3.Mag2();

  double EN = std::sqrt(std::max(0.0, mN * mN + p2));
  if (target->IsNucleus() && !cfg.fermi_keep_on_shell) {
    const genie::FermiMoverInteractionType_t fermi_type =
        cfg.nucl_model->GetFermiMoverInteractionType();
    if (fermi_type == genie::kFermiMoveEffectiveSF1p1h ||
        cfg.fermi_momentum_dependent_ermv) {
      EN = mN - removal_energy -
           p2 / (2.0 * (target->Mass() - mN));
    } else if (fermi_type == genie::kFermiMoveEffectiveSF2p2h_eject ||
               fermi_type == genie::kFermiMoveEffectiveSF2p2h_noeject) {
      const int other_pdg = genie::pdg::IsProton(target->HitNucPdg())
                                ? genie::kPdgNeutron
                                : genie::kPdgProton;
      TParticlePDG *other = genie::PDGLibrary::Instance()->Find(other_pdg);
      TParticlePDG *deuteron =
          genie::PDGLibrary::Instance()->Find(1000010020);
      if (!other || !deuteron) return false;
      EN = deuteron->Mass() - 2.0 * removal_energy -
           std::sqrt(std::max(0.0, p2 + other->Mass() * other->Mass()));
    } else {
      const bool is_proton = genie::pdg::IsProton(target->HitNucPdg());
      const int remnant_z = is_proton ? target->Z() - 1 : target->Z();
      const int remnant_a = target->A() - 1;
      TParticlePDG *remnant = genie::PDGLibrary::Instance()->Find(
          genie::pdg::IonPdgCode(remnant_a, remnant_z));
      if (!remnant) return false;
      EN = target->Mass() -
           std::sqrt(std::max(0.0, p2 + remnant->Mass() * remnant->Mass()));
    }
  }

  if (!std::isfinite(EN)) return false;
  TLorentzVector *hit_p4 = target->HitNucP4Ptr();
  hit_p4->SetPxPyPzE(p3.Px(), p3.Py(), p3.Pz(), EN);
  return true;
}

bool qel_lock_old_q2_kinematics(genie::Interaction *interaction,
                                double q2)
{
  TParticlePDG *recoil =
      TDatabasePDG::Instance()->GetParticle(interaction->RecoilNucleonPdg());
  if (!recoil) return false;

  const genie::InitialState &init = interaction->InitState();
  const double E_hit_rest = init.ProbeE(genie::kRfHitNucRest);
  const double M_hit = init.Tgt().HitNucP4().M();
  const double W = recoil->Mass();
  if (!std::isfinite(E_hit_rest) || !std::isfinite(M_hit) ||
      E_hit_rest <= 0.0 || M_hit <= 0.0 || W <= 0.0) {
    return false;
  }

  double x = 0.0;
  double y = 0.0;
  genie::utils::kinematics::WQ2toXY(E_hit_rest, M_hit, W, q2, x, y);
  if (!std::isfinite(x) || !std::isfinite(y)) return false;

  genie::Kinematics *k = interaction->KinePtr();
  k->SetQ2(q2, true);
  k->Setq2(-q2);
  k->SetW(W, true);
  k->Setx(x, true);
  k->Sety(y, true);
  k->ClearRunningValues();
  return true;
}

bool qel_set_old_primary_lepton(genie::Interaction *interaction,
                                double phi_hit_rest)
{
  const TLorentzVector &hit = interaction->InitState().Tgt().HitNucP4();
  const TVector3 beta = hit.BoostVector();
  std::unique_ptr<TLorentzVector> probe_lab(
      interaction->InitState().GetProbeP4(genie::kRfLab));
  if (!probe_lab) return false;

  TLorentzVector probe_hit_rest(*probe_lab);
  probe_hit_rest.Boost(-1.0 * beta);
  const double Ev = probe_hit_rest.E();
  const double Q2 = interaction->Kine().Q2(true);
  const double y = interaction->Kine().y(true);
  const double ml = final_lepton_mass(*interaction);
  const double ml2 = ml * ml;
  if (!std::isfinite(Ev) || !std::isfinite(Q2) || !std::isfinite(y) ||
      Ev <= 0.0) {
    return false;
  }

  const double El = (1.0 - y) * Ev;
  const double plp = El - 0.5 * (Q2 + ml2) / Ev;
  const double plt2 = El * El - plp * plp - ml2;
  if (!std::isfinite(El) || !std::isfinite(plp) || plt2 < -1e-10 ||
      El < ml) {
    return false;
  }
  const double plt = std::sqrt(std::max(0.0, plt2));
  TVector3 unit_probe = probe_hit_rest.Vect();
  if (unit_probe.Mag() <= 0.0) return false;
  unit_probe = unit_probe.Unit();

  TVector3 lepton3(plt * std::cos(phi_hit_rest),
                   plt * std::sin(phi_hit_rest), plp);
  lepton3.RotateUz(unit_probe);
  TLorentzVector lepton(lepton3, El);
  lepton.Boost(beta);
  if (!std::isfinite(lepton.E()) || lepton.E() <= ml) return false;

  genie::Kinematics *k = interaction->KinePtr();
  k->SetFSLeptonP4(lepton);

  TLorentzVector q4 = *probe_lab - lepton;
  k->SetHadSystP4(hit + q4);
  return true;
}

bool qel_set_old_generator_q2_kinematics(genie::Interaction *interaction,
                                         double q2,
                                         double phi_hit_rest)
{
  if (!qel_lock_old_q2_kinematics(interaction, q2)) return false;
  return qel_set_old_primary_lepton(interaction, phi_hit_rest);
}

bool qel_passes_old_pauli_blocker(const genie::Interaction &interaction,
                                  const QELFoldConfig &cfg,
                                  QELFoldDebugCounters *debug = nullptr)
{
  const genie::Target &target = interaction.InitState().Tgt();
  if (!target.IsNucleus()) {
    if (debug) ++debug->pauli_ok;
    return true;
  }

  const int recoil_pdg = interaction.RecoilNucleonPdg();
  if (!genie::pdg::IsProton(recoil_pdg) &&
      !genie::pdg::IsNeutron(recoil_pdg)) {
    if (debug) ++debug->pauli_ok;
    return true;
  }

  double kf = 0.0;
  if (cfg.nucl_model &&
      cfg.nucl_model->ModelType(target) == genie::kNucmLocalFermiGas) {
    const bool is_p = genie::pdg::IsProton(recoil_pdg);
    const double num_nucleons = is_p ? target.Z() : target.N();
    const double hbarc =
        genie::constants::kLightSpeed * genie::constants::kPlankConstant /
        genie::units::fermi;
    kf = std::pow(3.0 * genie::constants::kPi2 * num_nucleons *
                      genie::utils::nuclear::Density(target.HitNucPosition(),
                                                     target.A()),
                  1.0 / 3.0) *
         hbarc;
  } else if (cfg.nucl_model) {
    kf = cfg.nucl_model->FermiMomentum(target, recoil_pdg);
  }

  const double recoil_p = interaction.Kine().HadSystP4().P();
  const bool pass = !std::isfinite(kf) || kf <= 0.0 || recoil_p >= kf;
  if (debug) {
    if (pass) ++debug->pauli_ok;
    else ++debug->pauli_blocked;
  }
  return pass;
}

double qel_hit_nucleon_count(const genie::Target &target)
{
  if (!target.IsNucleus()) return 1.0;
  const int hit = target.HitNucPdg();
  if (genie::pdg::IsProton(hit)) return static_cast<double>(target.Z());
  if (genie::pdg::IsNeutron(hit)) return static_cast<double>(target.N());
  return 1.0;
}

double rosenbluth_qel_old_generator_q2_xsec(
    const genie::XSecAlgorithmI &model,
    genie::Interaction *interaction,
    const QELFoldConfig &cfg,
    double q2,
    double phi_hit_rest,
    QELFoldDebugCounters *debug = nullptr)
{
  interaction->KinePtr()->SetQ2(q2);
  const double xsec = model.XSec(interaction, genie::kPSQ2fE);
  if (!std::isfinite(xsec) || xsec <= 0.0) return 0.0;
  if (debug) ++debug->q2_ok;

  interaction->ResetBit(genie::kISkipProcessChk);
  interaction->ResetBit(genie::kISkipKinematicChk);
  interaction->ResetBit(genie::kIAssumeFreeNucleon);

  if (!qel_set_old_generator_q2_kinematics(interaction, q2,
                                           phi_hit_rest)) {
    return 0.0;
  }
  if (debug) ++debug->kin_ok;
  if (!qel_passes_old_pauli_blocker(*interaction, cfg, debug)) return 0.0;
  if (debug) ++debug->xsec_ok;
  return xsec;
}

double rosenbluth_qel_evgen_xsec(const genie::XSecAlgorithmI &model,
                                 const genie::Interaction &bound_interaction,
                                 double cos_theta0, double phi0,
                                 double cos_theta0_max,
                                 const QELFoldConfig &cfg,
                                 const Options &opt,
                                 QELFoldDebugCounters *debug = nullptr)
{
  genie::Interaction interaction(bound_interaction);
  if (!qel_set_evgen_kinematics(&interaction, cos_theta0, phi0,
                                cfg.min_angle_em)) {
    return 0.0;
  }
  if (debug) ++debug->kin_ok;
  if (!model.ValidProcess(&interaction)) return 0.0;

  const genie::Range1D_t q2lim = interaction.PhaseSpace().Q2Lim();
  const double q2 = interaction.Kine().Q2();
  if (q2 < q2lim.min || q2 > q2lim.max) return 0.0;

  const double dq2 = model.XSec(&interaction, genie::kPSQ2fE);
  if (!std::isfinite(dq2) || dq2 <= 0.0) return 0.0;
  if (debug) ++debug->q2_ok;

  double xsec = 0.0;
  if (opt.qel_fold_density == "qel-delta") {
    const double factor =
        genie::utils::EnergyDeltaFunctionSolutionQEL(interaction);
    if (!std::isfinite(factor) || factor <= 0.0) return 0.0;
    if (debug) ++debug->factor_ok;
    xsec = dq2 * factor;
  } else {
    const double jac = qel_evgen_jacobian_q2phi(bound_interaction, cos_theta0,
                                                phi0, cfg.min_angle_em,
                                                cos_theta0_max);
    if (!std::isfinite(jac) || jac <= 0.0) return 0.0;
    if (debug) ++debug->jac_ok;
    xsec = dq2 * jac / (2.0 * kPi);
  }

  if (!std::isfinite(xsec) || xsec <= 0.0) return 0.0;
  if (debug) ++debug->xsec_ok;
  return xsec;
}

bool qel_lab_lepton_in_bin(const genie::Interaction &interaction,
                           const KinePoint &center, const Options &opt)
{
  const TLorentzVector &lep = interaction.Kine().FSLeptonP4();
  const double p = lep.P();
  if (p <= 0.0) return false;
  const double Eprime = lep.E();
  const double costh = lep.Pz() / p;
  return std::fabs(Eprime - center.Eprime) <= 0.5 * opt.qel_bin_width_energy &&
         std::fabs(costh - center.costheta_l) <= 0.5 * opt.qel_bin_width_costh;
}

std::vector<double> scan_centers(const RangeSpec &range)
{
  std::vector<double> centers;
  centers.reserve(range.n);
  for (int i = 0; i < range.n; ++i) {
    const double t = (range.n == 1)
                         ? 0.0
                         : static_cast<double>(i) / static_cast<double>(range.n - 1);
    centers.push_back(range.min + t * (range.max - range.min));
  }
  return centers;
}

std::optional<int> scan_index_for_value(const RangeSpec &range, double value)
{
  if (range.n <= 0) return std::nullopt;
  if (range.n == 1) return 0;
  const double denom = range.max - range.min;
  if (denom == 0.0) return std::nullopt;
  const double t = (value - range.min) / denom;
  const int idx = static_cast<int>(
      std::llround(t * static_cast<double>(range.n - 1)));
  if (idx < 0 || idx >= range.n) return std::nullopt;
  const double center = range.min +
                        static_cast<double>(idx) /
                            static_cast<double>(range.n - 1) * denom;
  const double step = std::fabs(denom) / static_cast<double>(range.n - 1);
  if (std::fabs(value - center) > std::max(1e-10, 1e-6 * step)) {
    return std::nullopt;
  }
  return idx;
}

bool qel_generator_q2_scan_cache_supported(const Options &opt,
                                           const std::map<Var, double> &values)
{
  if (!opt.qel_fold_scan_cache) return false;
  if (opt.observable != Observable::D2) return false;
  if (!opt.integrate_over.empty() || !opt.points_path.empty()) return false;
  if (opt.qel_fold_method != "lattice" ||
      opt.qel_fold_density != "generator-q2") {
    return false;
  }
  if (!qel_bin_fold_supports_diff(opt.diff_vars)) return false;
  if (opt.scans.size() != 1) return false;
  const Var scan_var = opt.scans.front().var;
  if (scan_var != Var::Eprime && scan_var != Var::Omega) return false;
  return has(values, Var::E) && has(values, Var::CosThetaL) &&
         has(values, scan_var);
}

std::string qel_scan_fold_cache_key(const genie::XSecAlgorithmI &model,
                                    const genie::Interaction &base_interaction,
                                    const std::map<Var, double> &values,
                                    const Options &opt)
{
  const RangeSpec &range = opt.scans.front();
  std::ostringstream key;
  key << std::setprecision(17)
      << "generator-q2-scan|"
      << model.Id().Key() << "|"
      << component_name(base_interaction, model) << "|"
      << "target=" << base_interaction.InitState().Tgt().Pdg() << "|"
      << "hit=" << base_interaction.InitState().Tgt().HitNucPdg() << "|"
      << "probe=" << base_interaction.InitState().ProbePdg() << "|"
      << "E=" << get(values, Var::E) << "|"
      << "costh=" << get(values, Var::CosThetaL) << "|"
      << "scan=" << canonical_var_name(range.var) << ":"
      << range.min << ":" << range.max << ":" << range.n << "|"
      << "dE=" << opt.qel_bin_width_energy << "|"
      << "dcth=" << opt.qel_bin_width_costh << "|"
      << "samples=" << opt.qel_fold_samples << "|"
      << "nr=" << opt.qel_fold_nr << "|"
      << "np=" << opt.qel_fold_np << "|"
      << "ncp=" << opt.qel_fold_ncosth_p << "|"
      << "npp=" << opt.qel_fold_nphi_p << "|"
      << "nq2=" << opt.qel_fold_ncos0 << "|"
      << "nphi=" << opt.qel_fold_nphi0 << "|"
      << "pmax=" << opt.qel_fold_pmax << "|"
      << "kf=" << opt.qel_fold_kf << "|"
      << "removal=" << opt.qel_fold_removal_energy << "|"
      << "emq2=" << (opt.em_q2_min_set ? opt.em_q2_min : -1.0);
  return key.str();
}

bool qel_fill_scan_bins(std::vector<double> &folded,
                        const std::vector<double> &centers,
                        double observed_energy_var, double weight,
                        double half_width)
{
  bool filled = false;
  for (std::size_t i = 0; i < centers.size(); ++i) {
    if (std::fabs(observed_energy_var - centers[i]) <= half_width) {
      folded[i] += weight;
      filled = true;
    }
  }
  return filled;
}

std::optional<double> rosenbluth_qel_lattice_scan_cache_density_generator_q2(
    const genie::XSecAlgorithmI &model,
    const genie::Interaction &base_interaction,
    const std::map<Var, double> &values,
    const Options &opt)
{
  if (!qel_generator_q2_scan_cache_supported(opt, values)) {
    return std::nullopt;
  }

  const RangeSpec &range = opt.scans.front();
  const std::optional<int> point_idx =
      scan_index_for_value(range, get(values, range.var));
  if (!point_idx) return std::nullopt;

  static std::map<std::string, std::vector<double>> cache;
  const std::string key =
      qel_scan_fold_cache_key(model, base_interaction, values, opt);
  auto cached = cache.find(key);
  if (cached == cache.end()) {
    const double E = get(values, Var::E);
    const double costh_center = get(values, Var::CosThetaL);
    const double half_energy = 0.5 * opt.qel_bin_width_energy;
    const double half_costh = 0.5 * opt.qel_bin_width_costh;
    const std::vector<double> centers = scan_centers(range);

    genie::Interaction fold_base(base_interaction);
    set_probe_p4(fold_base.InitStatePtr(), E);

    const QELFoldConfig cfg = qel_fold_config(fold_base);
    const std::vector<WeightedValue> radii =
        qel_radius_grid(fold_base, cfg, opt);
    const double removal_energy = qel_fold_removal_energy(fold_base, opt);

    std::vector<std::vector<WeightedValue>> pgrids;
    pgrids.reserve(radii.size());
    for (const WeightedValue &rv : radii) {
      pgrids.push_back(qel_momentum_grid(fold_base, cfg, rv.value, opt));
    }

    std::vector<double> folded(centers.size(), 0.0);
    QELFoldDebugCounters debug;
    long scan_bin_fills = 0;
    for (int is = 0; is < opt.qel_fold_samples; ++is) {
      const std::uint64_t sample = static_cast<std::uint64_t>(is);
      const std::size_t ir = select_weighted_index(radii, qel_lattice_u(sample, 0));
      const WeightedValue &rv = radii[ir];
      const std::vector<WeightedValue> &pgrid = pgrids[ir];
      const std::size_t ip = select_weighted_index(pgrid, qel_lattice_u(sample, 1));
      const WeightedValue &pv = pgrid[ip];

      const double ctp = -1.0 + 2.0 * qel_lattice_u(sample, 2);
      const double stp = std::sqrt(std::max(0.0, 1.0 - ctp * ctp));
      const double phip = 2.0 * kPi * qel_lattice_u(sample, 3);
      TVector3 p3(pv.value * stp * std::cos(phip),
                  pv.value * stp * std::sin(phip),
                  pv.value * ctp);

      genie::Interaction bound(fold_base);
      bound.SetBit(genie::kISkipProcessChk);
      bound.SetBit(genie::kISkipKinematicChk);
      bound.SetBit(genie::kIAssumeFreeNucleon);
      bound.InitStatePtr()->TgtPtr()->SetHitNucPosition(rv.value);
      cfg.nucl_model->SetMomentum3(p3);
      cfg.nucl_model->SetRemovalEnergy(removal_energy);
      if (!qel_bind_old_fermi_mover(&bound, cfg, p3, removal_energy)) {
        continue;
      }
      ++debug.states;

      const double hit_e = bound.InitState().Tgt().HitNucP4().E();
      debug.hit_e_min = std::min(debug.hit_e_min, hit_e);
      debug.hit_e_max = std::max(debug.hit_e_max, hit_e);
      const double sqrt_s = bound.InitState().CMEnergy();
      debug.sqrt_s_min = std::min(debug.sqrt_s_min, sqrt_s);
      debug.sqrt_s_max = std::max(debug.sqrt_s_max, sqrt_s);

      const genie::Range1D_t q2lim = bound.PhaseSpace().Limits(genie::kKVQ2);
      if (!std::isfinite(q2lim.min) || !std::isfinite(q2lim.max) ||
          q2lim.max <= q2lim.min) {
        continue;
      }
      const double q2min = q2lim.min + genie::controls::kASmallNum;
      const double q2max = q2lim.max - genie::controls::kASmallNum;
      if (q2max <= q2min) continue;

      ++debug.q2_cells;
      const double q2 = q2min + (q2max - q2min) * qel_lattice_u(sample, 4);
      const double phi = 2.0 * kPi * qel_lattice_u(sample, 5);
      genie::Interaction generated(bound);
      const double xsec = rosenbluth_qel_old_generator_q2_xsec(
          model, &generated, cfg, q2, phi,
          opt.qel_fold_debug ? &debug : nullptr);
      if (xsec <= 0.0) continue;

      const TLorentzVector &lep = generated.Kine().FSLeptonP4();
      const double p = lep.P();
      if (p <= 0.0) continue;
      const double costh = lep.Pz() / p;
      if (std::fabs(costh - costh_center) > half_costh) continue;

      const double energy_var =
          (range.var == Var::Omega) ? (E - lep.E()) : lep.E();
      const double weight = xsec * (q2max - q2min);
      if (qel_fill_scan_bins(folded, centers, energy_var, weight, half_energy)) {
        ++debug.bin_hits;
        ++scan_bin_fills;
      }
    }

    const double norm = static_cast<double>(opt.qel_fold_samples) *
                        opt.qel_bin_width_energy * opt.qel_bin_width_costh;
    if (norm > 0.0) {
      for (double &x : folded) x /= norm;
    }

    if (opt.qel_fold_debug) {
      std::cerr << "qel-fold-debug method=lattice"
                << " density=generator-q2"
                << " cache=scan"
                << " component=" << component_name(base_interaction, model)
                << " samples=" << opt.qel_fold_samples
                << " scan_bins=" << centers.size()
                << " states=" << debug.states
                << " q2_cells=" << debug.q2_cells
                << " kin_ok=" << debug.kin_ok
                << " q2_ok=" << debug.q2_ok
                << " xsec_ok=" << debug.xsec_ok
                << " pauli_ok=" << debug.pauli_ok
                << " pauli_blocked=" << debug.pauli_blocked
                << " bin_hits=" << debug.bin_hits
                << " scan_bin_fills=" << scan_bin_fills
                << " hitE_min=" << debug.hit_e_min
                << " hitE_max=" << debug.hit_e_max
                << " sqrtS_min=" << debug.sqrt_s_min
                << " sqrtS_max=" << debug.sqrt_s_max << "\n";
    }

    cached = cache.emplace(key, std::move(folded)).first;
  }

  return cached->second.at(static_cast<std::size_t>(*point_idx));
}

TVector3 unit_perpendicular(const TVector3 &axis)
{
  TVector3 ref(1.0, 0.0, 0.0);
  if (std::fabs(axis.Dot(ref)) > 0.9) ref = TVector3(0.0, 1.0, 0.0);
  TVector3 out = ref - axis * ref.Dot(axis);
  if (out.Mag() <= 0.0) return TVector3(1.0, 0.0, 0.0);
  return out.Unit();
}

double rosenbluth_qel_exact_theta_density(
    const genie::XSecAlgorithmI &model,
    const genie::Interaction &base_interaction,
    const std::map<Var, double> &values,
    const Options &opt)
{
  if (!qel_bin_fold_supports_diff(opt.diff_vars)) {
    throw std::runtime_error(
        "--qel-bin-fold supports d2 in Eprime,costheta_l or omega,costheta_l only");
  }

  genie::Interaction center_interaction(base_interaction);
  const KinePoint center = solve_kinematics(values, center_interaction);
  genie::Interaction fold_base(base_interaction);
  set_probe_p4(fold_base.InitStatePtr(), center.E);

  const double mi = particle_mass(fold_base.InitState().ProbePdg());
  const double mf = final_lepton_mass(fold_base);
  const double pin = std::sqrt(std::max(0.0, center.E * center.E - mi * mi));
  const double pout =
      std::sqrt(std::max(0.0, center.Eprime * center.Eprime - mf * mf));
  if (pin <= 0.0 || pout <= 0.0) return 0.0;

  const double sinth =
      std::sqrt(std::max(0.0, 1.0 - center.costheta_l * center.costheta_l));
  const TVector3 kin(0.0, 0.0, pin);
  const TVector3 klep(pout * sinth, 0.0, pout * center.costheta_l);
  const TVector3 qvec = kin - klep;
  const double qmag = qvec.Mag();
  if (qmag <= 0.0) return 0.0;
  const TVector3 qhat = qvec.Unit();
  const TVector3 xhat = unit_perpendicular(qhat);
  const TVector3 yhat = qhat.Cross(xhat).Unit();

  TParticlePDG *recoil =
      TDatabasePDG::Instance()->GetParticle(base_interaction.RecoilNucleonPdg());
  if (!recoil) return 0.0;
  const double m_recoil = recoil->Mass();
  const double hit_count =
      qel_hit_nucleon_count(fold_base.InitState().Tgt());

  const QELFoldConfig cfg = qel_fold_config(fold_base);
  const std::vector<WeightedValue> radii =
      qel_radius_grid(fold_base, cfg, opt);
  const double removal_energy = qel_fold_removal_energy(fold_base, opt);

  double folded = 0.0;
  QELFoldDebugCounters debug;
  long delta_hits = 0;
  for (const WeightedValue &rv : radii) {
    std::vector<WeightedValue> pgrid =
        qel_momentum_grid(fold_base, cfg, rv.value, opt);
    for (const WeightedValue &pv : pgrid) {
      const double p = pv.value;
      if (p <= 0.0) continue;

      genie::Interaction energy_probe(fold_base);
      energy_probe.SetBit(genie::kISkipProcessChk);
      energy_probe.SetBit(genie::kISkipKinematicChk);
      energy_probe.SetBit(genie::kIAssumeFreeNucleon);
      energy_probe.InitStatePtr()->TgtPtr()->SetHitNucPosition(rv.value);
      const TVector3 p_for_energy = p * qhat;
      cfg.nucl_model->SetMomentum3(p_for_energy);
      cfg.nucl_model->SetRemovalEnergy(removal_energy);
      if (!qel_bind_old_fermi_mover(&energy_probe, cfg, p_for_energy,
                                    removal_energy)) {
        continue;
      }
      const double hit_e = energy_probe.InitState().Tgt().HitNucP4().E();
      if (!std::isfinite(hit_e)) continue;
      const double ef = hit_e + center.omega;
      if (ef <= m_recoil) continue;

      const double cos_pq =
          (ef * ef - m_recoil * m_recoil - p * p - qmag * qmag) /
          (2.0 * p * qmag);
      if (!std::isfinite(cos_pq) || cos_pq < -1.0 || cos_pq > 1.0) {
        continue;
      }
      const double sin_pq =
          std::sqrt(std::max(0.0, 1.0 - cos_pq * cos_pq));
      const double delta_factor = ef / (2.0 * p * qmag);
      if (!std::isfinite(delta_factor) || delta_factor <= 0.0) continue;

      double phi_sum = 0.0;
      int phi_ok = 0;
      for (int ia = 0; ia < opt.qel_fold_nphi_p; ++ia) {
        const double alpha =
            2.0 * kPi * (static_cast<double>(ia) + 0.5) /
            static_cast<double>(opt.qel_fold_nphi_p);
        const TVector3 transverse =
            std::cos(alpha) * xhat + std::sin(alpha) * yhat;
        const TVector3 p3 = p * (cos_pq * qhat + sin_pq * transverse);

        genie::Interaction fixed(fold_base);
        fixed.SetBit(genie::kISkipProcessChk);
        fixed.SetBit(genie::kISkipKinematicChk);
        fixed.SetBit(genie::kIAssumeFreeNucleon);
        fixed.InitStatePtr()->TgtPtr()->SetHitNucPosition(rv.value);
        cfg.nucl_model->SetMomentum3(p3);
        cfg.nucl_model->SetRemovalEnergy(removal_energy);
        if (!qel_bind_old_fermi_mover(&fixed, cfg, p3, removal_energy)) {
          continue;
        }
        ++debug.states;
        const double bound_hit_e = fixed.InitState().Tgt().HitNucP4().E();
        debug.hit_e_min = std::min(debug.hit_e_min, bound_hit_e);
        debug.hit_e_max = std::max(debug.hit_e_max, bound_hit_e);
        const double sqrt_s = fixed.InitState().CMEnergy();
        debug.sqrt_s_min = std::min(debug.sqrt_s_min, sqrt_s);
        debug.sqrt_s_max = std::max(debug.sqrt_s_max, sqrt_s);

        const genie::Range1D_t q2lim = fixed.PhaseSpace().Limits(genie::kKVQ2);
        if (center.Q2 < q2lim.min || center.Q2 > q2lim.max) continue;

        fixed.KinePtr()->SetQ2(center.Q2);
        const double dq2 = model.XSec(&fixed, genie::kPSQ2fE);
        if (!std::isfinite(dq2) || dq2 <= 0.0) continue;
        ++debug.q2_ok;

        TLorentzVector lepton_lab(klep.Px(), klep.Py(), klep.Pz(),
                                  center.Eprime);
        fixed.KinePtr()->SetFSLeptonP4(lepton_lab);
        const TVector3 pf3 = p3 + qvec;
        fixed.KinePtr()->SetHadSystP4(
            TLorentzVector(pf3.Px(), pf3.Py(), pf3.Pz(), ef));
        if (!qel_passes_old_pauli_blocker(fixed, cfg,
                                          opt.qel_fold_debug ? &debug
                                                             : nullptr)) {
          continue;
        }
        ++debug.xsec_ok;
        ++phi_ok;
        TLorentzVector lepton_hit_rest(lepton_lab);
        lepton_hit_rest.Boost(-1.0 * fixed.InitState().Tgt().HitNucP4().BoostVector());
        const double dcosth_density = dq2 * 2.0 * lepton_hit_rest.E() *
                                      lepton_hit_rest.E();
        if (std::isfinite(dcosth_density) && dcosth_density > 0.0) {
          phi_sum += dcosth_density;
        }
      }

      if (phi_ok > 0) {
        ++delta_hits;
        const double phi_avg =
            phi_sum / static_cast<double>(opt.qel_fold_nphi_p);
        folded += hit_count * rv.weight * pv.weight * delta_factor * phi_avg;
      }
    }
  }

  if (opt.qel_fold_debug) {
    std::cerr << "qel-fold-debug method=grid"
              << " density=exact-theta"
              << " component=" << component_name(base_interaction, model)
              << " states=" << debug.states
              << " delta_hits=" << delta_hits
              << " q2_ok=" << debug.q2_ok
              << " xsec_ok=" << debug.xsec_ok
              << " pauli_ok=" << debug.pauli_ok
              << " pauli_blocked=" << debug.pauli_blocked
              << " hitE_min=" << debug.hit_e_min
              << " hitE_max=" << debug.hit_e_max
              << " sqrtS_min=" << debug.sqrt_s_min
              << " sqrtS_max=" << debug.sqrt_s_max
              << " folded=" << folded << "\n";
  }
  return folded;
}

double rosenbluth_qel_grid_fold_density_generator_q2(
    const genie::XSecAlgorithmI &model,
    const genie::Interaction &base_interaction,
    const std::map<Var, double> &values,
    const Options &opt)
{
  if (!qel_bin_fold_supports_diff(opt.diff_vars)) {
    throw std::runtime_error(
        "--qel-bin-fold supports d2 in Eprime,costheta_l or omega,costheta_l only");
  }

  genie::Interaction center_interaction(base_interaction);
  const KinePoint center = solve_kinematics(values, center_interaction);
  genie::Interaction fold_base(base_interaction);
  set_probe_p4(fold_base.InitStatePtr(), center.E);

  const QELFoldConfig cfg = qel_fold_config(fold_base);
  const std::vector<WeightedValue> radii =
      qel_radius_grid(fold_base, cfg, opt);
  const double removal_energy = qel_fold_removal_energy(fold_base, opt);

  double folded = 0.0;
  QELFoldDebugCounters debug;
  for (const WeightedValue &rv : radii) {
    std::vector<WeightedValue> pgrid =
        qel_momentum_grid(fold_base, cfg, rv.value, opt);
    for (const WeightedValue &pv : pgrid) {
      for (int icp = 0; icp < opt.qel_fold_ncosth_p; ++icp) {
        const double ctp = -1.0 + 2.0 * (static_cast<double>(icp) + 0.5) /
                                      static_cast<double>(opt.qel_fold_ncosth_p);
        const double stp = std::sqrt(std::max(0.0, 1.0 - ctp * ctp));
        for (int ipp = 0; ipp < opt.qel_fold_nphi_p; ++ipp) {
          const double phip = 2.0 * kPi * (static_cast<double>(ipp) + 0.5) /
                              static_cast<double>(opt.qel_fold_nphi_p);
          TVector3 p3(pv.value * stp * std::cos(phip),
                      pv.value * stp * std::sin(phip),
                      pv.value * ctp);

          genie::Interaction bound(fold_base);
          bound.SetBit(genie::kISkipProcessChk);
          bound.SetBit(genie::kISkipKinematicChk);
          bound.SetBit(genie::kIAssumeFreeNucleon);
          bound.InitStatePtr()->TgtPtr()->SetHitNucPosition(rv.value);
          cfg.nucl_model->SetMomentum3(p3);
          cfg.nucl_model->SetRemovalEnergy(removal_energy);
          if (!qel_bind_old_fermi_mover(&bound, cfg, p3,
                                        removal_energy)) {
            continue;
          }
          ++debug.states;

          const double hit_e = bound.InitState().Tgt().HitNucP4().E();
          debug.hit_e_min = std::min(debug.hit_e_min, hit_e);
          debug.hit_e_max = std::max(debug.hit_e_max, hit_e);
          const double sqrt_s = bound.InitState().CMEnergy();
          debug.sqrt_s_min = std::min(debug.sqrt_s_min, sqrt_s);
          debug.sqrt_s_max = std::max(debug.sqrt_s_max, sqrt_s);

          const genie::Range1D_t q2lim =
              bound.PhaseSpace().Limits(genie::kKVQ2);
          if (!std::isfinite(q2lim.min) || !std::isfinite(q2lim.max) ||
              q2lim.max <= q2lim.min) {
            continue;
          }
          const double q2min = q2lim.min + genie::controls::kASmallNum;
          const double q2max = q2lim.max - genie::controls::kASmallNum;
          if (q2max <= q2min) continue;
          const double dq2 =
              (q2max - q2min) / static_cast<double>(opt.qel_fold_ncos0);
          const double dphi =
              2.0 * kPi / static_cast<double>(opt.qel_fold_nphi0);

          const double state_weight =
              rv.weight * pv.weight /
              static_cast<double>(opt.qel_fold_ncosth_p * opt.qel_fold_nphi_p);

          for (int iq2 = 0; iq2 < opt.qel_fold_ncos0; ++iq2) {
            const double q2 = q2min + (static_cast<double>(iq2) + 0.5) * dq2;
            for (int iphi = 0; iphi < opt.qel_fold_nphi0; ++iphi) {
              ++debug.q2_cells;
              const double phi =
                  2.0 * kPi * (static_cast<double>(iphi) + 0.5) /
                  static_cast<double>(opt.qel_fold_nphi0);
              genie::Interaction generated(bound);
              const double xsec = rosenbluth_qel_old_generator_q2_xsec(
                  model, &generated, cfg, q2, phi,
                  opt.qel_fold_debug ? &debug : nullptr);
              if (xsec <= 0.0) continue;
              if (!qel_lab_lepton_in_bin(generated, center, opt)) continue;
              ++debug.bin_hits;
              folded += state_weight * xsec * dq2 * dphi / (2.0 * kPi);
            }
          }
        }
      }
    }
  }
  if (opt.qel_fold_debug) {
    std::cerr << "qel-fold-debug method=grid"
              << " density=generator-q2"
              << " component=" << component_name(base_interaction, model)
              << " states=" << debug.states
              << " q2_cells=" << debug.q2_cells
              << " kin_ok=" << debug.kin_ok
              << " q2_ok=" << debug.q2_ok
              << " xsec_ok=" << debug.xsec_ok
              << " pauli_ok=" << debug.pauli_ok
              << " pauli_blocked=" << debug.pauli_blocked
              << " bin_hits=" << debug.bin_hits
              << " hitE_min=" << debug.hit_e_min
              << " hitE_max=" << debug.hit_e_max
              << " sqrtS_min=" << debug.sqrt_s_min
              << " sqrtS_max=" << debug.sqrt_s_max
              << " folded=" << folded << "\n";
  }
  return folded / (opt.qel_bin_width_energy * opt.qel_bin_width_costh);
}

double rosenbluth_qel_lattice_fold_density_generator_q2(
    const genie::XSecAlgorithmI &model,
    const genie::Interaction &base_interaction,
    const std::map<Var, double> &values,
    const Options &opt)
{
  if (!qel_bin_fold_supports_diff(opt.diff_vars)) {
    throw std::runtime_error(
        "--qel-bin-fold supports d2 in Eprime,costheta_l or omega,costheta_l only");
  }

  genie::Interaction center_interaction(base_interaction);
  const KinePoint center = solve_kinematics(values, center_interaction);
  genie::Interaction fold_base(base_interaction);
  set_probe_p4(fold_base.InitStatePtr(), center.E);

  const QELFoldConfig cfg = qel_fold_config(fold_base);
  const std::vector<WeightedValue> radii =
      qel_radius_grid(fold_base, cfg, opt);
  const double removal_energy = qel_fold_removal_energy(fold_base, opt);

  std::vector<std::vector<WeightedValue>> pgrids;
  pgrids.reserve(radii.size());
  for (const WeightedValue &rv : radii) {
    pgrids.push_back(qel_momentum_grid(fold_base, cfg, rv.value, opt));
  }

  double folded = 0.0;
  QELFoldDebugCounters debug;
  for (int is = 0; is < opt.qel_fold_samples; ++is) {
    const std::uint64_t sample = static_cast<std::uint64_t>(is);
    const std::size_t ir = select_weighted_index(radii, qel_lattice_u(sample, 0));
    const WeightedValue &rv = radii[ir];
    const std::vector<WeightedValue> &pgrid = pgrids[ir];
    const std::size_t ip = select_weighted_index(pgrid, qel_lattice_u(sample, 1));
    const WeightedValue &pv = pgrid[ip];

    const double ctp = -1.0 + 2.0 * qel_lattice_u(sample, 2);
    const double stp = std::sqrt(std::max(0.0, 1.0 - ctp * ctp));
    const double phip = 2.0 * kPi * qel_lattice_u(sample, 3);
    TVector3 p3(pv.value * stp * std::cos(phip),
                pv.value * stp * std::sin(phip),
                pv.value * ctp);

    genie::Interaction bound(fold_base);
    bound.SetBit(genie::kISkipProcessChk);
    bound.SetBit(genie::kISkipKinematicChk);
    bound.SetBit(genie::kIAssumeFreeNucleon);
    bound.InitStatePtr()->TgtPtr()->SetHitNucPosition(rv.value);
    cfg.nucl_model->SetMomentum3(p3);
    cfg.nucl_model->SetRemovalEnergy(removal_energy);
    if (!qel_bind_old_fermi_mover(&bound, cfg, p3, removal_energy)) {
      continue;
    }
    ++debug.states;

    const double hit_e = bound.InitState().Tgt().HitNucP4().E();
    debug.hit_e_min = std::min(debug.hit_e_min, hit_e);
    debug.hit_e_max = std::max(debug.hit_e_max, hit_e);
    const double sqrt_s = bound.InitState().CMEnergy();
    debug.sqrt_s_min = std::min(debug.sqrt_s_min, sqrt_s);
    debug.sqrt_s_max = std::max(debug.sqrt_s_max, sqrt_s);

    const genie::Range1D_t q2lim = bound.PhaseSpace().Limits(genie::kKVQ2);
    if (!std::isfinite(q2lim.min) || !std::isfinite(q2lim.max) ||
        q2lim.max <= q2lim.min) {
      continue;
    }
    const double q2min = q2lim.min + genie::controls::kASmallNum;
    const double q2max = q2lim.max - genie::controls::kASmallNum;
    if (q2max <= q2min) continue;

    ++debug.q2_cells;
    const double q2 = q2min + (q2max - q2min) * qel_lattice_u(sample, 4);
    const double phi = 2.0 * kPi * qel_lattice_u(sample, 5);
    genie::Interaction generated(bound);
    const double xsec = rosenbluth_qel_old_generator_q2_xsec(
        model, &generated, cfg, q2, phi,
        opt.qel_fold_debug ? &debug : nullptr);
    if (xsec <= 0.0) continue;
    if (!qel_lab_lepton_in_bin(generated, center, opt)) continue;
    ++debug.bin_hits;

    folded += xsec * (q2max - q2min);
  }

  folded /= static_cast<double>(opt.qel_fold_samples);
  if (opt.qel_fold_debug) {
    std::cerr << "qel-fold-debug method=lattice"
              << " density=generator-q2"
              << " component=" << component_name(base_interaction, model)
              << " samples=" << opt.qel_fold_samples
              << " states=" << debug.states
              << " q2_cells=" << debug.q2_cells
              << " kin_ok=" << debug.kin_ok
              << " q2_ok=" << debug.q2_ok
              << " xsec_ok=" << debug.xsec_ok
              << " pauli_ok=" << debug.pauli_ok
              << " pauli_blocked=" << debug.pauli_blocked
              << " bin_hits=" << debug.bin_hits
              << " hitE_min=" << debug.hit_e_min
              << " hitE_max=" << debug.hit_e_max
              << " sqrtS_min=" << debug.sqrt_s_min
              << " sqrtS_max=" << debug.sqrt_s_max
              << " folded=" << folded << "\n";
  }
  return folded / (opt.qel_bin_width_energy * opt.qel_bin_width_costh);
}

double rosenbluth_qel_grid_fold_density(
    const genie::XSecAlgorithmI &model,
    const genie::Interaction &base_interaction,
    const std::map<Var, double> &values,
    const Options &opt)
{
  if (!qel_bin_fold_supports_diff(opt.diff_vars)) {
    throw std::runtime_error(
        "--qel-bin-fold supports d2 in Eprime,costheta_l or omega,costheta_l only");
  }

  genie::Interaction center_interaction(base_interaction);
  const KinePoint center = solve_kinematics(values, center_interaction);
  genie::Interaction fold_base(base_interaction);
  set_probe_p4(fold_base.InitStatePtr(), center.E);

  const QELFoldConfig cfg = qel_fold_config(fold_base);
  const std::vector<WeightedValue> radii =
      qel_radius_grid(fold_base, cfg, opt);
  const double removal_energy = qel_fold_removal_energy(fold_base, opt);

  double folded = 0.0;
  QELFoldDebugCounters debug;
  for (const WeightedValue &rv : radii) {
    std::vector<WeightedValue> pgrid =
        qel_momentum_grid(fold_base, cfg, rv.value, opt);
    for (const WeightedValue &pv : pgrid) {
      for (int icp = 0; icp < opt.qel_fold_ncosth_p; ++icp) {
        const double ctp = -1.0 + 2.0 * (static_cast<double>(icp) + 0.5) /
                                      static_cast<double>(opt.qel_fold_ncosth_p);
        const double stp = std::sqrt(std::max(0.0, 1.0 - ctp * ctp));
        for (int ipp = 0; ipp < opt.qel_fold_nphi_p; ++ipp) {
          const double phip = 2.0 * kPi * (static_cast<double>(ipp) + 0.5) /
                              static_cast<double>(opt.qel_fold_nphi_p);
          TVector3 p3(pv.value * stp * std::cos(phip),
                      pv.value * stp * std::sin(phip),
                      pv.value * ctp);

          genie::Interaction bound(fold_base);
          bound.SetBit(genie::kISkipProcessChk);
          bound.SetBit(genie::kISkipKinematicChk);
          bound.InitStatePtr()->TgtPtr()->SetHitNucPosition(rv.value);
          cfg.nucl_model->SetMomentum3(p3);
          cfg.nucl_model->SetRemovalEnergy(removal_energy);
          double Eb = removal_energy;
          genie::utils::BindHitNucleon(bound, *cfg.nucl_model, Eb,
                                       cfg.binding_mode);
          ++debug.states;

          const double cos0_max =
              std::min(1.0, genie::utils::CosTheta0Max(bound));
          debug.cos0_max_min = std::min(debug.cos0_max_min, cos0_max);
          debug.cos0_max_max = std::max(debug.cos0_max_max, cos0_max);
          const double hit_e = bound.InitState().Tgt().HitNucP4().E();
          debug.hit_e_min = std::min(debug.hit_e_min, hit_e);
          debug.hit_e_max = std::max(debug.hit_e_max, hit_e);
          const double sqrt_s = bound.InitState().CMEnergy();
          debug.sqrt_s_min = std::min(debug.sqrt_s_min, sqrt_s);
          debug.sqrt_s_max = std::max(debug.sqrt_s_max, sqrt_s);
          if (!std::isfinite(cos0_max) || cos0_max <= -1.0) continue;
          const double dcos0 =
              (cos0_max + 1.0) / static_cast<double>(opt.qel_fold_ncos0);
          const double dphi0 = 2.0 * kPi /
                               static_cast<double>(opt.qel_fold_nphi0);

          const double state_weight =
              rv.weight * pv.weight /
              static_cast<double>(opt.qel_fold_ncosth_p * opt.qel_fold_nphi_p);

          for (int ic0 = 0; ic0 < opt.qel_fold_ncos0; ++ic0) {
            const double cos0 = -1.0 + (static_cast<double>(ic0) + 0.5) * dcos0;
            for (int ip0 = 0; ip0 < opt.qel_fold_nphi0; ++ip0) {
              ++debug.cos0_cells;
              const double phi0 =
                  2.0 * kPi * (static_cast<double>(ip0) + 0.5) /
                  static_cast<double>(opt.qel_fold_nphi0);
              const double xsec = rosenbluth_qel_evgen_xsec(
                  model, bound, cos0, phi0, cos0_max, cfg, opt,
                  opt.qel_fold_debug ? &debug : nullptr);
              if (xsec <= 0.0) continue;

              genie::Interaction generated(bound);
              if (!qel_set_evgen_kinematics(&generated, cos0, phi0,
                                            cfg.min_angle_em)) {
                continue;
              }
              if (!qel_lab_lepton_in_bin(generated, center, opt)) continue;
              ++debug.bin_hits;
              folded += state_weight * xsec * dcos0 * dphi0;
            }
          }
        }
      }
    }
  }
  if (opt.qel_fold_debug) {
    std::cerr << "qel-fold-debug component=" << component_name(base_interaction, model)
              << " states=" << debug.states
              << " cells=" << debug.cos0_cells
              << " kin_ok=" << debug.kin_ok
              << " q2_ok=" << debug.q2_ok
              << " jac_ok=" << debug.jac_ok
              << " factor_ok=" << debug.factor_ok
              << " xsec_ok=" << debug.xsec_ok
              << " bin_hits=" << debug.bin_hits
              << " cos0_max_min=" << debug.cos0_max_min
              << " cos0_max_max=" << debug.cos0_max_max
              << " hitE_min=" << debug.hit_e_min
              << " hitE_max=" << debug.hit_e_max
              << " sqrtS_min=" << debug.sqrt_s_min
              << " sqrtS_max=" << debug.sqrt_s_max
              << " folded=" << folded << "\n";
  }
  return folded / (opt.qel_bin_width_energy * opt.qel_bin_width_costh);
}

double rosenbluth_qel_lattice_fold_density(
    const genie::XSecAlgorithmI &model,
    const genie::Interaction &base_interaction,
    const std::map<Var, double> &values,
    const Options &opt)
{
  if (!qel_bin_fold_supports_diff(opt.diff_vars)) {
    throw std::runtime_error(
        "--qel-bin-fold supports d2 in Eprime,costheta_l or omega,costheta_l only");
  }

  genie::Interaction center_interaction(base_interaction);
  const KinePoint center = solve_kinematics(values, center_interaction);
  genie::Interaction fold_base(base_interaction);
  set_probe_p4(fold_base.InitStatePtr(), center.E);

  const QELFoldConfig cfg = qel_fold_config(fold_base);
  const std::vector<WeightedValue> radii =
      qel_radius_grid(fold_base, cfg, opt);
  const double removal_energy = qel_fold_removal_energy(fold_base, opt);

  std::vector<std::vector<WeightedValue>> pgrids;
  pgrids.reserve(radii.size());
  for (const WeightedValue &rv : radii) {
    pgrids.push_back(qel_momentum_grid(fold_base, cfg, rv.value, opt));
  }

  double folded = 0.0;
  QELFoldDebugCounters debug;
  for (int is = 0; is < opt.qel_fold_samples; ++is) {
    const std::uint64_t sample = static_cast<std::uint64_t>(is);
    const std::size_t ir = select_weighted_index(radii, qel_lattice_u(sample, 0));
    const WeightedValue &rv = radii[ir];
    const std::vector<WeightedValue> &pgrid = pgrids[ir];
    const std::size_t ip = select_weighted_index(pgrid, qel_lattice_u(sample, 1));
    const WeightedValue &pv = pgrid[ip];

    const double ctp = -1.0 + 2.0 * qel_lattice_u(sample, 2);
    const double stp = std::sqrt(std::max(0.0, 1.0 - ctp * ctp));
    const double phip = 2.0 * kPi * qel_lattice_u(sample, 3);
    TVector3 p3(pv.value * stp * std::cos(phip),
                pv.value * stp * std::sin(phip),
                pv.value * ctp);

    genie::Interaction bound(fold_base);
    bound.SetBit(genie::kISkipProcessChk);
    bound.SetBit(genie::kISkipKinematicChk);
    bound.InitStatePtr()->TgtPtr()->SetHitNucPosition(rv.value);
    cfg.nucl_model->SetMomentum3(p3);
    cfg.nucl_model->SetRemovalEnergy(removal_energy);
    double Eb = removal_energy;
    genie::utils::BindHitNucleon(bound, *cfg.nucl_model, Eb,
                                 cfg.binding_mode);
    ++debug.states;

    const double cos0_max = std::min(1.0, genie::utils::CosTheta0Max(bound));
    debug.cos0_max_min = std::min(debug.cos0_max_min, cos0_max);
    debug.cos0_max_max = std::max(debug.cos0_max_max, cos0_max);
    const double hit_e = bound.InitState().Tgt().HitNucP4().E();
    debug.hit_e_min = std::min(debug.hit_e_min, hit_e);
    debug.hit_e_max = std::max(debug.hit_e_max, hit_e);
    const double sqrt_s = bound.InitState().CMEnergy();
    debug.sqrt_s_min = std::min(debug.sqrt_s_min, sqrt_s);
    debug.sqrt_s_max = std::max(debug.sqrt_s_max, sqrt_s);
    if (!std::isfinite(cos0_max) || cos0_max <= -1.0) continue;

    ++debug.cos0_cells;
    const double cos0 = -1.0 + (cos0_max + 1.0) * qel_lattice_u(sample, 4);
    const double phi0 = 2.0 * kPi * qel_lattice_u(sample, 5);
    const double xsec = rosenbluth_qel_evgen_xsec(
        model, bound, cos0, phi0, cos0_max, cfg, opt,
        opt.qel_fold_debug ? &debug : nullptr);
    if (xsec <= 0.0) continue;

    genie::Interaction generated(bound);
    if (!qel_set_evgen_kinematics(&generated, cos0, phi0, cfg.min_angle_em)) {
      continue;
    }
    if (!qel_lab_lepton_in_bin(generated, center, opt)) continue;
    ++debug.bin_hits;

    folded += xsec * (cos0_max + 1.0) * 2.0 * kPi;
  }

  folded /= static_cast<double>(opt.qel_fold_samples);
  if (opt.qel_fold_debug) {
    std::cerr << "qel-fold-debug method=lattice"
              << " component=" << component_name(base_interaction, model)
              << " samples=" << opt.qel_fold_samples
              << " states=" << debug.states
              << " cells=" << debug.cos0_cells
              << " kin_ok=" << debug.kin_ok
              << " q2_ok=" << debug.q2_ok
              << " jac_ok=" << debug.jac_ok
              << " factor_ok=" << debug.factor_ok
              << " xsec_ok=" << debug.xsec_ok
              << " bin_hits=" << debug.bin_hits
              << " cos0_max_min=" << debug.cos0_max_min
              << " cos0_max_max=" << debug.cos0_max_max
              << " hitE_min=" << debug.hit_e_min
              << " hitE_max=" << debug.hit_e_max
              << " sqrtS_min=" << debug.sqrt_s_min
              << " sqrtS_max=" << debug.sqrt_s_max
              << " folded=" << folded << "\n";
  }
  return folded / (opt.qel_bin_width_energy * opt.qel_bin_width_costh);
}

double rosenbluth_qel_bin_fold_density(
    const genie::XSecAlgorithmI &model,
    const genie::Interaction &base_interaction,
    const std::map<Var, double> &values,
    const Options &opt)
{
  if (opt.qel_fold_density == "exact-theta") {
    return rosenbluth_qel_exact_theta_density(
        model, base_interaction, values, opt);
  }
  if (opt.qel_fold_density == "generator-q2") {
    if (opt.qel_fold_method == "lattice") {
      const std::optional<double> cached =
          rosenbluth_qel_lattice_scan_cache_density_generator_q2(
              model, base_interaction, values, opt);
      if (cached) return *cached;
    }
    if (opt.qel_fold_method == "grid") {
      return rosenbluth_qel_grid_fold_density_generator_q2(
          model, base_interaction, values, opt);
    }
    return rosenbluth_qel_lattice_fold_density_generator_q2(
        model, base_interaction, values, opt);
  }
  if (opt.qel_fold_method == "grid") {
    return rosenbluth_qel_grid_fold_density(model, base_interaction, values, opt);
  }
  return rosenbluth_qel_lattice_fold_density(model, base_interaction, values, opt);
}

std::vector<double> native_coords(NativePS native, const KinePoint &kp)
{
  switch (native) {
    case NativePS::WQ2: return {kp.W, kp.Q2};
    case NativePS::TlCtl: return {kp.Tl, kp.costheta_l};
    case NativePS::XY: return {kp.xB, kp.y};
    case NativePS::Q2: return {kp.Q2};
    case NativePS::W: return {kp.W};
    case NativePS::Auto: break;
  }
  throw std::runtime_error("Internal error: unresolved native coords");
}

std::vector<double> native_coords_for_values(
    NativePS native, const std::map<Var, double> &values,
    const genie::Interaction &interaction)
{
  return native_coords(native, solve_kinematics(values, interaction));
}

double finite_step(double x, double rel_step)
{
  return rel_step * std::max(1.0, std::fabs(x));
}

double derivative_1d_native_wrt_var(NativePS native,
                                    const std::map<Var, double> &values,
                                    Var var,
                                    const genie::Interaction &interaction,
                                    double rel_step)
{
  if (var == Var::E) {
    throw std::runtime_error("E is a conditioning variable, not a differential phase-space variable");
  }
  const double x0 = get(values, var);
  const double h = finite_step(x0, rel_step);
  auto plus = values;
  auto minus = values;
  plus[var] = x0 + h;
  minus[var] = x0 - h;
  try {
    const double cp = native_coords_for_values(native, plus, interaction).at(0);
    const double cm = native_coords_for_values(native, minus, interaction).at(0);
    return (cp - cm) / (2.0 * h);
  } catch (const OutOfPhaseSpaceError &) {
    const double c0 = native_coords_for_values(native, values, interaction).at(0);
    try {
      const double cp = native_coords_for_values(native, plus, interaction).at(0);
      return (cp - c0) / h;
    } catch (const OutOfPhaseSpaceError &) {
      const double cm = native_coords_for_values(native, minus, interaction).at(0);
      return (c0 - cm) / h;
    }
  }
}

double jacobian_native_wrt_vars(NativePS native,
                                const std::map<Var, double> &values,
                                const std::vector<Var> &vars,
                                const genie::Interaction &interaction,
                                double rel_step)
{
  if (vars.size() == 1) {
    return std::fabs(derivative_1d_native_wrt_var(native, values, vars[0],
                                                  interaction, rel_step));
  }
  if (vars.size() != 2) {
    throw std::runtime_error("Only one- and two-dimensional densities are supported");
  }
  if (vars[0] == Var::E || vars[1] == Var::E) {
    throw std::runtime_error("E is a conditioning variable, not a differential phase-space variable");
  }

  double deriv[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  for (int j = 0; j < 2; ++j) {
    const Var var = vars[j];
    const double x0 = get(values, var);
    const double h = finite_step(x0, rel_step);
    auto plus = values;
    auto minus = values;
    plus[var] = x0 + h;
    minus[var] = x0 - h;
    std::vector<double> cp;
    std::vector<double> cm;
    bool central = false;
    try {
      cp = native_coords_for_values(native, plus, interaction);
      cm = native_coords_for_values(native, minus, interaction);
      central = true;
    } catch (const OutOfPhaseSpaceError &) {
      const std::vector<double> c0 = native_coords_for_values(native, values, interaction);
      try {
        cp = native_coords_for_values(native, plus, interaction);
        cm = c0;
        central = false;
      } catch (const OutOfPhaseSpaceError &) {
        cp = c0;
        cm = native_coords_for_values(native, minus, interaction);
        central = false;
      }
    }
    if (cp.size() != 2 || cm.size() != 2) {
      throw std::runtime_error("Native phase space is not two-dimensional");
    }
    const double denom = central ? 2.0 * h : h;
    for (int row = 0; row < 2; ++row) deriv[row][j] = (cp[row] - cm[row]) / denom;
  }
  return std::fabs(deriv[0][0] * deriv[1][1] - deriv[0][1] * deriv[1][0]);
}

bool initial_state_fold_supports_diff(const std::vector<Var> &diff_vars)
{
  if (diff_vars.size() != 2) return false;
  bool has_energy = false;
  bool has_costh = false;
  for (Var var : diff_vars) {
    has_energy = has_energy || var == Var::Eprime || var == Var::Omega;
    has_costh = has_costh || var == Var::CosThetaL;
  }
  return has_energy && has_costh;
}

double raw_density_for_dims(const genie::XSecAlgorithmI &model,
                            const genie::Interaction &base_interaction,
                            const std::map<Var, double> &values,
                            const std::vector<Var> &diff_vars,
                            const Options &opt);

bool initial_state_fold_applies(const genie::Interaction &interaction,
                                const Options &opt)
{
  const genie::ProcessInfo &process = interaction.ProcInfo();
  if (process.IsResonant()) return opt.initial_state_fold.count("res") != 0;
  if (process.IsMEC()) return opt.initial_state_fold.count("mec") != 0;
  if (process.IsDeepInelastic()) return opt.initial_state_fold.count("dis") != 0;
  return false;
}

bool folded_wq2(const genie::Interaction &interaction,
                 const std::map<Var, double> &values,
                 double &W, double &Q2, KinePoint *resolved = nullptr,
                 double phi_lab = 0.0)
{
  KinePoint kp;
  try {
    kp = solve_kinematics(values, interaction);
  } catch (const OutOfPhaseSpaceError &) {
    return false;
  }

  const double mi = particle_mass(interaction.InitState().ProbePdg());
  const double mf = final_lepton_mass(interaction);
  const double pi = std::sqrt(std::max(0.0, kp.E * kp.E - mi * mi));
  const double pf = std::sqrt(std::max(0.0, kp.Eprime * kp.Eprime - mf * mf));
  const double sintheta =
      std::sqrt(std::max(0.0, 1.0 - kp.costheta_l * kp.costheta_l));
  TLorentzVector probe(0.0, 0.0, pi, kp.E);
  TLorentzVector lepton(pf * sintheta * std::cos(phi_lab),
                        pf * sintheta * std::sin(phi_lab),
                        pf * kp.costheta_l, kp.Eprime);
  const TLorentzVector q4 = probe - lepton;
  const TLorentzVector had4 = interaction.InitState().Tgt().HitNucP4() + q4;
  Q2 = -q4.Mag2();
  const double W2 = had4.Mag2();
  if (!std::isfinite(Q2) || !std::isfinite(W2) || Q2 < 0.0 || W2 <= 0.0) {
    return false;
  }
  W = std::sqrt(W2);
  if (resolved) *resolved = kp;
  return true;
}

bool apply_folded_kinematics(genie::Interaction *interaction,
                             const std::map<Var, double> &values,
                             double phi_lab = 0.0)
{
  double W = 0.0;
  double Q2 = 0.0;
  KinePoint kp;
  if (!folded_wq2(*interaction, values, W, Q2, &kp, phi_lab)) return false;
  set_probe_p4(interaction->InitStatePtr(), kp.E);

  const double E_hit_rest = interaction->InitState().ProbeE(genie::kRfHitNucRest);
  const double M_hit = interaction->InitState().Tgt().HitNucP4().M();
  if (!std::isfinite(E_hit_rest) || !std::isfinite(M_hit) ||
      E_hit_rest <= 0.0 || M_hit <= 0.0) {
    return false;
  }
  double x = 0.0;
  double y = 0.0;
  genie::utils::kinematics::WQ2toXY(E_hit_rest, M_hit, W, Q2, x, y);
  if (!std::isfinite(x) || !std::isfinite(y)) return false;

  const double mf = final_lepton_mass(*interaction);
  const double pf = std::sqrt(std::max(0.0, kp.Eprime * kp.Eprime - mf * mf));
  const double sintheta =
      std::sqrt(std::max(0.0, 1.0 - kp.costheta_l * kp.costheta_l));
  TLorentzVector lepton(pf * sintheta * std::cos(phi_lab),
                        pf * sintheta * std::sin(phi_lab),
                        pf * kp.costheta_l, kp.Eprime);
  std::unique_ptr<TLorentzVector> probe(
      interaction->InitState().GetProbeP4(genie::kRfLab));
  const TLorentzVector had4 =
      interaction->InitState().Tgt().HitNucP4() + *probe - lepton;

  genie::Kinematics *k = interaction->KinePtr();
  k->ClearRunningValues();
  k->SetW(W);
  k->SetW(W, true);
  k->SetQ2(Q2);
  k->SetQ2(Q2, true);
  k->Setq2(-Q2);
  k->Setq2(-Q2, true);
  k->Setx(x);
  k->Setx(x, true);
  k->Sety(y);
  k->Sety(y, true);
  k->SetKV(genie::kKVTl, kp.Tl);
  k->SetKV(genie::kKVctl, kp.costheta_l);
  k->SetKV(genie::kKVv, kp.omega);
  k->SetKV(genie::kKVQ0, kp.omega);
  k->SetKV(genie::kKVQ3, kp.q3);
  k->SetFSLeptonP4(lepton);
  k->SetHadSystP4(had4);
  return true;
}

double folded_wq2_jacobian(const genie::Interaction &interaction,
                           const std::map<Var, double> &values,
                           const std::vector<Var> &diff_vars,
                           double rel_step)
{
  // XSec(W,Q2) is integrated over the outgoing-lepton azimuth in the hit-
  // nucleon rest frame.  For a moving hit state, the correct lab density needs
  // the full map (E', cos(theta), phi_lab) -> (W, Q2, phi_hit), not merely the
  // upper-left 2x2 W,Q2 Jacobian.
  auto coords = [&](const std::map<Var, double> &point, double phi_lab,
                    double out[3]) -> bool {
    genie::Interaction state(interaction);
    if (!apply_folded_kinematics(&state, point, phi_lab)) return false;
    out[0] = state.Kine().W();
    out[1] = state.Kine().Q2();
    out[2] = scattering_phi_hit_rest(state);
    return std::isfinite(out[0]) && std::isfinite(out[1]) &&
           std::isfinite(out[2]);
  };

  double deriv[3][3] = {{0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
  for (int j = 0; j < 3; ++j) {
    auto plus = values;
    auto minus = values;
    double phi_plus = 0.0;
    double phi_minus = 0.0;
    double h = 1e-5;
    if (j < 2) {
      const Var var = diff_vars[j];
      const double x0 = get(values, var);
      h = finite_step(x0, rel_step);
      plus[var] = x0 + h;
      minus[var] = x0 - h;
    } else {
      phi_plus = h;
      phi_minus = -h;
    }
    double cp[3] = {0.0, 0.0, 0.0};
    double cm[3] = {0.0, 0.0, 0.0};
    if (!coords(plus, phi_plus, cp) || !coords(minus, phi_minus, cm)) {
      return 0.0;
    }
    deriv[0][j] = (cp[0] - cm[0]) / (2.0 * h);
    deriv[1][j] = (cp[1] - cm[1]) / (2.0 * h);
    deriv[2][j] = wrap_angle(cp[2] - cm[2]) / (2.0 * h);
  }
  const double det =
      deriv[0][0] * (deriv[1][1] * deriv[2][2] -
                     deriv[1][2] * deriv[2][1]) -
      deriv[0][1] * (deriv[1][0] * deriv[2][2] -
                     deriv[1][2] * deriv[2][0]) +
      deriv[0][2] * (deriv[1][0] * deriv[2][1] -
                     deriv[1][1] * deriv[2][0]);
  return std::fabs(det);
}

TVector3 isotropic_vector(double magnitude, double u_costh, double u_phi)
{
  const double costh = -1.0 + 2.0 * u_costh;
  const double sinth = std::sqrt(std::max(0.0, 1.0 - costh * costh));
  const double phi = 2.0 * kPi * u_phi;
  return TVector3(magnitude * sinth * std::cos(phi),
                  magnitude * sinth * std::sin(phi), magnitude * costh);
}

double folded_state_density(const genie::XSecAlgorithmI &model,
                            genie::Interaction &interaction,
                            const std::map<Var, double> &values,
                            const std::vector<Var> &diff_vars,
                            const Options &opt,
                            bool *below_threshold = nullptr,
                            bool *outside_phase_space = nullptr)
{
  if (below_threshold) *below_threshold = false;
  if (outside_phase_space) *outside_phase_space = false;
  if (!apply_folded_kinematics(&interaction, values)) return 0.0;
  if (opt.initial_state_fold_event_phase_space) {
    interaction.ResetBit(genie::kISkipKinematicChk);
    if (!interaction.PhaseSpace().IsAboveThreshold()) {
      if (below_threshold) *below_threshold = true;
      return 0.0;
    }
    if (!interaction.PhaseSpace().IsAllowed()) {
      if (outside_phase_space) *outside_phase_space = true;
      return 0.0;
    }
  } else {
    interaction.SetBit(genie::kISkipKinematicChk);
  }
  if (!model.ValidProcess(&interaction)) return 0.0;
  const double value = model.XSec(&interaction, genie::kPSWQ2fE);
  if (!std::isfinite(value) || value <= 0.0) return 0.0;
  const double jac = folded_wq2_jacobian(interaction, values, diff_vars,
                                         opt.jac_step);
  if (!std::isfinite(jac) || jac <= 0.0) return 0.0;
  return value * jac;
}

std::pair<int, int> mec_cluster_constituents(int cluster_pdg)
{
  if (cluster_pdg == genie::kPdgClusterNN) {
    return {genie::kPdgNeutron, genie::kPdgNeutron};
  }
  if (cluster_pdg == genie::kPdgClusterNP) {
    return {genie::kPdgNeutron, genie::kPdgProton};
  }
  if (cluster_pdg == genie::kPdgClusterPP) {
    return {genie::kPdgProton, genie::kPdgProton};
  }
  throw std::runtime_error("MEC initial-state fold found unknown cluster PDG " +
                           std::to_string(cluster_pdg));
}

double initial_state_fold_density(const genie::XSecAlgorithmI &model,
                                  const genie::Interaction &base_interaction,
                                  const std::map<Var, double> &values,
                                  const std::vector<Var> &diff_vars,
                                  const Options &opt)
{
  if (!initial_state_fold_supports_diff(diff_vars)) {
    throw std::runtime_error(
        "--initial-state-fold supports only d2 Eprime/omega,costheta_l");
  }
  if (!base_interaction.InitState().Tgt().IsNucleus()) {
    return raw_density_for_dims(model, base_interaction, values, diff_vars, opt);
  }

  Options grid_opt = opt;
  grid_opt.qel_fold_nr = opt.initial_state_fold_nr;
  grid_opt.qel_fold_np = opt.initial_state_fold_np;
  const QELFoldConfig cfg = qel_fold_config(base_interaction);
  double sum = 0.0;
  long accepted = 0;
  long below_threshold = 0;
  long outside_phase_space = 0;

  if (base_interaction.ProcInfo().IsMEC()) {
    const int cluster_pdg = base_interaction.InitState().Tgt().HitNucPdg();
    const auto constituents = mec_cluster_constituents(cluster_pdg);
    genie::Interaction first(base_interaction);
    genie::Interaction second(base_interaction);
    first.InitStatePtr()->TgtPtr()->SetHitNucPdg(constituents.first);
    second.InitStatePtr()->TgtPtr()->SetHitNucPdg(constituents.second);
    const std::vector<WeightedValue> pgrid1 =
        qel_momentum_grid(first, cfg, 0.0, grid_opt);
    const std::vector<WeightedValue> pgrid2 =
        qel_momentum_grid(second, cfg, 0.0, grid_opt);
    const double cluster_mass = particle_mass(cluster_pdg);

    for (int is = 0; is < opt.initial_state_fold_samples; ++is) {
      const std::uint64_t sample = static_cast<std::uint64_t>(is);
      const double p1 = pgrid1[select_weighted_index(
          pgrid1, qel_lattice_u(sample, 0))].value;
      const double p2 = pgrid2[select_weighted_index(
          pgrid2, qel_lattice_u(sample, 3))].value;
      const TVector3 cluster_p3 =
          isotropic_vector(p1, qel_lattice_u(sample, 1),
                           qel_lattice_u(sample, 2)) +
          isotropic_vector(p2, qel_lattice_u(sample, 4),
                           qel_lattice_u(sample, 5));
      genie::Interaction state(base_interaction);
      state.InitStatePtr()->TgtPtr()->SetHitNucP4(TLorentzVector(
          cluster_p3, std::sqrt(cluster_mass * cluster_mass +
                                cluster_p3.Mag2())));
      bool sample_below_threshold = false;
      bool sample_outside_phase_space = false;
      const double density = folded_state_density(
          model, state, values, diff_vars, opt, &sample_below_threshold,
          &sample_outside_phase_space);
      sum += density;
      below_threshold += sample_below_threshold ? 1 : 0;
      outside_phase_space += sample_outside_phase_space ? 1 : 0;
      if (density > 0.0) ++accepted;
    }
  } else {
    const std::vector<WeightedValue> radii =
        qel_radius_grid(base_interaction, cfg, grid_opt);
    std::vector<std::vector<WeightedValue>> pgrids;
    for (const WeightedValue &radius : radii) {
      pgrids.push_back(
          qel_momentum_grid(base_interaction, cfg, radius.value, grid_opt));
    }
    const double removal_energy =
        qel_fold_removal_energy(base_interaction, grid_opt);

    for (int is = 0; is < opt.initial_state_fold_samples; ++is) {
      const std::uint64_t sample = static_cast<std::uint64_t>(is);
      const std::size_t ir =
          select_weighted_index(radii, qel_lattice_u(sample, 0));
      const std::vector<WeightedValue> &pgrid = pgrids[ir];
      const double p = pgrid[select_weighted_index(
          pgrid, qel_lattice_u(sample, 1))].value;
      const TVector3 p3 = isotropic_vector(
          p, qel_lattice_u(sample, 2), qel_lattice_u(sample, 3));
      genie::Interaction state(base_interaction);
      if (!qel_bind_old_fermi_mover(&state, cfg, p3, removal_energy)) continue;
      bool sample_below_threshold = false;
      bool sample_outside_phase_space = false;
      const double density = folded_state_density(
          model, state, values, diff_vars, opt, &sample_below_threshold,
          &sample_outside_phase_space);
      sum += density;
      below_threshold += sample_below_threshold ? 1 : 0;
      outside_phase_space += sample_outside_phase_space ? 1 : 0;
      if (density > 0.0) ++accepted;
    }
  }

  const double folded = sum / static_cast<double>(opt.initial_state_fold_samples);
  if (opt.initial_state_fold_debug) {
    std::cerr << "initial-state-fold component="
              << component_name(base_interaction, model)
              << " samples=" << opt.initial_state_fold_samples
              << " accepted=" << accepted
              << " below-threshold=" << below_threshold
              << " outside-phase-space=" << outside_phase_space
              << " density=" << folded << "\n";
  }
  return folded;
}

double raw_density_for_dims(const genie::XSecAlgorithmI &model,
                            const genie::Interaction &base_interaction,
                            const std::map<Var, double> &values,
                            const std::vector<Var> &diff_vars,
                            const Options &opt)
{
  if (diff_vars.empty()) throw std::runtime_error("No differential variables requested");

  genie::Interaction interaction(base_interaction);
  const KinePoint kp = solve_kinematics(values, interaction);
  apply_kinematics(&interaction, kp);

  NativePS native = opt.native_ps;
  if (native == NativePS::Auto) {
    native = auto_native_ps(interaction, model, static_cast<int>(diff_vars.size()));
  }
  if (diff_vars.size() == 1 && (native == NativePS::WQ2 || native == NativePS::TlCtl || native == NativePS::XY)) {
    throw std::runtime_error(
        "d1 from a two-dimensional native model needs --integrate-over");
  }
  if (diff_vars.size() == 2 && (native == NativePS::Q2 || native == NativePS::W)) {
    throw std::runtime_error("d2 needs a two-dimensional native phase space");
  }

  if (!model.ValidProcess(&interaction)) return 0.0;
  const double native_value = model.XSec(&interaction, phase_space(native));
  if (!std::isfinite(native_value)) throw std::runtime_error("Model returned non-finite cross section");
  const double jac = jacobian_native_wrt_vars(native, values, diff_vars,
                                              interaction, opt.jac_step);
  return native_value * jac;
}

std::string shape_norm_cache_key(const genie::XSecAlgorithmI &model,
                                 const genie::Interaction &interaction,
                                 double E,
                                 const Options &opt)
{
  std::ostringstream key;
  key << std::setprecision(17)
      << model.Id().Key() << "|"
      << component_name(interaction, model) << "|"
      << "target=" << interaction.InitState().Tgt().Pdg() << "|"
      << "probe=" << interaction.InitState().ProbePdg() << "|"
      << "E=" << E << "|"
      << "mode=" << shape_norm_name(opt.shape_norm) << "|"
      << "ne=" << opt.shape_norm_ne << "|"
      << "ncth=" << opt.shape_norm_ncosth << "|"
      << "threshold=" << opt.shape_norm_auto_threshold << "|"
      << "native=" << native_name(opt.native_ps) << "|"
      << "emq2=" << (opt.em_q2_min_set ? opt.em_q2_min : -1.0);
  return key.str();
}

bool shape_norm_applies_to_component(const genie::Interaction &interaction,
                                     const genie::XSecAlgorithmI &model,
                                     const Options &opt)
{
  switch (opt.shape_norm) {
    case ShapeNormMode::Off:
      return false;
    case ShapeNormMode::EmpiricalMEC:
      return is_empirical_mec(interaction, model);
    case ShapeNormMode::Auto:
    case ShapeNormMode::All:
      return true;
  }
  return false;
}

double shape_norm_factor(const genie::XSecAlgorithmI &model,
                         const genie::Interaction &base_interaction,
                         const std::map<Var, double> &values,
                         const Options &opt)
{
  if (!shape_norm_applies_to_component(base_interaction, model, opt) ||
      !has(values, Var::E)) {
    return 1.0;
  }

  const double E = get(values, Var::E);
  static std::map<std::string, double> cache;
  const std::string key = shape_norm_cache_key(model, base_interaction, E, opt);
  auto cached = cache.find(key);
  if (cached != cache.end()) return cached->second;

  genie::Interaction norm_interaction(base_interaction);
  set_probe_p4(norm_interaction.InitStatePtr(), E);
  if (!model.ValidProcess(&norm_interaction)) {
    cache[key] = 1.0;
    return 1.0;
  }

  const double total = model.Integral(&norm_interaction);
  if (!std::isfinite(total) || total <= 0.0) {
    cache[key] = 1.0;
    return 1.0;
  }

  const double mf = final_lepton_mass(norm_interaction);
  const double eps = 1e-6;
  const double ep_min = mf + eps;
  const double ep_max = E - eps;
  if (ep_max <= ep_min) {
    throw std::runtime_error("Shape normalization has empty Eprime range");
  }

  const int ne = opt.shape_norm_ne;
  const int ncth = opt.shape_norm_ncosth;
  const double dep = (ep_max - ep_min) / static_cast<double>(ne);
  const double dcth = 2.0 / static_cast<double>(ncth);
  const std::vector<Var> density_vars = {Var::Eprime, Var::CosThetaL};

  double shape = 0.0;
  for (int ie = 0; ie < ne; ++ie) {
    const double ep = ep_min + (ie + 0.5) * dep;
    for (int ic = 0; ic < ncth; ++ic) {
      const double cth = -1.0 + (ic + 0.5) * dcth;
      std::map<Var, double> p;
      p[Var::E] = E;
      p[Var::Eprime] = ep;
      p[Var::CosThetaL] = cth;
      try {
        const double density = raw_density_for_dims(model, norm_interaction,
                                                    p, density_vars, opt);
        if (std::isfinite(density) && density > 0.0) {
          shape += density * dep * dcth;
        }
      } catch (const OutOfPhaseSpaceError &) {
        // Kinematically forbidden midpoints contribute zero.
      }
    }
  }

  if (!std::isfinite(shape) || shape <= 0.0) {
    if (opt.shape_norm == ShapeNormMode::Auto) {
      cache[key] = 1.0;
      std::cerr << "shape-norm component="
                << component_name(base_interaction, model)
                << " mode=" << shape_norm_name(opt.shape_norm)
                << " applied=no"
                << " reason=zero-shape"
                << " E=" << E
                << " total=" << total
                << " grid=" << ne << "x" << ncth << "\n";
      return 1.0;
    }
    throw std::runtime_error("Shape normalization integral is zero");
  }

  const double factor = total / shape;
  if (!std::isfinite(factor) || factor <= 0.0) {
    if (opt.shape_norm == ShapeNormMode::Auto) {
      cache[key] = 1.0;
      return 1.0;
    }
    throw std::runtime_error("Shape normalization factor is non-finite");
  }

  const bool apply =
      opt.shape_norm == ShapeNormMode::All ||
      opt.shape_norm == ShapeNormMode::EmpiricalMEC ||
      factor > opt.shape_norm_auto_threshold ||
      factor < 1.0 / opt.shape_norm_auto_threshold;
  cache[key] = apply ? factor : 1.0;
  std::cerr << "shape-norm component="
            << component_name(base_interaction, model)
            << " mode=" << shape_norm_name(opt.shape_norm)
            << " applied=" << (apply ? "yes" : "no")
            << " E=" << E
            << " total=" << total
            << " shape=" << shape
            << " factor=" << factor
            << " used_factor=" << cache[key]
            << " grid=" << ne << "x" << ncth << "\n";
  return cache[key];
}

double apply_shape_norm(const genie::XSecAlgorithmI &model,
                        const genie::Interaction &interaction,
                        const std::map<Var, double> &values,
                        const Options &opt,
                        double value)
{
  if (value == 0.0) return value;
  return value * shape_norm_factor(model, interaction, values, opt);
}

double integrate_density(const genie::XSecAlgorithmI &model,
                         const genie::Interaction &interaction,
                         const std::map<Var, double> &base_values,
                         const std::vector<Var> &density_vars,
                         const std::vector<RangeSpec> &ranges,
                         const Options &opt)
{
  double sum = 0.0;
  std::size_t physical_midpoints = 0;
  std::map<Var, double> current = base_values;
  std::function<void(size_t, double)> rec = [&](size_t idx, double weight) {
    if (idx == ranges.size()) {
      try {
        sum += weight * raw_density_for_dims(model, interaction, current,
                                             density_vars, opt);
        ++physical_midpoints;
      } catch (const OutOfPhaseSpaceError &) {
        // Integrals are over the intersection with physical phase space.
      }
      return;
    }
    const RangeSpec &r = ranges[idx];
    const double width = (r.max - r.min) / static_cast<double>(r.n);
    for (int i = 0; i < r.n; ++i) {
      current[r.var] = r.min + (static_cast<double>(i) + 0.5) * width;
      rec(idx + 1, weight * width);
    }
  };
  rec(0, 1.0);
  if (physical_midpoints == 0) {
    throw OutOfPhaseSpaceError("Integration range contains no physical midpoints");
  }
  return sum;
}

double evaluate_model(const genie::XSecAlgorithmI &model,
                      const genie::Interaction &base_interaction,
                      const std::map<Var, double> &values,
                      const Options &opt)
{
  if (opt.observable == Observable::D2) {
    if (!opt.integrate_over.empty()) {
      throw std::runtime_error("d2 does not support extra integration variables");
    }
    if (opt.qel_bin_fold && is_rosenbluth_qel(base_interaction, model)) {
      return rosenbluth_qel_bin_fold_density(model, base_interaction, values, opt);
    }
    if (initial_state_fold_applies(base_interaction, opt)) {
      const double value = initial_state_fold_density(
          model, base_interaction, values, opt.diff_vars, opt);
      return apply_shape_norm(model, base_interaction, values, opt, value);
    }
    const double value =
        raw_density_for_dims(model, base_interaction, values, opt.diff_vars, opt);
    return apply_shape_norm(model, base_interaction, values, opt, value);
  }

  if (opt.observable == Observable::D1) {
    double value = 0.0;
    if (opt.integrate_over.empty()) {
      value = raw_density_for_dims(model, base_interaction, values, opt.diff_vars, opt);
    } else {
      if (opt.integrate_over.size() != 1) {
        throw std::runtime_error("d1 projection needs exactly one --integrate-over range");
      }
      std::vector<Var> density_vars = {opt.diff_vars[0], opt.integrate_over[0].var};
      value = integrate_density(model, base_interaction, values, density_vars,
                                opt.integrate_over, opt);
    }
    return apply_shape_norm(model, base_interaction, values, opt, value);
  }

  if (opt.integrate_over.empty()) {
    genie::Interaction interaction(base_interaction);
    set_probe_p4(interaction.InitStatePtr(), get(values, Var::E));
    if (!model.ValidProcess(&interaction)) return 0.0;
    return model.Integral(&interaction);
  }
  if (opt.integrate_over.size() == 1) {
    std::vector<Var> density_vars = {opt.integrate_over[0].var};
    const double value = integrate_density(model, base_interaction, values, density_vars,
                                           opt.integrate_over, opt);
    return apply_shape_norm(model, base_interaction, values, opt, value);
  }
  if (opt.integrate_over.size() == 2) {
    std::vector<Var> density_vars = {opt.integrate_over[0].var,
                                     opt.integrate_over[1].var};
    const double value = integrate_density(model, base_interaction, values, density_vars,
                                           opt.integrate_over, opt);
    return apply_shape_norm(model, base_interaction, values, opt, value);
  }
  throw std::runtime_error("total projections support at most two --integrate-over ranges");
}

ComponentValue evaluate_component(const std::string &name,
                                  const genie::XSecAlgorithmI &model,
                                  const genie::Interaction &interaction,
                                  const std::map<Var, double> &values,
                                  const Options &opt)
{
  ComponentValue out;
  out.component = name;
  try {
    if (opt.observable == Observable::Total && opt.integrate_over.empty()) {
      if (has(values, Var::E)) out.kin_values[Var::E] = get(values, Var::E);
    } else if (opt.integrate_over.empty()) {
      genie::Interaction tmp(interaction);
      out.kin_values = output_map_from_kine(solve_kinematics(values, tmp));
    }
    out.value = evaluate_model(model, interaction, values, opt);
    if (!std::isfinite(out.value)) throw std::runtime_error("non-finite result");
  } catch (const std::exception &e) {
    if (opt.strict) throw;
    out.status = "invalid";
    out.value = 0.0;
    out.message = e.what();
  }
  return out;
}

std::vector<ComponentValue> evaluate_tune_point(const Point &point,
                                                const Options &opt)
{
  const double E = get(point.values, Var::E);
  genie::InitialState init = make_initial_state(opt, E);

  genie::GEVGDriver driver;
  driver.SetEventGeneratorList(opt.event_generator_list);
  driver.Configure(init);
  const genie::InteractionList *interactions = driver.Interactions();
  if (!interactions || interactions->empty()) {
    throw std::runtime_error("GEVGDriver returned no interactions");
  }

  std::vector<ComponentValue> components;
  for (const genie::Interaction *interaction : *interactions) {
    const genie::EventGeneratorI *generator = driver.FindGenerator(interaction);
    if (!generator || !generator->CrossSectionAlg()) continue;
    const genie::XSecAlgorithmI *model = generator->CrossSectionAlg();
    components.push_back(evaluate_component(component_name(*interaction, *model),
                                            *model, *interaction, point.values, opt));
  }
  return components;
}

std::vector<ComponentValue> evaluate_alg_point(const Point &point,
                                               const Options &opt,
                                               const genie::XSecAlgorithmI &model)
{
  const double E = get(point.values, Var::E);
  genie::InitialState init = make_initial_state(opt, E);
  std::unique_ptr<genie::InteractionList> list =
      make_algorithm_interactions(opt, init);

  std::set<std::string> seen;
  std::vector<ComponentValue> components;
  for (const genie::Interaction *interaction : *list) {
    if (opt.interaction_sum == "hit-state") {
      const std::string key = hit_state_key(*interaction);
      if (seen.count(key)) continue;
      seen.insert(key);
    }
    components.push_back(evaluate_component(component_name(*interaction, model),
                                            model, *interaction, point.values, opt));
  }
  return components;
}

void write_metadata(std::ostream &out, int argc, char **argv,
                    const Options &opt)
{
  out << "# command";
  for (int i = 0; i < argc; ++i) out << "," << csv_escape(argv[i]);
  out << "\n";
  out << "# mode," << (opt.mode == Mode::Tune ? "tune" : "alg") << "\n";
  out << "# tune," << opt.tune << "\n";
  out << "# event_generator_list," << opt.event_generator_list << "\n";
  out << "# xsec_alg," << opt.xsec_alg << "\n";
  out << "# process," << opt.process << "\n";
  out << "# interaction_sum," << opt.interaction_sum << "\n";
  out << "# probe," << opt.probe << "\n";
  out << "# target," << opt.target << "\n";
  out << "# observable," << observable_name(opt.observable) << "\n";
  out << "# diff," << diff_label(opt.diff_vars) << "\n";
  out << "# scans," << range_label(opt.scans) << "\n";
  out << "# integrate_over," << range_label(opt.integrate_over) << "\n";
  if (!opt.integrate_over.empty()) {
    out << "# integration_rule,midpoint\n";
    out << "# integration_phase_space,physical-intersection\n";
  }
  out << "# xsec_unit," << value_unit_label(opt) << "\n";
  out << "# native_phase_space," << native_name(opt.native_ps) << "\n";
  out << "# em_q2_min," << (opt.em_q2_min_set ? std::to_string(opt.em_q2_min) : "") << "\n";
  out << "# em_q2_min_status," << g_em_q2_min_status << "\n";
  if (!g_em_q2_min_message.empty()) {
    out << "# em_q2_min_message," << csv_escape(g_em_q2_min_message) << "\n";
  }
  out << "# qel_bin_fold," << (opt.qel_bin_fold ? "on" : "off") << "\n";
  out << "# initial_state_fold,";
  if (opt.initial_state_fold.empty()) {
    out << "off\n";
  } else {
    bool first = true;
    for (const std::string &process : opt.initial_state_fold) {
      if (!first) out << "+";
      out << process;
      first = false;
    }
    out << "\n";
    out << "# initial_state_fold_samples,"
        << opt.initial_state_fold_samples << "\n";
    out << "# initial_state_fold_nr," << opt.initial_state_fold_nr << "\n";
    out << "# initial_state_fold_np," << opt.initial_state_fold_np << "\n";
    out << "# initial_state_fold_event_phase_space,"
        << (opt.initial_state_fold_event_phase_space ? "on" : "off") << "\n";
  }
  out << "# shape_norm," << shape_norm_name(opt.shape_norm) << "\n";
  if (opt.shape_norm != ShapeNormMode::Off) {
    out << "# shape_norm_ne," << opt.shape_norm_ne << "\n";
    out << "# shape_norm_ncosth," << opt.shape_norm_ncosth << "\n";
    out << "# shape_norm_auto_threshold,"
        << opt.shape_norm_auto_threshold << "\n";
  }
  if (opt.qel_bin_fold) {
    out << "# qel_bin_width_energy," << opt.qel_bin_width_energy << "\n";
    out << "# qel_bin_width_costheta_l," << opt.qel_bin_width_costh << "\n";
    out << "# qel_fold_method," << opt.qel_fold_method << "\n";
    out << "# qel_fold_density," << opt.qel_fold_density << "\n";
    out << "# qel_fold_samples," << opt.qel_fold_samples << "\n";
    out << "# qel_fold_nr," << opt.qel_fold_nr << "\n";
    out << "# qel_fold_np," << opt.qel_fold_np << "\n";
    out << "# qel_fold_ncosth_p," << opt.qel_fold_ncosth_p << "\n";
    out << "# qel_fold_nphi_p," << opt.qel_fold_nphi_p << "\n";
    out << "# qel_fold_ncos0," << opt.qel_fold_ncos0 << "\n";
    out << "# qel_fold_nphi0," << opt.qel_fold_nphi0 << "\n";
    out << "# qel_fold_scan_cache,"
        << (opt.qel_fold_scan_cache ? "on" : "off") << "\n";
    out << "# qel_fold_pmax,"
        << (opt.qel_fold_pmax > 0.0 ? std::to_string(opt.qel_fold_pmax) : "auto")
        << "\n";
    out << "# qel_fold_kf,"
        << (opt.qel_fold_kf >= 0.0 ? std::to_string(opt.qel_fold_kf) : "auto")
        << "\n";
    out << "# qel_fold_removal_energy,"
        << (opt.qel_fold_removal_energy >= 0.0
                ? std::to_string(opt.qel_fold_removal_energy)
                : "auto")
        << "\n";
  }
}

std::vector<Var> output_vars()
{
  return {Var::E, Var::Eprime, Var::Q2, Var::W, Var::Omega,
          Var::Q3, Var::XB, Var::Y, Var::CosThetaL, Var::ThetaLDeg, Var::Tl};
}

void write_header(std::ostream &out)
{
  out << "row_id,component,observable,value,unit,status,message";
  for (Var v : output_vars()) out << "," << canonical_var_name(v);
  out << "\n";
}

void write_row(std::ostream &out, const Point &point,
               const ComponentValue &value, const Options &opt,
               double unit_factor)
{
  out << point.row_id << ","
      << csv_escape(value.component) << ","
      << observable_name(opt.observable) << ","
      << std::setprecision(12) << (value.value / unit_factor) << ","
      << csv_escape(value_unit_label(opt)) << ","
      << value.status << ","
      << csv_escape(value.message);

  for (Var v : output_vars()) {
    out << ",";
    auto kit = value.kin_values.find(v);
    if (kit != value.kin_values.end()) {
      out << std::setprecision(12) << kit->second;
      continue;
    }
    auto it = point.values.find(v);
    if (it != point.values.end()) out << std::setprecision(12) << it->second;
  }
  out << "\n";
}

int run(int argc, char **argv)
{
  Options opt = parse_options(argc, argv);
  configure_runtime_xml_path(argv[0]);
  if (opt.message_thresholds == "config/Messenger_whisper.xml") {
    if (const char *genie = std::getenv("GENIE")) {
      const std::filesystem::path default_thresholds =
          std::filesystem::path(genie) / "config" / "Messenger_whisper.xml";
      if (std::filesystem::exists(default_thresholds)) {
        opt.message_thresholds = default_thresholds.string();
      }
    }
  }
  const std::vector<Point> points = make_points(opt);
  validate_projection_point_values(points, opt);
  const double unit_factor = xsec_unit_factor(opt.xsec_unit);

  genie::RunOpt::Instance()->SetTuneName(opt.tune);
  genie::RunOpt::Instance()->SetEventGeneratorList(opt.event_generator_list);
  genie::RunOpt::Instance()->BuildTune();
  apply_em_q2_min_override(opt);
  if (!opt.message_thresholds.empty()) {
    genie::utils::app_init::MesgThresholds(opt.message_thresholds);
  }

  const genie::XSecAlgorithmI *direct_model = nullptr;
  if (opt.mode == Mode::Alg) direct_model = get_xsec_alg(opt.xsec_alg);

  std::filesystem::path out_path(opt.output_path);
  if (!out_path.parent_path().empty()) {
    std::filesystem::create_directories(out_path.parent_path());
  }
  std::ofstream out(opt.output_path);
  if (!out) throw std::runtime_error("Could not open output: " + opt.output_path);

  write_metadata(out, argc, argv, opt);
  write_header(out);

  for (const Point &point : points) {
    std::vector<ComponentValue> components =
        (opt.mode == Mode::Tune)
            ? evaluate_tune_point(point, opt)
            : evaluate_alg_point(point, opt, *direct_model);

    ComponentValue total;
    total.component = "total";
    bool all_components_ok = !components.empty();
    bool have_kinematics = false;
    if (components.empty()) {
      total.status = "invalid";
      total.message = "No cross-section components were evaluated";
    }
    for (const ComponentValue &component : components) {
      if (component.status != "ok") {
        all_components_ok = false;
      }
      if (component.status != "ok" && total.status == "ok") {
        total.status = component.status;
        total.message = "Component " + component.component + ": " +
                        component.message;
      } else if (component.status == "ok") {
        total.value += component.value;
      }
      if (!have_kinematics && !component.kin_values.empty()) {
        total.kin_values = component.kin_values;
        have_kinematics = true;
      }
    }
    if (all_components_ok) {
      total.status = "ok";
      total.message.clear();
    } else {
      total.value = 0.0;
    }
    write_row(out, point, total, opt, unit_factor);
    if (opt.components) {
      for (const ComponentValue &component : components) {
        write_row(out, point, component, opt, unit_factor);
      }
    }
  }

  out.close();
  std::cerr << "wrote " << opt.output_path << "\n";
  std::_Exit(0);
  return 0;
}

} // namespace

int main(int argc, char **argv)
{
  try {
    return run(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << "xsec_scan: " << e.what() << "\n\n";
    print_usage(std::cerr);
    std::_Exit(1);
  }
}
