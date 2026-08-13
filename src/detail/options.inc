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
  std::string batch;
  int target = 0;
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
  bool auto_fold = false;
  bool auto_qel_fold = false;
  bool auto_initial_state_fold = false;
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
  int jobs = 1;
  bool batch_mode = false;
  std::string batch_config_path;
  std::string message_thresholds = "config/Messenger_whisper.xml";
};

struct BatchScan {
  std::string name;
  std::optional<int> target;
  std::vector<std::string> fixed;
  std::vector<std::string> ranges;
};

struct BatchDocument {
  std::string path;
  std::vector<std::string> global_args;
  std::vector<BatchScan> scans;
};

struct Invocation {
  Options options;
  BatchDocument batch;
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
      << "  --fold auto model-aware event-generator-path folding with phase-space vetoes\n"
      << "  --jobs N|auto evaluate scan-point chunks in forked workers (default: 1)\n"
      << "  --batch-config FILE run named scans from an INI-style configuration\n"
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
      std::cout.flush();
      std::_Exit(0);
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
    } else if (arg == "--fold") {
      const std::string mode = lower(trim(need_value(arg)));
      if (mode == "auto") {
        opt.auto_fold = true;
      } else if (mode == "manual" || mode == "off" || mode == "none") {
        opt.auto_fold = false;
      } else {
        throw std::runtime_error("--fold accepts only auto, manual, or off");
      }
    } else if (arg == "--jobs") {
      const std::string value = lower(trim(need_value(arg)));
      if (value == "auto") {
        const long detected = sysconf(_SC_NPROCESSORS_ONLN);
        opt.jobs = detected > 0 ? static_cast<int>(detected) : 1;
      } else {
        opt.jobs = std::stoi(value);
      }
    } else if (arg == "--initial-state-fold") {
      for (std::string process : split(need_value(arg), ',')) {
        process = lower(trim(process));
        if (process.empty() || process == "off" || process == "none") continue;
        if (process == "auto") {
          opt.auto_fold = true;
          continue;
        }
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
  if (opt.mode == Mode::Alg && opt.auto_fold) {
    throw std::runtime_error(
        "--fold auto requires --mode tune so the resolved event-generator "
        "module chain can be inspected");
  }
  if (opt.auto_fold && opt.qel_bin_fold) {
    throw std::runtime_error(
        "--qel-bin-fold is a manual override and cannot be combined with "
        "--fold auto");
  }
  if (opt.auto_fold && !opt.initial_state_fold.empty()) {
    throw std::runtime_error(
        "an explicit --initial-state-fold process list cannot be combined "
        "with --fold auto");
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
  if (opt.jobs <= 0) {
    throw std::runtime_error("--jobs must be positive or auto");
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

std::string unquote_config_value(const std::string &value)
{
  if (value.size() >= 2 &&
      ((value.front() == '"' && value.back() == '"') ||
       (value.front() == '\'' && value.back() == '\''))) {
    return value.substr(1, value.size() - 2);
  }
  return value;
}

bool config_bool(const std::string &value, const std::string &where)
{
  const std::string normalized = lower(trim(value));
  if (normalized == "true" || normalized == "yes" || normalized == "on" ||
      normalized == "1") {
    return true;
  }
  if (normalized == "false" || normalized == "no" || normalized == "off" ||
      normalized == "0") {
    return false;
  }
  throw std::runtime_error(where + " expects a boolean value");
}

std::string config_option_name(std::string key)
{
  key = lower(trim(key));
  std::replace(key.begin(), key.end(), '_', '-');
  return key;
}

BatchDocument read_batch_config(const std::string &path)
{
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Could not open batch config: " + path);

  BatchDocument document;
  document.path = path;
  enum class Section { None, Global, Scan } section = Section::None;
  BatchScan *scan = nullptr;
  std::set<std::string> scan_names;
  const std::set<std::string> boolean_options = {
      "components", "strict", "qel-bin-fold", "qel-fold-debug",
      "qel-fold-scan-cache", "initial-state-fold-event-phase-space",
      "initial-state-fold-debug"};

  std::string line;
  int line_number = 0;
  while (std::getline(in, line)) {
    ++line_number;
    line = trim(line);
    if (line.empty() || line[0] == '#' || line[0] == ';') continue;
    const std::string where = path + ":" + std::to_string(line_number);
    if (line.front() == '[' && line.back() == ']') {
      const std::string heading = trim(line.substr(1, line.size() - 2));
      if (lower(heading) == "global") {
        section = Section::Global;
        scan = nullptr;
        continue;
      }
      const std::string prefix = "scan ";
      if (lower(heading).rfind(prefix, 0) == 0) {
        const std::string name = trim(heading.substr(prefix.size()));
        if (name.empty()) throw std::runtime_error(where + " has an empty scan name");
        if (!scan_names.insert(name).second) {
          throw std::runtime_error(where + " repeats scan name " + name);
        }
        document.scans.push_back(BatchScan{name, std::nullopt, {}, {}});
        scan = &document.scans.back();
        section = Section::Scan;
        continue;
      }
      throw std::runtime_error(where + " has unknown section [" + heading + "]");
    }

    const std::size_t equals = line.find('=');
    if (equals == std::string::npos) {
      throw std::runtime_error(where + " expects KEY = VALUE");
    }
    const std::string key = config_option_name(line.substr(0, equals));
    const std::string value =
        unquote_config_value(trim(line.substr(equals + 1)));
    if (key.empty() || value.empty()) {
      throw std::runtime_error(where + " has an empty key or value");
    }
    if (section == Section::Scan) {
      if (key == "fixed") {
        scan->fixed.push_back(value);
      } else if (key == "range" || key == "scan") {
        scan->ranges.push_back(value);
      } else if (key == "target") {
        std::size_t consumed = 0;
        long parsed = 0;
        try {
          parsed = std::stol(value, &consumed);
        } catch (const std::exception &) {
          throw std::runtime_error(where + " has an invalid target PDG");
        }
        if (consumed != value.size() || parsed <= 0 ||
            parsed > std::numeric_limits<int>::max()) {
          throw std::runtime_error(where + " has an invalid target PDG");
        }
        scan->target = static_cast<int>(parsed);
      } else {
        throw std::runtime_error(
            where + " scan sections accept only target, fixed, or range");
      }
      continue;
    }
    if (section != Section::Global) {
      throw std::runtime_error(where + " entry appears before a section");
    }
    if (key == "batch-config") {
      throw std::runtime_error(where + " cannot recursively set batch-config");
    }
    if (boolean_options.count(key)) {
      const bool enabled = config_bool(value, where);
      if (key == "qel-fold-scan-cache" && !enabled) {
        document.global_args.push_back("--no-qel-fold-scan-cache");
      } else if (enabled) {
        document.global_args.push_back("--" + key);
      }
    } else {
      document.global_args.push_back("--" + key);
      document.global_args.push_back(value);
    }
  }
  if (document.scans.empty()) {
    throw std::runtime_error("Batch config has no [scan NAME] sections: " + path);
  }
  for (const BatchScan &entry : document.scans) {
    if (entry.ranges.empty()) {
      throw std::runtime_error("Batch scan " + entry.name + " has no range entries");
    }
  }
  return document;
}

Invocation parse_invocation(int argc, char **argv)
{
  std::string batch_path;
  std::vector<std::string> cli_args;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--batch-config") {
      if (i + 1 >= argc) {
        throw std::runtime_error("Missing value for --batch-config");
      }
      batch_path = argv[++i];
      continue;
    }
    cli_args.push_back(arg);
  }

  Invocation invocation;
  if (!batch_path.empty()) invocation.batch = read_batch_config(batch_path);
  std::vector<std::string> storage;
  storage.push_back(argv[0]);
  storage.insert(storage.end(), invocation.batch.global_args.begin(),
                 invocation.batch.global_args.end());
  storage.insert(storage.end(), cli_args.begin(), cli_args.end());
  std::vector<char *> merged_argv;
  merged_argv.reserve(storage.size());
  for (std::string &arg : storage) merged_argv.push_back(arg.data());
  invocation.options =
      parse_options(static_cast<int>(merged_argv.size()), merged_argv.data());
  if (!batch_path.empty()) {
    invocation.options.batch_mode = true;
    invocation.options.batch_config_path = batch_path;
    if (!invocation.options.scans.empty() ||
        !invocation.options.points_path.empty()) {
      throw std::runtime_error(
          "--batch-config cannot be combined with global --scan or --points");
    }
  }
  return invocation;
}
