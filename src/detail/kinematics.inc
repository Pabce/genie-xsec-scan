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
