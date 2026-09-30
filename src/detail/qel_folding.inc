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
  return xsec * qel_hit_nucleon_count(interaction->InitState().Tgt());
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

#include "qel_adaptive.inc"

double rosenbluth_qel_bin_fold_density(
    const genie::XSecAlgorithmI &model,
    const genie::Interaction &base_interaction,
    const std::map<Var, double> &values,
    const Options &opt)
{
  if (!qel_bin_fold_supports_diff(opt.diff_vars)) {
    throw std::runtime_error("QE fold supports d2 Eprime/omega,costheta_l only");
  }
  if (opt.qel_fold_density == "native-q2-reference") {
    return rosenbluth_qel_native_reference(model, base_interaction, values, opt);
  }
  if (opt.qel_fold_density == "auto" || opt.qel_fold_density == "adaptive-theta") {
    return rosenbluth_qel_adaptive_theta(model, base_interaction, values, opt);
  }
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
