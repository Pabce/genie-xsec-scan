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
  if (opt.auto_initial_state_fold) return true;
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
                           double /* rel_step */)
{
  // d^3k/E = p dE dcos(theta) dphi is Lorentz invariant.  Combining that
  // measure with W^2=M^2+2M(E_i^*-E_f^*)-Q^2 gives the exact full determinant
  //
  // |d(W,Q2,phi*) / d(E'_lab,cos(theta_lab),phi_lab)|
  //       = 2 M p_i^* p'_lab / W.
  //
  // This is the same 3x3 transformation previously evaluated with six
  // finite-difference kinematic solves per Fermi-motion sample.
  if (!initial_state_fold_supports_diff(diff_vars)) return 0.0;
  const double W = interaction.Kine().W();
  const double M = interaction.InitState().Tgt().HitNucP4().M();
  const double Ei_star =
      interaction.InitState().ProbeE(genie::kRfHitNucRest);
  const double mi = particle_mass(interaction.InitState().ProbePdg());
  const double pi_star2 = Ei_star * Ei_star - mi * mi;
  const double mf = final_lepton_mass(interaction);
  const double Eprime = interaction.Kine().FSLeptonP4().E();
  const double pf_lab2 = Eprime * Eprime - mf * mf;
  if (!std::isfinite(W) || !std::isfinite(M) || !std::isfinite(Ei_star) ||
      W <= 0.0 || M <= 0.0 || pi_star2 <= 0.0 || pf_lab2 <= 0.0) {
    return 0.0;
  }
  return 2.0 * M * std::sqrt(pi_star2) * std::sqrt(pf_lab2) / W;
}

TVector3 isotropic_vector(double magnitude, double u_costh, double u_phi)
{
  const double costh = -1.0 + 2.0 * u_costh;
  const double sinth = std::sqrt(std::max(0.0, 1.0 - costh * costh));
  const double phi = 2.0 * kPi * u_phi;
  return TVector3(magnitude * sinth * std::cos(phi),
                  magnitude * sinth * std::sin(phi), magnitude * costh);
}

double prepare_folded_state(genie::Interaction &interaction,
                            const std::map<Var, double> &values,
                            const std::vector<Var> &diff_vars,
                            const Options &opt,
                            bool *below_threshold = nullptr,
                            bool *outside_phase_space = nullptr)
{
  if (below_threshold) *below_threshold = false;
  if (outside_phase_space) *outside_phase_space = false;
  if (!apply_folded_kinematics(&interaction, values)) return 0.0;
  if (opt.auto_fold || opt.initial_state_fold_event_phase_space) {
    interaction.ResetBit(genie::kISkipKinematicChk);
    if (!interaction.PhaseSpace().IsAboveThreshold()) {
      if (below_threshold) *below_threshold = true;
      return 0.0;
    }
    if (!interaction.PhaseSpace().IsAllowed()) {
      if (outside_phase_space) *outside_phase_space = true;
      return 0.0;
    }
    // The event-chain checks above are authoritative. Avoid making every
    // cross-section model repeat its own equivalent kinematic-limit lookup.
    interaction.SetBit(genie::kISkipKinematicChk);
  } else {
    interaction.SetBit(genie::kISkipKinematicChk);
  }
  const double jac = folded_wq2_jacobian(interaction, values, diff_vars,
                                         opt.jac_step);
  if (!std::isfinite(jac) || jac <= 0.0) return 0.0;
  return jac;
}

double folded_state_density(const genie::XSecAlgorithmI &model,
                            genie::Interaction &interaction,
                            const std::map<Var, double> &values,
                            const std::vector<Var> &diff_vars,
                            const Options &opt,
                            bool *below_threshold = nullptr,
                            bool *outside_phase_space = nullptr)
{
  const double jac = prepare_folded_state(
      interaction, values, diff_vars, opt, below_threshold,
      outside_phase_space);
  if (jac <= 0.0 || !model.ValidProcess(&interaction)) return 0.0;
  const double value = model.XSec(&interaction, genie::kPSWQ2fE);
  if (!std::isfinite(value) || value <= 0.0) return 0.0;
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

struct InitialStateFoldGrid {
  std::vector<WeightedValue> radii;
  std::vector<std::vector<WeightedValue>> momenta;
  double removal_energy = 0.0;
};

std::string initial_state_fold_grid_key(const genie::Interaction &interaction,
                                        const QELFoldConfig &cfg,
                                        const Options &opt)
{
  std::ostringstream key;
  key << std::setprecision(17)
      << "model=" << reinterpret_cast<std::uintptr_t>(cfg.nucl_model)
      << "|target=" << interaction.InitState().Tgt().Pdg()
      << "|hit=" << interaction.InitState().Tgt().HitNucPdg()
      << "|nr=" << opt.qel_fold_nr << "|np=" << opt.qel_fold_np
      << "|pmax=" << opt.qel_fold_pmax << "|kf=" << opt.qel_fold_kf
      << "|ermv=" << opt.qel_fold_removal_energy;
  return key.str();
}

const InitialStateFoldGrid &cached_initial_state_fold_grid(
    const genie::Interaction &interaction, const QELFoldConfig &cfg,
    const Options &opt)
{
  static std::map<std::string, InitialStateFoldGrid> cache;
  const std::string key = initial_state_fold_grid_key(interaction, cfg, opt);
  auto found = cache.find(key);
  if (found != cache.end()) return found->second;

  InitialStateFoldGrid grid;
  grid.radii = qel_radius_grid(interaction, cfg, opt);
  grid.momenta.reserve(grid.radii.size());
  for (const WeightedValue &radius : grid.radii) {
    grid.momenta.push_back(
        qel_momentum_grid(interaction, cfg, radius.value, opt));
  }
  grid.removal_energy = qel_fold_removal_energy(interaction, opt);
  return cache.emplace(key, std::move(grid)).first->second;
}

const std::vector<TVector3> &cached_initial_state_fold_momenta(
    const genie::Interaction &interaction, const QELFoldConfig &cfg,
    const Options &opt)
{
  static std::map<std::string, std::vector<TVector3>> cache;
  const std::string key = initial_state_fold_grid_key(interaction, cfg, opt) +
                          "|samples=" +
                          std::to_string(opt.initial_state_fold_samples);
  auto found = cache.find(key);
  if (found != cache.end()) return found->second;

  const InitialStateFoldGrid &grid =
      cached_initial_state_fold_grid(interaction, cfg, opt);
  std::vector<TVector3> momenta;
  momenta.reserve(opt.initial_state_fold_samples);
  for (int is = 0; is < opt.initial_state_fold_samples; ++is) {
    const std::uint64_t sample = static_cast<std::uint64_t>(is);
    const std::size_t ir =
        select_weighted_index(grid.radii, qel_lattice_u(sample, 0));
    const std::vector<WeightedValue> &pgrid = grid.momenta[ir];
    const double p = pgrid[select_weighted_index(
        pgrid, qel_lattice_u(sample, 1))].value;
    momenta.push_back(isotropic_vector(
        p, qel_lattice_u(sample, 2), qel_lattice_u(sample, 3)));
  }
  return cache.emplace(key, std::move(momenta)).first->second;
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
    genie::Interaction selected(base_interaction);
    set_probe_p4(selected.InitStatePtr(),get(values,Var::E));
    const auto *normalization = opt.fold_normalization == "event" ?
        &mec_state_normalization(model,selected,opt) : nullptr;
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
    genie::Interaction state(base_interaction);

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
      state.InitStatePtr()->TgtPtr()->SetHitNucP4(TLorentzVector(
          cluster_p3, std::sqrt(cluster_mass * cluster_mass +
                                cluster_p3.Mag2())));
      bool sample_below_threshold = false;
      bool sample_outside_phase_space = false;
      double density = folded_state_density(
          model, state, values, diff_vars, opt, &sample_below_threshold,
          &sample_outside_phase_space);
      // Match SelectEmpiricalKinematics' proposal support as well as its
      // physical phase space. The lower W edge is not the cluster mass.
      if(normalization && density>0 &&
         (state.Kine().W()<1.88 || state.Kine().W()>3. ||
          state.Kine().Q2()<std::max(.01,fold_em_q2_min(state)) ||
          state.Kine().Q2()>8.)) {
        density=0;
        sample_outside_phase_space=true;
      }
      if(normalization && density>0) {
        double er=state.InitState().ProbeE(genie::kRfHitNucRest);
        double norm=(er>=normalization->emin && er<=normalization->emax) ?
            normalization->sigma(er) : mec_state_shape_integral(model,selected,er,64);
        if(!(norm>0))throw std::runtime_error("Nonzero MEC density has zero state normalization");
        density *= normalization->total/norm;
      }
      sum += density;
      below_threshold += sample_below_threshold ? 1 : 0;
      outside_phase_space += sample_outside_phase_space ? 1 : 0;
      if (density > 0.0) ++accepted;
    }
  } else {
    const InitialStateFoldGrid &grid =
        cached_initial_state_fold_grid(base_interaction, cfg, grid_opt);
    const double removal_energy = grid.removal_energy;
    const std::vector<TVector3> &sample_momenta =
        cached_initial_state_fold_momenta(base_interaction, cfg, grid_opt);
    genie::Interaction state(base_interaction);

    for (int is = 0; is < opt.initial_state_fold_samples; ++is) {
      const TVector3 &p3 = sample_momenta[is];
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
