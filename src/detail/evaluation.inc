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
    if ((opt.auto_qel_fold || opt.qel_bin_fold) &&
        is_rosenbluth_qel(base_interaction, model)) {
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
    g_qe_accuracy_message.clear();
    out.value = evaluate_model(model, interaction, values, opt);
    if (!g_qe_accuracy_message.empty()) out.message = g_qe_accuracy_message;
    if (!std::isfinite(out.value)) throw std::runtime_error("non-finite result");
  } catch (const std::exception &e) {
    if (opt.strict) throw;
    out.message = e.what();
    out.status = out.message.rfind("QE nonconverged:", 0) == 0
                     ? "nonconverged" : "invalid";
    out.value = 0.0;
  }
  return out;
}

int event_generator_module_index(const genie::EventGeneratorI &generator,
                                 const std::string &name_fragment)
{
  const genie::Registry &config = generator.GetConfig();
  if (!config.Exists("NModules")) return -1;
  const int nmodules = config.GetInt("NModules");
  for (int index = 0; index < nmodules; ++index) {
    const std::string key = "Module-" + std::to_string(index);
    if (!config.Exists(key)) continue;
    const RgAlg module = config.GetAlg(key);
    if (module.name.find(name_fragment) != std::string::npos) return index;
  }
  return -1;
}

Options auto_fold_options(const Options &base,
                          const genie::EventGeneratorI &generator,
                          const genie::Interaction &interaction,
                          const genie::XSecAlgorithmI &model)
{
  Options resolved = base;
  resolved.auto_qel_fold = false;
  resolved.auto_initial_state_fold = false;
  if (!base.auto_fold) return resolved;

  const genie::ProcessInfo &process = interaction.ProcInfo();
  const int fermi = event_generator_module_index(generator, "FermiMover");
  if (process.IsQuasiElastic()) {
    const int kinematics =
        event_generator_module_index(generator, "QELKinematicsGenerator");
    resolved.auto_qel_fold = is_rosenbluth_qel(interaction, model) &&
                             fermi >= 0 && kinematics > fermi;
  } else if (process.IsResonant()) {
    const int kinematics =
        event_generator_module_index(generator, "RESKinematicsGenerator");
    resolved.auto_initial_state_fold = fermi >= 0 && kinematics > fermi;
  } else if (process.IsDeepInelastic()) {
    const int kinematics =
        event_generator_module_index(generator, "DISKinematicsGenerator");
    resolved.auto_initial_state_fold = fermi >= 0 && kinematics > fermi;
  } else if (process.IsMEC()) {
    const int mec_generator =
        event_generator_module_index(generator, "MECGenerator");
    resolved.auto_initial_state_fold =
        mec_generator >= 0 && is_empirical_mec(interaction, model);
  }
  static std::set<std::string> reported;
  const std::string process_name =
      genie::ScatteringType::AsString(process.ScatteringTypeId());
  const std::string report_key = process_name + "|" + model.Id().Key() +
                                 "|" + generator.Id().Key();
  if (reported.insert(report_key).second) {
    std::cerr << "auto-fold process=" << process_name
              << " model=" << model.Id().Key()
              << " generator=" << generator.Id().Key()
              << " qel=" << (resolved.auto_qel_fold ? "on" : "off")
              << " initial-state="
              << (resolved.auto_initial_state_fold ? "on" : "off")
              << " event-phase-space="
              << ((resolved.auto_qel_fold || resolved.auto_initial_state_fold)
                      ? "on" : "not-applicable")
              << "\n";
  }
  return resolved;
}

struct TuneComponentSpec {
  std::string name;
  std::unique_ptr<genie::Interaction> interaction;
  const genie::XSecAlgorithmI *model = nullptr;
  Options options;
  bool exact_zero = false;
  std::string equivalent_key;
};

struct TuneContext {
  std::unique_ptr<genie::GEVGDriver> driver;
  std::vector<TuneComponentSpec> components;
};

std::unique_ptr<TuneContext> make_tune_context(double E, const Options &opt)
{
  genie::InitialState init = make_initial_state(opt, E);
  auto context = std::make_unique<TuneContext>();
  context->driver = std::make_unique<genie::GEVGDriver>();
  context->driver->SetEventGeneratorList(opt.event_generator_list);
  context->driver->Configure(init);
  const genie::InteractionList *interactions = context->driver->Interactions();
  if (!interactions || interactions->empty()) {
    throw std::runtime_error("GEVGDriver returned no interactions");
  }
  for (const genie::Interaction *interaction : *interactions) {
    const genie::EventGeneratorI *generator =
        context->driver->FindGenerator(interaction);
    if (!generator || !generator->CrossSectionAlg()) continue;
    const genie::XSecAlgorithmI *model = generator->CrossSectionAlg();
    TuneComponentSpec spec;
    spec.name = component_name(*interaction, *model);
    spec.interaction = std::make_unique<genie::Interaction>(*interaction);
    spec.model = model;
    spec.options = auto_fold_options(opt, *generator, *interaction, *model);
    const genie::Target &target = interaction->InitState().Tgt();
    const genie::Resonance_t resonance = interaction->ExclTag().Resonance();
    bool uses_rs_em_proton_amplitudes = true;
    const genie::Registry &model_config = model->GetConfig();
    if (model_config.Exists("HelicityAmplEMpAlg")) {
      uses_rs_em_proton_amplitudes =
          model_config.GetAlg("HelicityAmplEMpAlg").name ==
          "genie::RSHelicityAmplModelEMp";
    }
    spec.exact_zero =
        interaction->ProcInfo().IsResonant() &&
        interaction->ProcInfo().IsEM() &&
        target.HitNucPdg() == genie::kPdgProton &&
        model->Id().Name() == "genie::BergerSehgalRESPXSec2014" &&
        uses_rs_em_proton_amplitudes &&
        (resonance == genie::kS11_1650 ||
         resonance == genie::kD13_1700 ||
         resonance == genie::kD15_1675 ||
         resonance == genie::kF17_1970);
    if (interaction->ProcInfo().IsDeepInelastic() &&
        interaction->ProcInfo().IsEM() && target.HitSeaQrk() &&
        model->Id().Name() == "genie::KNOTunedQPMDISPXSec") {
      // Electromagnetic DIS weights a sea quark and its antiquark identically.
      // GENIE nevertheless emits both interaction-list entries.
      spec.equivalent_key = model->Id().Key() + "|hit=" +
                            std::to_string(target.HitNucPdg()) + "|sea=" +
                            std::to_string(std::abs(target.HitQrkPdg()));
    }
    context->components.push_back(std::move(spec));
  }
  return context;
}

void prime_tune_context_caches(const TuneContext &context, double E)
{
  std::map<Var, double> energy_values;
  energy_values[Var::E] = E;
  for (const TuneComponentSpec &spec : context.components) {
    const genie::Interaction &interaction = *spec.interaction;
    if (spec.options.auto_initial_state_fold &&
        interaction.InitState().Tgt().IsNucleus() &&
        !interaction.ProcInfo().IsMEC()) {
      Options grid_opt = spec.options;
      grid_opt.qel_fold_nr = spec.options.initial_state_fold_nr;
      grid_opt.qel_fold_np = spec.options.initial_state_fold_np;
      const QELFoldConfig cfg = qel_fold_config(interaction);
      (void)cached_initial_state_fold_momenta(interaction, cfg, grid_opt);
    }
    if (shape_norm_applies_to_component(interaction, *spec.model,
                                        spec.options)) {
      (void)shape_norm_factor(*spec.model, interaction, energy_values,
                              spec.options);
    }
  }
}

std::vector<ComponentValue> evaluate_initial_state_fold_group(
    const Point &point, const std::vector<const TuneComponentSpec *> &specs)
{
  std::vector<ComponentValue> results(specs.size());
  for (std::size_t i = 0; i < specs.size(); ++i) {
    results[i].component = specs[i]->name;
  }
  if (specs.empty()) return results;

  const TuneComponentSpec &first = *specs.front();
  const Options &opt = first.options;
  try {
    if (!initial_state_fold_supports_diff(opt.diff_vars)) {
      throw std::runtime_error(
          "--initial-state-fold supports only d2 Eprime/omega,costheta_l");
    }
    const std::map<Var, double> kine_values =
        output_map_from_kine(solve_kinematics(point.values, *first.interaction));
    for (ComponentValue &result : results) result.kin_values = kine_values;

    Options grid_opt = opt;
    grid_opt.qel_fold_nr = opt.initial_state_fold_nr;
    grid_opt.qel_fold_np = opt.initial_state_fold_np;
    const QELFoldConfig cfg = qel_fold_config(*first.interaction);
    const InitialStateFoldGrid &grid =
        cached_initial_state_fold_grid(*first.interaction, cfg, grid_opt);
    const std::vector<TVector3> &sample_momenta =
        cached_initial_state_fold_momenta(*first.interaction, cfg, grid_opt);
    const double removal_energy = grid.removal_energy;

    std::vector<double> sums(specs.size(), 0.0);
    std::vector<long> accepted(specs.size(), 0);
    long below_threshold = 0;
    long outside_phase_space = 0;
    genie::Interaction state(*first.interaction);
    for (int is = 0; is < opt.initial_state_fold_samples; ++is) {
      if (!qel_bind_old_fermi_mover(
              &state, cfg, sample_momenta[is], removal_energy)) {
        continue;
      }
      bool sample_below_threshold = false;
      bool sample_outside_phase_space = false;
      const double jac = prepare_folded_state(
          state, point.values, opt.diff_vars, opt, &sample_below_threshold,
          &sample_outside_phase_space);
      below_threshold += sample_below_threshold ? 1 : 0;
      outside_phase_space += sample_outside_phase_space ? 1 : 0;
      if (jac <= 0.0) continue;

      for (std::size_t i = 0; i < specs.size(); ++i) {
        const genie::Interaction &source = *specs[i]->interaction;
        state.ExclTagPtr()->Copy(source.ExclTag());
        state.InitStatePtr()->TgtPtr()->SetHitQrkPdg(
            source.InitState().Tgt().HitQrkPdg());
        state.InitStatePtr()->TgtPtr()->SetHitSeaQrk(
            source.InitState().Tgt().HitSeaQrk());
        const genie::XSecAlgorithmI &model = *specs[i]->model;
        if (!model.ValidProcess(&state)) continue;
        const double value = model.XSec(&state, genie::kPSWQ2fE);
        if (!std::isfinite(value) || value <= 0.0) continue;
        sums[i] += value * jac;
        ++accepted[i];
      }
    }

    for (std::size_t i = 0; i < specs.size(); ++i) {
      const double folded =
          sums[i] / static_cast<double>(opt.initial_state_fold_samples);
      results[i].value = apply_shape_norm(
          *specs[i]->model, *specs[i]->interaction, point.values,
          specs[i]->options, folded);
      if (opt.initial_state_fold_debug) {
        std::cerr << "initial-state-fold component=" << specs[i]->name
                  << " samples=" << opt.initial_state_fold_samples
                  << " accepted=" << accepted[i]
                  << " below-threshold=" << below_threshold
                  << " outside-phase-space=" << outside_phase_space
                  << " density=" << folded << "\n";
      }
    }
  } catch (const std::exception &e) {
    if (opt.strict) throw;
    for (ComponentValue &result : results) {
      result.status = "invalid";
      result.value = 0.0;
      result.message = e.what();
    }
  }
  return results;
}

std::vector<ComponentValue> evaluate_tune_point(const Point &point,
                                                const TuneContext &context)
{
  std::vector<ComponentValue> components(context.components.size());
  std::vector<bool> ready(context.components.size(), false);
  std::vector<std::size_t> alias_owner(context.components.size(),
                                       context.components.size());
  std::map<std::string, std::size_t> equivalent_owner;
  std::map<std::string, std::vector<std::size_t>> fold_groups;

  for (std::size_t index = 0; index < context.components.size(); ++index) {
    const TuneComponentSpec &spec = context.components[index];
    const genie::Interaction &interaction = *spec.interaction;
    if (spec.exact_zero) {
      ComponentValue zero;
      zero.component = spec.name;
      try {
        if (spec.options.observable == Observable::Total &&
            spec.options.integrate_over.empty()) {
          if (has(point.values, Var::E)) {
            zero.kin_values[Var::E] = get(point.values, Var::E);
          }
        } else if (spec.options.integrate_over.empty()) {
          genie::Interaction tmp(interaction);
          zero.kin_values =
              output_map_from_kine(solve_kinematics(point.values, tmp));
        }
      } catch (const std::exception &e) {
        if (spec.options.strict) throw;
        zero.status = "invalid";
        zero.message = e.what();
      }
      components[index] = std::move(zero);
      ready[index] = true;
      continue;
    }
    if (!spec.equivalent_key.empty()) {
      auto owner = equivalent_owner.find(spec.equivalent_key);
      if (owner != equivalent_owner.end()) {
        alias_owner[index] = owner->second;
        continue;
      }
      equivalent_owner[spec.equivalent_key] = index;
    }
    if (initial_state_fold_applies(interaction, spec.options) &&
        !interaction.ProcInfo().IsMEC() &&
        interaction.InitState().Tgt().IsNucleus()) {
      std::ostringstream key;
      key << reinterpret_cast<std::uintptr_t>(spec.model)
          << "|proc=" << interaction.ProcInfo().ScatteringTypeId()
          << "|hit=" << interaction.InitState().Tgt().HitNucPdg()
          << "|samples=" << spec.options.initial_state_fold_samples
          << "|nr=" << spec.options.initial_state_fold_nr
          << "|np=" << spec.options.initial_state_fold_np
          << "|diff=" << diff_label(spec.options.diff_vars);
      fold_groups[key.str()].push_back(index);
    }
  }

  for (const auto &entry : fold_groups) {
    if (entry.second.size() < 2) continue;
    std::vector<const TuneComponentSpec *> specs;
    specs.reserve(entry.second.size());
    for (std::size_t index : entry.second) {
      specs.push_back(&context.components[index]);
    }
    std::vector<ComponentValue> grouped =
        evaluate_initial_state_fold_group(point, specs);
    for (std::size_t i = 0; i < entry.second.size(); ++i) {
      components[entry.second[i]] = std::move(grouped[i]);
      ready[entry.second[i]] = true;
    }
  }

  for (std::size_t index = 0; index < context.components.size(); ++index) {
    if (ready[index] || alias_owner[index] != context.components.size()) {
      continue;
    }
    const TuneComponentSpec &spec = context.components[index];
    components[index] = evaluate_component(
        spec.name, *spec.model, *spec.interaction, point.values, spec.options);
    ready[index] = true;
  }
  for (std::size_t index = 0; index < context.components.size(); ++index) {
    if (alias_owner[index] == context.components.size()) continue;
    components[index] = components[alias_owner[index]];
    components[index].component = context.components[index].name;
    ready[index] = true;
  }
  return components;
}

std::vector<ComponentValue> evaluate_alg_point(const Point &point,
                                               const Options &opt,
                                               const genie::XSecAlgorithmI &model)
{
  const double E = get(point.values, Var::E);
  Options point_opt = opt;
  point_opt.target = point.target;
  genie::InitialState init = make_initial_state(point_opt, E);
  std::unique_ptr<genie::InteractionList> list =
      make_algorithm_interactions(point_opt, init);

  std::set<std::string> seen;
  std::vector<ComponentValue> components;
  for (const genie::Interaction *interaction : *list) {
    if (point_opt.interaction_sum == "hit-state") {
      const std::string key = hit_state_key(*interaction);
      if (seen.count(key)) continue;
      seen.insert(key);
    }
    components.push_back(evaluate_component(component_name(*interaction, model),
                                            model, *interaction, point.values,
                                            point_opt));
  }
  return components;
}
