void write_metadata(std::ostream &out, int argc, char **argv,
                    const Options &opt)
{
  out << "# command";
  for (int i = 0; i < argc; ++i) out << "," << csv_escape(argv[i]);
  out << "\n";
  out << "# mode," << (opt.mode == Mode::Tune ? "tune" : "alg") << "\n";
  out << "# batch_config," << csv_escape(opt.batch_config_path) << "\n";
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
  out << "# folding_mode," << (opt.auto_fold ? "auto" : "manual") << "\n";
  out << "# jobs," << opt.jobs << "\n";
  out << "# qel_bin_fold," << (opt.qel_bin_fold ? "on" : "off") << "\n";
  out << "# qel_bin_fold_auto,"
      << (opt.auto_fold ? "model-dependent" : "off") << "\n";
  out << "# initial_state_fold,";
  if (opt.auto_fold) {
    out << "auto\n";
    out << "# initial_state_fold_samples,"
        << opt.initial_state_fold_samples << "\n";
    out << "# initial_state_fold_nr," << opt.initial_state_fold_nr << "\n";
    out << "# initial_state_fold_np," << opt.initial_state_fold_np << "\n";
    out << "# initial_state_fold_event_phase_space,on\n";
  } else if (opt.initial_state_fold.empty()) {
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
  if (opt.auto_fold || opt.qel_bin_fold) {
    out << "# qel_bin_width_energy," << opt.qel_bin_width_energy << "\n";
    out << "# qel_bin_width_costheta_l," << opt.qel_bin_width_costh << "\n";
    out << "# qel_fold_method," << opt.qel_fold_method << "\n";
    out << "# qel_reference_seed," << opt.qel_reference_seed << "\n";
    out << "# qel_rel_tol," << opt.qel_rel_tol << "\n";
    out << "# qel_abs_tol_genie," << opt.qel_abs_tol << "\n";
    out << "# qel_max_eval," << opt.qel_max_eval << "\n";
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

void write_header(std::ostream &out, const Options &opt)
{
  out << "row_id,";
  if (opt.batch_mode) out << "batch,target,";
  out << "component,observable,value,unit,status,message";
  for (Var v : output_vars()) out << "," << canonical_var_name(v);
  out << "\n";
}

void write_row(std::ostream &out, const Point &point,
               const ComponentValue &value, const Options &opt,
               double unit_factor)
{
  out << point.row_id << ",";
  if (opt.batch_mode) {
    out << csv_escape(point.batch) << "," << point.target << ",";
  }
  out << csv_escape(value.component) << ","
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

void write_point_result(std::ostream &out, const Point &point,
                        const std::vector<ComponentValue> &components,
                        const Options &opt, double unit_factor)
{
  ComponentValue total;
  total.component = "total";
  bool all_components_ok = !components.empty();
  bool have_kinematics = false;
  if (components.empty()) {
    total.status = "invalid";
    total.message = "No cross-section components were evaluated";
  }
  for (const ComponentValue &component : components) {
    if (component.status != "ok") all_components_ok = false;
    if (component.status != "ok" && total.status == "ok") {
      total.status = component.status;
      total.message =
          "Component " + component.component + ": " + component.message;
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

int run(int argc, char **argv)
{
  Invocation invocation = parse_invocation(argc, argv);
  Options opt = invocation.options;
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
  const std::vector<Point> points =
      opt.batch_mode ? make_batch_points(opt, invocation.batch)
                     : make_points(opt);
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

  std::map<std::pair<int, double>, std::unique_ptr<TuneContext>> tune_contexts;
  std::vector<const TuneContext *> point_tune_context(points.size(), nullptr);
  if (opt.mode == Mode::Tune) {
    for (std::size_t index = 0; index < points.size(); ++index) {
      const double E = get(points[index].values, Var::E);
      const std::pair<int, double> key = {points[index].target, E};
      auto found = tune_contexts.find(key);
      if (found == tune_contexts.end()) {
        Options context_opt = opt;
        context_opt.target = points[index].target;
        found = tune_contexts
                    .emplace(key, make_tune_context(E, context_opt))
                    .first;
        prime_tune_context_caches(*found->second, E);
      }
      point_tune_context[index] = found->second.get();
    }
  }

  std::filesystem::path out_path(opt.output_path);
  if (!out_path.parent_path().empty()) {
    std::filesystem::create_directories(out_path.parent_path());
  }
  std::ofstream out(opt.output_path);
  if (!out) throw std::runtime_error("Could not open output: " + opt.output_path);

  write_metadata(out, argc, argv, opt);
  write_header(out, opt);

  auto evaluate_and_write = [&](std::ostream &stream, std::size_t begin,
                                std::size_t end) {
    for (std::size_t index = begin; index < end; ++index) {
      const std::vector<ComponentValue> components =
          (opt.mode == Mode::Tune)
              ? evaluate_tune_point(points[index], *point_tune_context[index])
              : evaluate_alg_point(points[index], opt, *direct_model);
      write_point_result(stream, points[index], components, opt, unit_factor);
    }
  };

  const int njobs = std::min<int>(opt.jobs, static_cast<int>(points.size()));
  if (njobs <= 1) {
    evaluate_and_write(out, 0, points.size());
  } else {
    out.flush();
    const std::string worker_prefix =
        opt.output_path + ".workers-" + std::to_string(getpid()) + "-";
    std::vector<pid_t> pids;
    std::vector<std::filesystem::path> worker_paths(points.size());
    for (std::size_t index = 0; index < points.size(); ++index) {
      worker_paths[index] =
          worker_prefix + std::to_string(index) + ".csv.tmp";
    }
    for (int job = 0; job < njobs; ++job) {
      const pid_t pid = fork();
      if (pid < 0) throw std::runtime_error("fork failed for scan worker");
      if (pid == 0) {
        try {
          for (std::size_t index = static_cast<std::size_t>(job);
               index < points.size(); index += static_cast<std::size_t>(njobs)) {
            std::ofstream worker_out(worker_paths[index]);
            if (!worker_out) {
              throw std::runtime_error("could not open worker output");
            }
            evaluate_and_write(worker_out, index, index + 1);
            worker_out.close();
          }
          std::_Exit(0);
        } catch (const std::exception &e) {
          std::cerr << "xsec_scan worker " << job << ": " << e.what() << "\n";
          std::_Exit(1);
        }
      }
      pids.push_back(pid);
    }

    bool worker_failed = false;
    for (pid_t pid : pids) {
      int status = 0;
      if (waitpid(pid, &status, 0) < 0 || !WIFEXITED(status) ||
          WEXITSTATUS(status) != 0) {
        worker_failed = true;
      }
    }
    if (!worker_failed) {
      for (const std::filesystem::path &worker_path : worker_paths) {
        std::ifstream worker_in(worker_path);
        if (!worker_in) {
          worker_failed = true;
          break;
        }
        out << worker_in.rdbuf();
      }
    }
    for (const std::filesystem::path &worker_path : worker_paths) {
      std::error_code ec;
      std::filesystem::remove(worker_path, ec);
    }
    if (worker_failed) {
      throw std::runtime_error("one or more scan workers failed");
    }
  }

  out.close();
  std::cerr << "wrote " << opt.output_path << "\n";
  std::_Exit(0);
  return 0;
}
