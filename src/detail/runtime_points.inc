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

std::vector<Point> make_batch_points(const Options &base,
                                     const BatchDocument &document)
{
  std::vector<Point> all_points;
  for (const BatchScan &batch_scan : document.scans) {
    Options scan_opt = base;
    scan_opt.scans.clear();
    scan_opt.points_path.clear();
    for (const std::string &fixed_text : batch_scan.fixed) {
      const auto fixed = parse_fixed(fixed_text);
      scan_opt.fixed[fixed.first] = fixed.second;
    }
    std::set<Var> scan_vars;
    for (const std::string &range_text : batch_scan.ranges) {
      const RangeSpec range = parse_range(range_text);
      if (!scan_vars.insert(range.var).second) {
        throw std::runtime_error("Batch scan " + batch_scan.name +
                                 " repeats range variable " +
                                 canonical_var_name(range.var));
      }
      if (scan_opt.fixed.count(range.var)) {
        throw std::runtime_error("Batch scan " + batch_scan.name +
                                 " fixes and ranges " +
                                 canonical_var_name(range.var));
      }
      scan_opt.scans.push_back(range);
    }
    std::vector<Point> points = make_points(scan_opt);
    for (Point &point : points) {
      point.batch = batch_scan.name;
      point.row_id = static_cast<int>(all_points.size());
      all_points.push_back(std::move(point));
    }
  }
  return all_points;
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
