std::string g_qe_accuracy_message;
// Controlled LocalFGM QE integration. Radius is marginalized analytically in
// the sense of piecewise native momentum PDFs plus a scalar radial quadrature.
// This specialization deliberately rejects models whose binding depends on r.
struct QELRadialMeasure {
  std::vector<double> edges, mass, branches;
  std::vector<std::vector<double>> pdf, suffix;
  double dp = 0, rmax = 0, norm = 0;
  int bins = 0;
};

double qel_radial_mass(int A, double a, double b) {
  long budget = 10000;
  return qe_numeric::integrate(
             [&](double r) {
               return r * r * genie::utils::nuclear::Density(r, A);
             },
             a, b, 1e-10, 1e-15, budget)
      .value;
}

const QELRadialMeasure &
qel_native_radial_measure(const genie::Interaction &base,
                          const QELFoldConfig &cfg) {
  const auto &t = base.InitState().Tgt();
  const auto *lfg = dynamic_cast<const genie::LocalFGM *>(cfg.nucl_model);
  if (!lfg)
    throw std::runtime_error(
        "adaptive-theta requires LocalFGM");
  const auto &reg = lfg->GetConfig();
  auto flag = [&](const char *k) { return reg.Exists(k) && reg.GetBool(k); };
  if ((reg.Exists("SRC-Fraction") && reg.GetDouble("SRC-Fraction") != 0) ||
      flag("LFG-MomentumDependentErmv") || flag("LFG-ForcePositiveErmv") ||
      flag("LFG-UseMBDistForNegativeErmv") || cfg.fermi_momentum_dependent_ermv)
    throw std::runtime_error("adaptive-theta unsupported "
                             "correlated/r-dependent removal-energy model");
  std::ostringstream key;
  key << lfg << '|' << t.Pdg() << '|' << t.HitNucPdg() << '|'
      << std::setprecision(17) << cfg.nucl_r0;
  static std::map<std::string, QELRadialMeasure> cache;
  auto it = cache.find(key.str());
  if (it != cache.end())
    return it->second;
  QELRadialMeasure m;
  double pmax =
      reg.Exists("LFG-MomentumMax") ? reg.GetDouble("LFG-MomentumMax") : 1.;
  m.bins = static_cast<int>(1000 * pmax);
  if (m.bins < 2)
    throw std::runtime_error("LocalFGM histogram has fewer than two bins");
  m.dp = pmax / m.bins;
  m.rmax = 3 * cfg.nucl_r0 * std::cbrt(t.A());
  auto kf = [&](double r) {
    return lfg->LocalFermiMomentum(t, t.HitNucPdg(), r);
  };
  // GENIE's modified harmonic-oscillator densities can have a central dip.
  // Locate the single interior maximum, then retain both monotone branches.
  // Keep the original heavy-nucleus path (and its numerical results) unchanged.
  m.branches = {0};
  double peak = 0;
  for (int i = 1; i < 1024; ++i)
    if (kf(m.rmax * i / 1024) > kf(peak)) peak = m.rmax * i / 1024;
  if (peak > 0) {
    double lo = std::max(0., peak - m.rmax / 1024);
    double hi = std::min(m.rmax, peak + m.rmax / 1024);
    for (int i = 0; i < 100; ++i) {
      double a = lo + (hi-lo)/3, b = hi - (hi-lo)/3;
      if (kf(a) < kf(b)) lo = a; else hi = b;
    }
    m.branches.push_back((lo+hi)/2);
  }
  m.branches.push_back(m.rmax);
  for (size_t b = 1; b < m.branches.size(); ++b) {
    double lo = m.branches[b-1], hi = m.branches[b];
    bool rising = kf(hi) > kf(lo);
    double previous = kf(lo);
    for (int i = 1; i <= 1024; ++i) {
      double k = kf(lo + (hi-lo)*i/1024);
      if ((rising ? previous-k : k-previous) > 1e-12)
        throw std::runtime_error("Unsupported multi-extremum radial density");
      previous = k;
    }
  }
  m.edges = m.branches;
  // Native histogram fills at i*pmax/(nbins-1), not at bin centres.
  double step = pmax / (m.bins - 1);
  for (size_t b = 1; b < m.branches.size(); ++b) {
    double lo = m.branches[b-1], hi = m.branches[b];
    for (int i = 1; i < m.bins; ++i) {
      double k = i * step;
      if (k > std::min(kf(lo),kf(hi)) && k < std::max(kf(lo),kf(hi)))
        m.edges.push_back(qe_numeric::root(
            [&](double r) { return kf(r) - k; }, lo, hi));
    }
  }
  std::sort(m.edges.begin(), m.edges.end());
  m.norm = qel_radial_mass(t.A(), 0, m.rmax);
  for (size_t s = 0; s + 1 < m.edges.size(); ++s) {
    double r = (m.edges[s] + m.edges[s + 1]) / 2, K = kf(r);
    std::vector<double> pdf(m.bins, 0);
    double sum = 0;
    for (int i = 1; i < m.bins; ++i) {
      double p = i * step;
      if (p > K)
        break;
      int b = std::min(m.bins - 1, int(p / m.dp));
      pdf[b] += p * p;
      sum += p * p;
    }
    if (sum > 0)
      for (double &v : pdf)
        v /= sum * m.dp;
    double mass = qel_radial_mass(t.A(), m.edges[s], m.edges[s + 1]) / m.norm;
    if (sum == 0) {
      std::cerr << "qe-native-empty-tail probability=" << mass << "\n";
      if (mass > 1e-5)
        throw std::runtime_error(
            "Native LocalFGM empty-distribution radial mass exceeds 1e-5");
    }
    // Check cloned piecewise measure against the active GENIE implementation.
    if (sum > 0)
      for (int b : {1, int(.5 * K / m.dp), int(.95 * K / m.dp)}) {
        b = std::max(0, std::min(b, m.bins - 1));
        double native = lfg->Prob((b + .5) * m.dp, -1, t, r);
        if (std::abs(native - pdf[b] * m.dp) > 1e-11)
          throw std::runtime_error(
              "Native LocalFGM histogram adapter mismatch");
      }
    m.pdf.push_back(std::move(pdf));
    m.mass.push_back(mass);
  }
  m.suffix.assign(m.pdf.size() + 1, std::vector<double>(m.bins, 0));
  for (int s = int(m.pdf.size()) - 1; s >= 0; --s)
    for (int j = 0; j < m.bins; ++j)
      m.suffix[s][j] = m.suffix[s + 1][j] + m.mass[s] * m.pdf[s][j];
  std::cerr << "qe-adaptive-cache target=" << t.Pdg()
            << " hit=" << t.HitNucPdg() << " radial_shells=" << m.pdf.size()
            << " momentum_bins=" << m.bins << "\n";
  return cache.emplace(key.str(), std::move(m)).first->second;
}

double rosenbluth_qel_adaptive_theta(const genie::XSecAlgorithmI &model,
                                     const genie::Interaction &original,
                                     const std::map<Var, double> &values,
                                     const Options &opt) {
  if (opt.qel_fold_pmax > 0 || opt.qel_fold_kf > 0 ||
      opt.qel_fold_removal_energy >= 0)
    throw std::runtime_error(
        "adaptive-theta requires unmodified native nuclear-state parameters");
  genie::Interaction base(original), center_interaction(original);
  auto center = solve_kinematics(values, center_interaction);
  set_probe_p4(base.InitStatePtr(), center.E);
  auto cfg = qel_fold_config(base);
  const auto &measure = qel_native_radial_measure(base, cfg);
  const auto *lfg = dynamic_cast<const genie::LocalFGM *>(cfg.nucl_model);
  const auto &target = base.InitState().Tgt();
  double mi = particle_mass(base.InitState().ProbePdg()),
         ml = final_lepton_mass(base);
  double pin = std::sqrt(std::max(0., center.E * center.E - mi * mi)),
         pout =
             std::sqrt(std::max(0., center.Eprime * center.Eprime - ml * ml));
  if (pin == 0 || pout == 0)
    return 0;
  TVector3 kin(0, 0, pin),
      klep(pout * std::sqrt(
                      std::max(0., 1 - center.costheta_l * center.costheta_l)),
           0, pout * center.costheta_l);
  TVector3 q = kin - klep;
  double qm = q.Mag();
  if (qm == 0)
    return 0;
  auto qhat = q.Unit();
  TVector3 xhat = kin - kin.Dot(qhat) * qhat;
  if (xhat.Mag() < 1e-14)
    xhat = unit_perpendicular(qhat);
  else
    xhat = xhat.Unit();
  auto yhat = qhat.Cross(xhat);
  double recoil = particle_mass(base.RecoilNucleonPdg()),
         removal = qel_fold_removal_energy(base, opt),
         count = qel_hit_nucleon_count(target);
  long budget = opt.qel_max_eval, calls = 0;
  double outer_error = 0, inner_error = 0;
  auto energy = [&](double p) {
    genie::Interaction state(base);
    if (!qel_bind_old_fermi_mover(&state, cfg, p * qhat, removal))
      return -1.;
    return state.InitState().Tgt().HitNucP4().E();
  };
  auto cosine = [&](double p) {
    double ef = energy(p) + center.omega;
    return (ef * ef - recoil * recoil - p * p - qm * qm) / (2 * p * qm);
  };
  auto recoil_kf = [&](double r) {
    return lfg->LocalFermiMomentum(target, base.RecoilNucleonPdg(), r);
  };
  auto radial_tail = [&](int bin, double r) {
    if (r >= measure.rmax) return 0.;
    size_t s = std::min(
        measure.pdf.size() - 1,
        size_t(std::upper_bound(measure.edges.begin(), measure.edges.end(), r) -
               measure.edges.begin() - 1));
    return measure.suffix[s + 1][bin] +
           measure.pdf[s][bin] *
               qel_radial_mass(target.A(), r, measure.edges[s + 1]) /
               measure.norm;
  };
  auto radial_pdf = [&](int bin, double pf) {
    double result = 0;
    for (size_t b = 1; b < measure.branches.size(); ++b) {
      double lo = measure.branches[b-1], hi = measure.branches[b];
      double kl = recoil_kf(lo), kh = recoil_kf(hi);
      if (pf <= std::min(kl,kh)) continue;
      if (pf < std::max(kl,kh)) {
        double r = qe_numeric::root(
            [&](double x) { return recoil_kf(x)-pf; }, lo, hi);
        if (kl < kh) hi = r; else lo = r;
      }
      result += radial_tail(bin,lo)-radial_tail(bin,hi);
    }
    return std::max(0., result);
  };
  auto kernel = [&](double p, int bin) {
    double en = energy(p), ef = en + center.omega, c = cosine(p);
    if (ef <= recoil || !std::isfinite(c) || std::abs(c) >= 1)
      return 0.;
    double pf = std::sqrt(ef * ef - recoil * recoil),
           density = radial_pdf(bin, pf);
    if (density <= 0)
      return 0.;
    double st = std::sqrt(1 - c * c);
    auto eval = [&](double phi, int which) {
      if (++calls > opt.qel_max_eval)
        throw std::runtime_error(
            "QE nonconverged: model evaluation budget exhausted");
      if (calls % 100000 == 0)
        std::cerr << "qe-adaptive-progress omega=" << center.omega
                  << " calls=" << calls << "\n";
      TVector3 p3 =
          p * (c * qhat + st * (std::cos(phi) * xhat + std::sin(phi) * yhat));
      genie::Interaction fixed(base);
      fixed.SetBit(genie::kISkipProcessChk);
      fixed.SetBit(genie::kISkipKinematicChk);
      fixed.SetBit(genie::kIAssumeFreeNucleon);
      qel_bind_old_fermi_mover(&fixed, cfg, p3, removal);
      const auto lim = fixed.PhaseSpace().Limits(genie::kKVQ2);
      if (which == 1)
        return center.Q2 - lim.min;
      if (which == 2)
        return lim.max - center.Q2;
      if (center.Q2 < lim.min || center.Q2 > lim.max)
        return 0.;
      fixed.KinePtr()->SetQ2(center.Q2);
      double v = model.XSec(&fixed, genie::kPSQ2fE);
      if (!std::isfinite(v) || v < 0)
        throw std::runtime_error("QE invalid model cross section");
      const auto &hit = fixed.InitState().Tgt().HitNucP4();
      double pk = hit.E() * center.E - hit.Vect().Dot(kin);
      // Invariant two-body measure: dQ2 dphi = 4*sqrt((p.k)^2-p^2*k^2)
      // *delta((p+q)^2-Mf^2)*d^3k'/E'. Averaging initial cos(p,q)
      // contributes 1/(4*p*|q|); final lab azimuth cancels 1/(2*pi).
      double flux = std::sqrt(std::max(0., pk * pk - hit.M2() * mi * mi));
      return v * flux * pout / (p * qm);
    };
    std::vector<double> edges = {0, kPi};
    // In this basis p.k is monotonic on [0,pi]. Locate phase-limit crossings.
    for (int w = 1; w <= 2; ++w)
      for (int i = 0; i < 8; ++i) {
        double a = kPi * i / 8, b = kPi * (i + 1) / 8, fa = eval(a, w),
               fb = eval(b, w);
        if (fa * fb < 0)
          edges.push_back(
              qe_numeric::root([&](double x) { return eval(x, w); }, a, b));
      }
    std::sort(edges.begin(), edges.end());
    double value = 0, angular_error = 0;
    for (size_t i = 1; i < edges.size(); ++i) {
      auto v = qe_numeric::integrate(
          [&](double phi) { return eval(phi, 0); }, edges[i - 1], edges[i],
          opt.qel_rel_tol * .1, opt.qel_abs_tol * .01, budget);
      value += v.value;
      angular_error += v.error;
    }
    inner_error = std::max(inner_error, angular_error);
    return density * count * value / kPi;
  };
  // Every change in the active radial-shell PDF becomes an explicit p break.
  std::vector<double> pauli_breaks;
  double pmax = measure.bins * measure.dp;
  for (double r : measure.edges) {
    double k = recoil_kf(r);
    auto f = [&](double p) {
      return energy(p) + center.omega - std::sqrt(recoil * recoil + k * k);
    };
    if (f(0) * f(pmax) < 0)
      pauli_breaks.push_back(qe_numeric::root(f, 0, pmax));
  }
  // The lower triangle boundary can enclose a narrow allowed island. Split
  // at its maximum before searching for its two roots, even if both roots
  // lie inside a single native histogram bin.
  auto lower_derivative = [&](double p) {
    double h = 1e-5, lo = std::max(0., p - h), hi = p + h;
    double de = (energy(hi) - energy(lo)) / (hi - lo);
    return 2 * (energy(p) + center.omega) * de - 2 * (p - qm);
  };
  if (lower_derivative(0) * lower_derivative(pmax) < 0)
    pauli_breaks.push_back(qe_numeric::root(lower_derivative, 0, pmax));
  std::sort(pauli_breaks.begin(), pauli_breaks.end());
  double total = 0;
  for (int j = 0; j < measure.bins; ++j) {
    if (measure.suffix[0][j] == 0)
      continue;
    double a = std::max(1e-12, j * measure.dp), b = (j + 1) * measure.dp;
    std::vector<double> edges = {a, b};
    for (auto it =
             std::upper_bound(pauli_breaks.begin(), pauli_breaks.end(), a);
         it != pauli_breaks.end() && *it < b; ++it)
      edges.push_back(*it);
    std::sort(edges.begin(), edges.end());
    auto partitions = edges;
    for (size_t k = 1; k < partitions.size(); ++k) {
      double lo = partitions[k - 1], hi = partitions[k];
      for (int sign : {-1, 1}) {
        auto f = [&](double p) {
          double ef = energy(p) + center.omega;
          return ef * ef - recoil * recoil - (p + sign * qm) * (p + sign * qm);
        };
        if (f(lo) * f(hi) < 0)
          edges.push_back(qe_numeric::root(f, lo, hi));
      }
      auto threshold = [&](double p) {
        return energy(p) + center.omega - recoil;
      };
      if (threshold(lo) * threshold(hi) < 0)
        edges.push_back(qe_numeric::root(threshold, lo, hi));
    }
    std::sort(edges.begin(), edges.end());
    for (size_t i = 1; i < edges.size(); ++i) {
      auto v = qe_numeric::integrate(
          [&](double p) { return kernel(p, j); }, edges[i - 1], edges[i],
          opt.qel_rel_tol * .3, opt.qel_abs_tol / measure.bins, budget);
      total += v.value;
      outer_error += v.error;
    }
  }
  // The marginalized momentum PDF integrates to at most one. Therefore the
  // maximum angular absolute error times count/pi bounds its p integral.
  double estimated_error = outer_error + count * inner_error / kPi;
  if (estimated_error >
      std::max(opt.qel_abs_tol, opt.qel_rel_tol * std::abs(total)))
    throw std::runtime_error(
        "QE nonconverged: combined quadrature error exceeds tolerance");
  std::ostringstream diagnostic;
  diagnostic << std::setprecision(10)
             << "adaptive-theta; estimated_abs_error_genie=" << estimated_error
             << "; outer_abs_error_genie=" << outer_error
             << "; max_angular_abs_error_genie=" << inner_error
             << "; calls=" << calls << "; rel_tol=" << opt.qel_rel_tol;
  g_qe_accuracy_message = diagnostic.str();
  std::cerr << "qe-adaptive-result target=" << target.Pdg()
            << " hit=" << target.HitNucPdg() << " omega=" << center.omega
            << " value=" << total << " outer_error=" << outer_error
            << " max_angular_error=" << inner_error << " calls=" << calls
            << " status=converged\n";
  return total;
}

// Validation backend: continuous native histogram sampling and old generator
// Q2 -> lab lepton mapping. No exact-theta Jacobian or radial Pauli marginal.
double rosenbluth_qel_native_reference(const genie::XSecAlgorithmI &model,
                                       const genie::Interaction &original,
                                       const std::map<Var, double> &values,
                                       const Options &opt) {
  genie::Interaction base(original), tmp(original);
  auto center = solve_kinematics(values, tmp);
  set_probe_p4(base.InitStatePtr(), center.E);
  auto cfg = qel_fold_config(base);
  const auto &m = qel_native_radial_measure(base, cfg);
  const auto &t = base.InitState().Tgt();
  const int nr = 8192;
  std::vector<double> rcdf(nr + 1, 0);
  for (int i = 1; i <= nr; ++i)
    rcdf[i] =
        rcdf[i - 1] +
        qel_radial_mass(t.A(), m.rmax * (i - 1) / nr, m.rmax * i / nr) / m.norm;
  for (double &v : rcdf)
    v /= rcdf.back();
  std::vector<std::vector<double>> pcdf(m.pdf.size(),
                                        std::vector<double>(m.bins + 1, 0));
  for (size_t s = 0; s < m.pdf.size(); ++s)
    for (int j = 0; j < m.bins; ++j)
      pcdf[s][j + 1] = pcdf[s][j] + m.pdf[s][j] * m.dp;
  auto invert = [](const std::vector<double> &cdf, double u) {
    size_t j = std::min(
        cdf.size() - 2,
        size_t(std::upper_bound(cdf.begin(), cdf.end(), u) - cdf.begin() - 1));
    return j + (u - cdf[j]) / (cdf[j + 1] - cdf[j]);
  };
  double shift[6];
  std::uint64_t state = static_cast<std::uint64_t>(opt.qel_reference_seed);
  for (double &v : shift) {
    state += 0x9e3779b97f4a7c15ULL;
    auto z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z ^= z >> 31;
    v = (z >> 11) * 0x1.0p-53;
  }
  auto u = [&](int i, int d) {
    double x = qel_lattice_u(i, d) + shift[d];
    return x - std::floor(x);
  };
  double sum = 0, removal = qel_fold_removal_energy(base, opt);
  long hits = 0;
  for (int i = 0; i < opt.qel_fold_samples; ++i) {
    double r = invert(rcdf, u(i, 0)) * m.rmax / nr;
    size_t shell =
        std::min(m.pdf.size() - 1,
                 size_t(std::upper_bound(m.edges.begin(), m.edges.end(), r) -
                        m.edges.begin() - 1));
    if (pcdf[shell].back() == 0)
      continue;
    double p = invert(pcdf[shell], u(i, 1)) * m.dp, c = 2 * u(i, 2) - 1,
           s = std::sqrt(1 - c * c), phi = 2 * kPi * u(i, 3);
    TVector3 p3(p * s * std::cos(phi), p * s * std::sin(phi), p * c);
    genie::Interaction fixed(base);
    fixed.SetBit(genie::kISkipProcessChk);
    fixed.SetBit(genie::kISkipKinematicChk);
    fixed.SetBit(genie::kIAssumeFreeNucleon);
    fixed.InitStatePtr()->TgtPtr()->SetHitNucPosition(r);
    if (!qel_bind_old_fermi_mover(&fixed, cfg, p3, removal))
      continue;
    auto lim = fixed.PhaseSpace().Limits(genie::kKVQ2);
    if (lim.max <= lim.min)
      continue;
    double q2 = lim.min + (lim.max - lim.min) * u(i, 4);
    double v = rosenbluth_qel_old_generator_q2_xsec(model, &fixed, cfg, q2,
                                                    2 * kPi * u(i, 5));
    if (v > 0 && qel_lab_lepton_in_bin(fixed, center, opt)) {
      sum += v * (lim.max - lim.min);
      ++hits;
    }
    if ((i + 1) % 1000000 == 0)
      std::cerr << "qe-reference-progress samples=" << i + 1 << "/"
                << opt.qel_fold_samples << " hits=" << hits << "\n";
  }
  std::ostringstream msg;
  msg << "native-q2-reference; seed=" << opt.qel_reference_seed
      << "; samples=" << opt.qel_fold_samples << "; hits=" << hits;
  g_qe_accuracy_message = msg.str();
  return sum / (opt.qel_fold_samples * opt.qel_bin_width_energy *
                opt.qel_bin_width_costh);
}
