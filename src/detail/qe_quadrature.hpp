#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace qe_numeric {
struct Result {
  double value = 0, error = 0;
  long evaluations = 0;
};
// Embedded Gauss 7 / Kronrod 15. Caller supplies known support boundaries.
template <class F> Result rule(F &&f, double a, double b) {
  static constexpr double x[] = {.9914553711208126, .9491079123427585,
                                 .8648644233597691, .7415311855993945,
                                 .5860872354676911, .4058451513773972,
                                 .2077849550078985, 0};
  static constexpr double w[] = {.02293532201052922, .06309209262997855,
                                 .1047900103222502,  .1406532597155259,
                                 .1690047266392679,  .1903505780647854,
                                 .2044329400752989,  .2094821410847278};
  static constexpr double g[] = {.1294849661688697, .2797053914892767,
                                 .3818300505051189, .4179591836734694};
  double mid = (a + b) / 2, h = (b - a) / 2, c = f(mid), k = w[7] * c,
         v = g[3] * c, ab = w[7] * std::abs(c);
  for (int i = 0; i < 7; ++i) {
    double l = f(mid - h * x[i]), r = f(mid + h * x[i]);
    k += w[i] * (l + r);
    ab += w[i] * (std::abs(l) + std::abs(r));
    if (i % 2)
      v += g[i / 2] * (l + r);
  }
  if (!std::isfinite(k) || !std::isfinite(v))
    throw std::runtime_error("QE nonfinite quadrature integrand");
  return {h * k,
          std::max(std::abs(h * (k - v)), 50 * 2.220446049250313e-16 * h * ab),
          15};
}
template <class F>
Result integrate(F &&f, double a, double b, double rel, double abs,
                 long &budget, int depth = 0) {
  if (b <= a)
    return {};
  if (budget < 15)
    throw std::runtime_error("QE nonconverged: evaluation budget exhausted");
  budget -= 15;
  auto q = rule(f, a, b);
  if (q.error <= std::max(abs, rel * std::abs(q.value)))
    return q;
  if (depth >= 48)
    throw std::runtime_error("QE nonconverged: subdivision depth exhausted");
  auto l = integrate(f, a, (a + b) / 2, rel, abs / 2, budget, depth + 1),
       r = integrate(f, (a + b) / 2, b, rel, abs / 2, budget, depth + 1);
  return {l.value + r.value, l.error + r.error,
          q.evaluations + l.evaluations + r.evaluations};
}
template <class F> double root(F &&f, double a, double b) {
  double fa = f(a), fb = f(b);
  if (fa == 0)
    return a;
  if (fb == 0)
    return b;
  if (fa * fb > 0)
    throw std::runtime_error("QE boundary root not bracketed");
  for (int i = 0; i < 60; ++i) {
    double m = (a + b) / 2, fm = f(m);
    if (fa * fm <= 0) {
      b = m;
      fb = fm;
    } else {
      a = m;
      fa = fm;
    }
  }
  return (a + b) / 2;
}
} // namespace qe_numeric
