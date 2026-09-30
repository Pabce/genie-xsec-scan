#include "../src/detail/qe_quadrature.hpp"
#include <iostream>
int main() {
  auto check = [](double a, double b, double tol) {
    if (std::abs(a - b) > tol)
      throw std::runtime_error("quadrature test failed");
  };
  long n = 100000;
  auto q = qe_numeric::integrate([](double x) { return std::exp(x); }, 0., 1.,
                                 1e-10, 1e-13, n);
  check(q.value, std::exp(1.) - 1, 1e-11);
  n = 100000;
  q = qe_numeric::integrate([](double x) { return std::sqrt(x); }, 0., 1., 1e-7,
                            1e-10, n);
  check(q.value, 2. / 3, 1e-8);
  check(qe_numeric::root([](double x) { return x * x - 2; }, 0., 2.),
        std::sqrt(2.), 1e-14);
  bool failed = false;
  try {
    n = 15;
    qe_numeric::integrate([](double x) { return std::exp(20 * x); }, 0., 1.,
                          1e-12, 1e-15, n);
  } catch (const std::exception &) {
    failed = true;
  }
  if (!failed)
    throw std::runtime_error("budget did not fail");
  std::cout << "quadrature tests passed\n";
}
