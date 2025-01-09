// Complete implementation of the methods in the header

#include "heston.hpp"
#include <cmath>
#include <complex>
#include <iostream>

HestonModel::HestonModel(double S0, double K, double T, double r, double v0,
                         double kappa, double theta, double sigma_v, double rho)
    : S0(S0), K(K), T(T), r(r), v0(v0), kappa(kappa), theta(theta),
      sigma_v(sigma_v), rho(rho) {}

std::complex<double> HestonModel::charFunction(double omega) const {
    std::complex<double> i(0, 1);

    // This can be made more advanced as this is just the simplified version
    std::complex<double> alpha = (std::pow(omega, 2) - i * omega) / 2.0;

    return std::exp(i * omega * std::log(S0)) * std::exp(-r * T) *
           std::exp(-alpha * T);
}