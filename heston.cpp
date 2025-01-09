// Complete implementation of the methods in the header

#include "heston.hpp"
#include <cmath>
#include <complex>
#include <iostream>
#include <random>

HestonModel::HestonModel(double S0, double K, double T, double r, double v0,
                         double kappa, double theta, double sigma_v, double rho)
    : S0(S0), K(K), T(T), r(r), v0(v0), kappa(kappa), theta(theta),
      sigma_v(sigma_v), rho(rho) {}

// Cholesky decomposition
void HestonModel::generateCorNorms(std::vector<double> &z1,
                                   std::vector<double> &z2, int size) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, 1.0);

    z1.resize(size);
    z2.resize(size);

    for (int i = 0; i < size; i++) {
        double z = dist(gen);
        z1[i] = dist(gen);
        z2[i] = rho * z + std::sqrt(1 - rho * rho) * dist(gen);
    }
}

double HestonModel::monteCarloPrice(int paths, int steps) const {
    double dt = T / steps;
    double discountFactor = std::exp(-r * T);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, 1.0);

    double sumPayoffs = 0;

    // This can be optimized
    for (int path = 0; path < paths; ++path) {
        double S = S0;
        double v = v0;

        for (int step = 0; step < steps; ++step) {
            double z1 = dist(gen);
            double z2 = rho * z1 + std::sqrt(1 - rho * rho) * dist(gen);

            v = std::max(0.0, v + kappa * (theta - v) * dt +
                                  sigma_v * std::sqrt(v * dt) * z2);

            S *= std::exp((r - 0.5 * v) * dt + std::sqrt(v * dt) * z1);
        }

        sumPayoffs += std::max(S - K, 0.0);
    }
    return discountFactor * (sumPayoffs / paths);
}