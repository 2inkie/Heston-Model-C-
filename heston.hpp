#pragma once

#include <complex>
#include <string>
#include <vector>

class HestonModel {
  public:
    double S0;      // stock price
    double K;       // strike
    double T;       // maturity
    double r;       // risk free rate
    double v0;      // variance
    double kappa;   // mean reversion
    double theta;   // variance in long-term
    double sigma_v; // the volatility of the volatility
    double rho;     // asset and colatility correlation

    HestonModel(double S0, double K, double T, double r, double v0,
                double kappa, double theta, double sigma_v, double rho);

    double monteCarloPrice(int paths, int steps) const;

  private:
    void generateCorNorms(std::vector<double> &z1, std::vector<double> &z2,
                          int size) const;
};