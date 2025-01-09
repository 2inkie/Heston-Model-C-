#pragma once

#include <complex>
#include <string>

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

    std::complex<double> charFunction(double omega) const;

    // Fourier transform for the calculation of the option price
    double optionPrice() const;

    double calculateGreeks(const std::string &greek) const;

    void display() const;
};