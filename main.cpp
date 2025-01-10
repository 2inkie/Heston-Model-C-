#include "heston.hpp"
#include <iostream>

int main() {
    double S0, K, T, r, v0, kappa, theta, sigma_v, rho;
    int paths = 10000;
    int steps = 252;

    std::cin >> S0 >> K >> T >> r >> v0 >> kappa >> theta >> sigma_v >> rho;

    HestonModel model(S0, K, T, r, v0, kappa, theta, sigma_v, rho);

    double price = model.monteCarloPrice(paths, steps);
    std::cout << price << std::endl;

    double noClose;
    std::cin >> noClose;
}