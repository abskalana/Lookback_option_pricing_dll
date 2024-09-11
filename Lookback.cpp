#include "pch.h"
#include "Lookback.h"
#include <iostream>
#include <cmath> 
#include <vector>
#include <numeric>
#include <random>
#include <algorithm>


std::vector<std::vector<double>> generatePath(double S, double K, double T, double r, double sigma, int numSimulations, double nStep) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0, 1);

    double dt = T / nStep;
    std::vector<std::vector<double>> stockPath(numSimulations, std::vector<double>(static_cast<int>(T / dt) + 1, 0.0));

    for (int i = 0; i < numSimulations; ++i) {
        stockPath[i][0] = S;

        for (int t = 1; t < stockPath[i].size(); ++t) {
            double z = dist(gen);
            stockPath[i][t] = stockPath[i][t - 1] * exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * z);
        }
    }
    return stockPath;
}

double priceLookbackoption(double S, double K, double T, double r, double sigma, int numSimulations = 1000, int nStep = 100, int optionType = 0, int lookbackStrikeType = 0) {
    std::vector<std::vector<double>> simulations = generatePath(S, K, T, r, sigma, numSimulations, nStep);


    std::vector<double> payoffs(numSimulations, 0.0);


    for (int i = 0; i < numSimulations; ++i) {
        double lastPrice = simulations[i][simulations[i].size() - 1];
        double maxPrice = *std::max_element(simulations[i].begin(), simulations[i].end());
        double minPrice = *std::min_element(simulations[i].begin(), simulations[i].end());

        if (lookbackStrikeType == 1) { // floating

            if (optionType == 1) {

                payoffs[i] = lastPrice > minPrice ? lastPrice - minPrice : 0.0;
            }
            else {

                payoffs[i] = lastPrice > maxPrice ? lastPrice - maxPrice : 0.0;
            }

        }
        else { //Fixed
            if (optionType == 1) {

                payoffs[i] = K > minPrice ? K - minPrice : 0.0;
            }
            else {

                payoffs[i] = maxPrice > K ? maxPrice - K : 0.0;
            }
        }

    }

    double discountFactor = exp(-r * T);
    double lookbackOptionPrice = discountFactor * std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / numSimulations;

    return lookbackOptionPrice;
}

double computeLookbackoption(double S, double K, double T, double r, double sigma, int numSimulations = 1000, int nStep = 100, int optionType = 0, int lookbackStrikeType = 0, int indicatorType = 0, double bumps = 0.01) {

    if (indicatorType == 1) {
        double optionPriceHigh = priceLookbackoption(S * (1 + bumps), K, T, r, sigma, numSimulations, nStep, optionType, lookbackStrikeType);
        double optionPriceLow = priceLookbackoption(S * (1 - bumps), K, T, r, sigma, numSimulations, nStep, optionType, lookbackStrikeType);
        return (optionPriceHigh - optionPriceLow) / (2 * bumps * S);
    }

    if (indicatorType == 2) { //theta
        double optionPriceHigh = priceLookbackoption(S, K, T * (1 + bumps), r, sigma, numSimulations, nStep, optionType, lookbackStrikeType);
        double optionPriceLow = priceLookbackoption(S, K, T * (1 - bumps), r, sigma, numSimulations, nStep, optionType, lookbackStrikeType);
        return (optionPriceHigh - optionPriceLow) / (2 * bumps * T);
    }

    if (indicatorType == 3) { //rho
        double optionPriceHigh = priceLookbackoption(S, K, T, r * (1 + bumps), sigma, numSimulations, nStep, optionType, lookbackStrikeType);
        double optionPriceLow = priceLookbackoption(S, K, T, r * (1 - bumps), sigma, numSimulations, nStep, optionType, lookbackStrikeType);
        return (optionPriceHigh - optionPriceLow) / (2 * bumps * r);
    }
    if (indicatorType == 4) { //vega
        double optionPriceHigh = priceLookbackoption(S, K, T, r, sigma * (1 + bumps), numSimulations, nStep, optionType, lookbackStrikeType);
        double optionPriceLow = priceLookbackoption(S, K, T, r, sigma * (1 - bumps), numSimulations, nStep, optionType, lookbackStrikeType);
        return (optionPriceHigh - optionPriceLow) / (2 * bumps * sigma);
    }

    if (indicatorType == 11) { //gamma
        double deltaHigh = computeLookbackoption(S * (1 + bumps), K, T, r, sigma, numSimulations, nStep, optionType, lookbackStrikeType, 1, bumps);
        double deltaLow = computeLookbackoption(S * (1 - bumps), K, T, r, sigma, numSimulations, nStep, optionType, lookbackStrikeType, 1, bumps);
        return (deltaHigh - deltaLow) / (2 * bumps * S);
    }

    if (indicatorType == 14) { //vanna
        double deltaHigh = computeLookbackoption(S, K, T, r, sigma * (1 + bumps), numSimulations, nStep, optionType, lookbackStrikeType, 1, bumps);
        double deltaLow = computeLookbackoption(S, K, T, r, sigma * (1 - bumps), numSimulations, nStep, optionType, lookbackStrikeType, 1, bumps);
        return (deltaHigh - deltaLow) / (2 * bumps * sigma);
    }
    return priceLookbackoption(S, K, T, r, sigma, numSimulations, nStep, optionType, lookbackStrikeType);
}

int main() {

    double S = 100;
    double K = 100;
    double T = 1;
    double r = 0.1;
    double sigma = 0.2;
    int Nsim = 100;
    int Nstep = 100;
    int optionType = 0;
    int strikeType = 0;
    int indicator = 1;
    double bumps = 1.0;

    double result = computeLookbackoption(S, K, T, r, sigma, Nsim, Nstep, optionType, strikeType, indicator, bumps);
    std::cout<< result;

}