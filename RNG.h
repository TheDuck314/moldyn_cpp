#ifndef RNG_H
#define RNG_H

#include <random>

class RNG
{
private:
    static std::random_device rd;
    static std::mt19937 mt;
    static std::normal_distribution<double> normal;
    static std::uniform_real_distribution<double> rand01;

public:
    static double Gaussian(double sigma = 1.0) { return sigma * normal(mt); }
    static double UnitInterval() { return rand01(mt); }
};

std::random_device RNG::rd;
std::mt19937 RNG::mt(rd()); // mersenne twister
std::normal_distribution<double> RNG::normal;
std::uniform_real_distribution<double> RNG::rand01;


#endif