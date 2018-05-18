#include "phasespace.hpp"
#include "BG.hpp"
#include <fstream>
#include <iostream>
#include <random>

using namespace std;

void test_BG_CO()
{
    std::mt19937 eng(std::random_device{}());
    std::uniform_real_distribution<double> urd(0, 1);
    const int n = 101;
    //Part 1: generate a new physical kinematics
    std::vector<L4v<double>> p(n);
    //z+ corresponding to 0, z- corresponding to inf
    //we use the convention that outgoing are positive
    double scale = 1e12;
    p[n - 1] = { -scale, 0, 0, -scale };
    p[0] = { -scale, 0, 0, 1 };
    L4v<double> Q{ 2 * scale, 0, 0, 0 };
    {
        double r[3 * (n - 2) - 4];
        for (int i = 0; i < 3 * (n - 2) - 4; ++i)
        {
            r[i] = urd(eng);
        }
        array<double, n> vm;
        for (int i = 0; i < n; ++i)
        {
            vm[i] = 0;
        }
        double dps;
        athl::rcm(n - 2, Q, r, &vm[0], &p[1], dps);
    }
    std::vector<int> vh(n);
    for (size_t i = 0; i < n; ++i)
    {
        vh[i] = (i % 2) * 2 - 1;
    }
    std::complex<double> res = amp_CO_BG(p, vh);
    std::cout << res << std::endl;
}
void test_BG_CO_arra()
{
    athl::arra_default_precision = 128;
    std::mt19937 eng(std::random_device{}());
    std::uniform_real_distribution<double> urd(0, 1);
    const int n = 8;
    //Part 1: generate a new physical kinematics
    std::vector<L4v<arra>> p(n);
    //z+ corresponding to 0, z- corresponding to inf
    //we use the convention that outgoing are positive
    arra scale = 1e12;
    p[n - 1] = { -scale, arra(0.), arra(0.), -scale };
    p[0] = { -scale, arra(0.), arra(0.), arra(scale) };
    L4v<arra> Q({ 2.*scale, arra(0.), arra(0.), arra(0.) });
    //std::cout<<Q[0]<<" "<<arra(2.*scale)<<std::endl;
    {
        arra r[3 * (n - 2) - 4];
        for (int i = 0; i < 3 * (n - 2) - 4; ++i)
        {
            r[i] = urd(eng);
        }
        array<arra, n> vm;
        for (int i = 0; i < n; ++i)
        {
            vm[i] = 0;
        }
        arra dps;
        athl::rcm(n - 2, Q, r, &vm[0], &p[1], dps);
    }
    std::vector<int> vh(n);
    for (size_t i = 0; i < n; ++i)
    {
        vh[i] = (i % 2) * 2 - 1;
    }
    acca res = amp_CO_BG(p, vh);
    std::cout << res << std::endl;
}
void benchmark_BG_CO()
{
    std::mt19937 eng(std::random_device{}());
    std::uniform_real_distribution<double> urd(0, 1);
    const int n = 9;
    //Part 1: generate a new physical kinematics
    std::vector<L4v<double>> p(n);
    //z+ corresponding to 0, z- corresponding to inf
    //we use the convention that outgoing are positive
    double scale = 1e12;
    p[n - 1] = { -scale, 0, 0, -scale };
    p[0] = { -scale, 0, 0, 1 };
    L4v<double> Q{ 2 * scale, 0, 0, 0 };
    {
        double r[3 * (n - 2) - 4];
        for (int i = 0; i < 3 * (n - 2) - 4; ++i)
        {
            r[i] = urd(eng);
        }
        array<double, n> vm;
        for (int i = 0; i < n; ++i)
        {
            vm[i] = 0;
        }
        double dps;
        athl::rcm(n - 2, Q, r, &vm[0], &p[1], dps);
    }
    std::vector<int> vh(n);
    for (uint64_t x = 0; x < (1ull << n); ++x)
    {
        for (size_t i = 0; i < n; ++i)
        {
            vh[i] = ((x >> i) & 0x1) * 2 - 1;
        }
        std::complex<double> res = amp_CO_BG(p, vh);
        std::cout << res << std::endl;
    }
}
int main()
{
    benchmark_BG_CO();
    //test_BG_CO();
    //test_BG_CO_arra();
    return 0;
}
