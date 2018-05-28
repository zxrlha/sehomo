#include "lorentz.hpp"
#include "rcm.hpp"
#include "comimp.hpp"
#include "posimp.hpp"
#include "phyimp.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <boost/lexical_cast.hpp>

using namespace std;

void prep_ic(int n)
{
    //std::mt19937 eng(std::random_device{}());
    std::mt19937 eng;
    std::uniform_real_distribution<double> urd(0, 1);
    //Part 1: generate the K+ kinematics
    vector<double> v0((n - 1) * (n - 2) / 2);
    for (int i = 1; i < v0.size(); ++i)
    {
        v0[i] = urd(eng) * 0.3 + 0.5;
    }
    {
        //conclude v0[0]
        double sum = 0;
        for (int i = 1; i < v0.size(); ++i)
        {
            sum += v0[i];
        }
        v0[0] = -sum;
    }
    //Part 2: generate the physical kinematics
    std::vector<L4v<double>> p(n);
    //z+ corresponding to 0, z- corresponding to inf
    //we use the convention that outgoing are positive
    p[n - 1] = { -1, 0, 0, -1 };
    p[0] = { -1, 0, 0, 1 };
    L4v<double> Q{ 2, 0, 0, 0 };
    while (true)
    {
        double r[3 * (n - 2) - 4];
        for (int i = 0; i < 3 * (n - 2) - 4; ++i)
        {
            r[i] = urd(eng);
        }
        vector<double> vm(n);
        for (int i = 0; i < n; ++i)
        {
            vm[i] = 0;
        }
        double dps;
        rcm(n - 2, Q, r, &p[1], dps);
        double minv = 1;
        for (size_t i = 0; i < p.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                //E_1E_2(1-cos(theta))=p_1 dot p_2
                minv = std::min(minv, (std::abs(dot_product(p[i], p[j]) / p[i][0] / p[j][0])));
            }
        }
        if (minv > 0.2) { break; }
    }
    vector<double> sij;
    double sum = 0;
    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            double tmp = 2.0 * dot_product(p[i], p[j]);
            sij.push_back(tmp);
            sum += tmp;
        }
    }
    std::cout << "SUM:" << sum << std::endl;
    comimp lps(n);
    auto tv = sij;
    lps.init(std::move(sij), std::move(v0));
    lps.solve();
    phyimp tp;
    tp.set_init_condition(n, tv, lps.solutions());
    std::ofstream ofs("ic-" + boost::lexical_cast<std::string>(n));
    boost::archive::text_oarchive oa(ofs);
    oa << tp;
}
int main()
{
    for (int n = 5; n <= 9; ++n)
    prep_ic(n);
    return 0;
}
