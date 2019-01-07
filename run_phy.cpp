#include "rcm.hpp"
#include "phyimp.hpp"
#include "timer.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <boost/lexical_cast.hpp>

using namespace std;

void test_phyimp(int n)
{
    phyimp tpi;
    {
        std::ifstream ifs("ic-" + boost::lexical_cast<std::string>(n));
        boost::archive::binary_iarchive ia(ifs);
        ia >> tpi;
    }
    std::mt19937 eng(std::random_device{}());
    //std::mt19937 eng;
    std::uniform_real_distribution<double> urd(0, 1);
    //Part 1: generate a new physical kinematics
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
    int nsol = 1;
    for (size_t j = 1; j <= n - 3; ++j)
    { nsol *= j; }
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
    //std::cout << "SUM:" << sum << std::endl;
    tpi.set(std::move(sij));
    std::vector<std::complex<double>> res;
    timer ti;
    ti.start();
    for (size_t j = 0; j < nsol; ++j)
    {
        tpi.solve(j, res);
    }
    ti.stop();
    std::cout<<"total: "<<ti.time()<<"s"<<std::endl;
}
int main()
{
    test_phyimp(9);
    return 0;
}
