#include "rcm.hpp"
#include "phyimp.hpp"
#include "timer.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <boost/lexical_cast.hpp>

using namespace std;

//Calculating n-point physical kinematics, using initial condition given in file "ic-*"
void test_phyimp(int n)
{
    phyimp tpi;
    {
        //reading the initial condition
        std::ifstream ifs("ic-" + boost::lexical_cast<std::string>(n));
        boost::archive::binary_iarchive ia(ifs);
        ia >> tpi;
    }
    std::mt19937 eng(std::random_device{}());
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
    tpi.set(std::move(sij));
    std::vector<std::vector<std::complex<double>>> vres(nsol);
    timer ti;
    ti.start();
    for (size_t j = 0; j < nsol; ++j)
    {
        tpi.solve(j, vres[j]);
        ti.stop();
        std::cout << "Step 2:" << j + 1 << " solutions obtained in " << ti.time() << "s\r";
    }
    std::cout << std::endl;
    ti.stop();
    //Now check all solutions are distinct
    std::sort(vres.begin(), vres.end(), [](const std::vector<std::complex<double>>& a, const std::vector<std::complex<double>>& b)
    {
        return a[0].real() < b[0].real();
    });
    double maxdiff = 0.0;
    for (int j = 0; j < nsol; j += 2)
    {
        maxdiff = std::max(maxdiff, std::abs(vres[j][0].real() - vres[j + 1][0].real()));
        maxdiff = std::max(maxdiff, std::abs(vres[j][0].imag() + vres[j + 1][0].imag()));
    }
    std::cout << "COMPLEX CONJUGATE CHECK: " << maxdiff << std::endl;
    double mindiff = 1e10;
    int mini;
    for (int j = 0; j < nsol - 2; j += 2)
    {
        double diff = std::abs(vres[j][0].real() - vres[j + 2][0].real());
        if (diff < mindiff)
        {
            mindiff = diff;
            mini = j;
        }
    }
    std::cout << "minimal diff: " << mindiff << std::endl;
	std::cout.precision(17);
    for (int i = 0; i < n - 3; ++i)
    {
        std::cout << vres[mini][i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < n - 3; ++i)
    {
        std::cout << vres[mini + 2][i] << " ";
    }
    std::cout << std::endl;
}
int main()
{
    //Calculating a random chosen 9-point physical kinematics corresponding to 2->7 scattering
    test_phyimp(9);
    return 0;
}
