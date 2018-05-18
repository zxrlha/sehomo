#include "finalstate.hpp"
#include "lastimp.hpp"
#include "phyimp.hpp"
#include "posimp.hpp"
#include <fstream>
#include <iostream>
#include <random>

using namespace std;
using wamam::momentum;

void test_posimp()
{
    //std::mt19937 eng(std::random_device{}());
    std::mt19937 eng;
    std::uniform_real_distribution<double> urd(0, 1);
    int n = 10;
    //cin>>n;
    vector<double> v0((n - 1) * (n - 2) / 2);
    for (int i = 1; i < v0.size(); ++i) {
        v0[i] = urd(eng);
    }
    //we require s(A,2)>s(A,3)>....>s(A,n-2), corresponding to v[i(i-1)/2] > v(j(j-1)/2] if i > j
    {
        vector<double> vr;
        for (int i = 2; i <= n - 2; ++i) {
            vr.push_back(urd(eng));
        }
        std::sort(vr.begin(), vr.end());
        for (size_t i = 2; i <= n - 2; ++i) {
            v0[i * (i - 1) / 2] = vr[vr.size() + 1 - i];
        }
        //conclude v0[0]
        double sum = 0;
        for (int i = 1; i < v0.size(); ++i) {
            sum += v0[i];
        }
        v0[0] = -sum;
    }
    posimp ps(n);
    ps.init(std::move(v0));
    ps.solve();
}
void test_lastimp()
{
    //std::mt19937 eng(std::random_device{}());
    std::mt19937 eng;
    std::uniform_real_distribution<double> urd(0, 1);
    const int n = 9;
    //cin>>n;
    //Part 1: generate the K+ kinematics
    vector<double> v0((n - 1) * (n - 2) / 2);
    for (int i = 1; i < v0.size(); ++i) {
        v0[i] = urd(eng) * 0.1 + 0.5;
    }
    //we require s(A,2)>s(A,3)>....>s(A,n-2), corresponding to v[i(i-1)/2] > v(j(j-1)/2] if i > j
    {
#if 0
        vector<double> vr;
        for (int i = 2; i <= n - 2; ++i)
        {
            vr.push_back(urd(eng)*0.1+0.5);
        }
        std::sort(vr.begin(), vr.end());
        for (size_t i = 2; i <= n - 2; ++i)
        {
            v0[i * (i - 1) / 2] = vr[vr.size() + 1 - i];
        }
#endif
        //conclude v0[0]
        double sum = 0;
        for (int i = 1; i < v0.size(); ++i) {
            sum += v0[i];
        }
        v0[0] = -sum;
    }
    //Part 2: generate the physical kinematics
    array<momentum<double>, n> p;
    auto pswap_ij = [&p](size_t i, size_t j) {
        std::swap(p[i], p[j]);
    };
    //z+ corresponding to 0, z- corresponding to inf
    //we use the convention that outgoing are positive
    p[n - 1] = { -1, 0, 0, -1 };
    p[0] = { -1, 0, 0, 1 };
    momentum<double> Q{ 2, 0, 0, 0 };
    while (true) {
        double r[3 * (n - 2) - 4];
        for (int i = 0; i < 3 * (n - 2) - 4; ++i) {
            r[i] = urd(eng);
        }
        array<double, n> vm;
        for (int i = 0; i < n; ++i) {
            vm[i] = 0;
        }
        double dps;
        wamam::rcm(n - 2, Q, r, &vm[0], &p[1], dps);
        bool flag = true;
        for (size_t i = 0; i < p.size(); ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (std::abs(p[i] * p[j]) < 0.001) {
                    flag = false;
                }
            }
        }
        if (flag)
            break;
    }
//let's change the order
#if 1
    {
        //Step 1: choose B from 1~n-2 satisfy s_{BC} > s_{iC}
        size_t mi = 0;
        size_t maxs = 0;
        for (size_t i = 1; i <= n - 2; ++i) {
            double s = -(p[i] * p[n - 1]);
            std::cerr << s << " " << i << std::endl;
            if (s > maxs) {
                maxs = s;
                mi = i;
            }
        }
        pswap_ij(mi, 1);
        //Step 2: sort the remain according to s_Ai
        std::sort(&p[1], &p[n - 1], [&](const momentum<double>& pa, const momentum<double>& pb) {
            return (-(pa * p[0])) >= (-(pb * p[0]));
        });
    }
#endif
    vector<double> sij;
    double sum = 0;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < i; ++j) {
            double tmp = 2.0 * (p[i] * p[j]);
            sij.push_back(tmp);
            sum += tmp;
        }
    }
    std::cout << "SUM:" << sum << std::endl;
    lastimp lps(n);
    auto tv = sij;
    lps.init(std::move(sij), std::move(v0));
    lps.solve();
    phyimp tp;
    tp.set_init_condition(n, tv, lps.solutions());
    std::ofstream ofs("SOL");
    boost::archive::text_oarchive oa(ofs);
    oa << tp;
}
void test_phyimp()
{
    phyimp tpi;
    {
        std::ifstream ifs("SOL");
        boost::archive::text_iarchive ia(ifs);
        ia >> tpi;
    }
    std::mt19937 eng(std::random_device{}());
    //std::mt19937 eng;
    std::uniform_real_distribution<double> urd(0, 1);
    const int n = 9;
    //Part 1: generate a new physical kinematics
    array<momentum<double>, n> p;
    auto pswap_ij = [&p](size_t i, size_t j) {
        std::swap(p[i], p[j]);
    };
    //z+ corresponding to 0, z- corresponding to inf
    //we use the convention that outgoing are positive
    p[n - 1] = { -1, 0, 0, -1 };
    p[0] = { -1, 0, 0, 1 };
    momentum<double> Q{ 2, 0, 0, 0 };
    while (true) {
        double r[3 * (n - 2) - 4];
        for (int i = 0; i < 3 * (n - 2) - 4; ++i) {
            r[i] = urd(eng);
        }
        array<double, n> vm;
        for (int i = 0; i < n; ++i) {
            vm[i] = 0;
        }
        double dps;
        wamam::rcm(n - 2, Q, r, &vm[0], &p[1], dps);
        bool flag = true;
        for (size_t i = 0; i < p.size(); ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (std::abs(p[i] * p[j]) < 0.01) {
                    flag = false;
                }
            }
        }
        if (flag)
            break;
    }
//let's change the order
#if 1
    {
        //Step 1: choose B from 1~n-2 satisfy s_{BC} > s_{iC}
        size_t mi = 0;
        size_t maxs = 0;
        for (size_t i = 1; i <= n - 2; ++i) {
            double s = -(p[i] * p[n - 1]);
            std::cerr << s << " " << i << std::endl;
            if (s > maxs) {
                maxs = s;
                mi = i;
            }
        }
        pswap_ij(mi, 1);
        //Step 2: sort the remain according to s_Ai
        std::sort(&p[1], &p[n - 1], [&](const momentum<double>& pa, const momentum<double>& pb) {
            return (-(pa * p[0])) >= (-(pb * p[0]));
        });
    }
#endif
    vector<double> sij;
    double sum = 0;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < i; ++j) {
            double tmp = 2.0 * (p[i] * p[j]);
            sij.push_back(tmp);
            sum += tmp;
        }
    }
    std::cout << "SUM:" << sum << std::endl;
    int nsol = 1;
    for (size_t j = 1; j <= n - 3; ++j)
        nsol *= j;
    tpi.set(std::move(sij));
    std::vector<std::complex<double> > res;
    for (size_t j = 0; j < nsol; ++j) {
        tpi.solve(j, res);
    }
}
int main()
{
    test_phyimp();
    return 0;
}
