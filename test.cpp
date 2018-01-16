#include <iostream>
#include <random>
#include "posimp.hpp"

using namespace std;

int main()
{
    //std::mt19937 eng(std::random_device{}());
    std::mt19937 eng;
    std::uniform_real_distribution<double> urd(0, 1);
    int n = 10;
    //cin>>n;
    double sum = 0;
    vector<double> v((n - 1) * (n - 2) / 2);
    for (int i = 1; i < v.size(); ++i)
    {
        v[i] = urd(eng);
        sum += v[i];
    }
    v[0] = -sum;
    posimp ps(n);
    ps.init(std::move(v));
    ps.solve();
    return 0;
}
