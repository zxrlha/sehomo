#include "BG.hpp"
#include <iostream>

#if 0
std::complex<double> amp_CO_BG(const std::vector<L4v<double>>& p, std::vector<int> h)
{
    int n = p.size();
    //first calculate all polarization vector
    std::vector<L4v<std::complex<double>>> eps;
    for (size_t i = 0; i < p.size(); ++i)
    {
        if (h[i] == -1)
        {
            eps.push_back(athl::wf_one_neg<std::complex<double>>(p[i]));
        }
        else
        {
            eps.push_back(conj(athl::wf_one_neg<std::complex<double>>(p[i])));
        }
    }
    std::vector<L4v<double>> vP(n * (n - 1) / 2);
    //we always have j>=i
    auto P = [&](int i, int j) -> L4v<double>&
    {
        return vP[j * (j + 1) / 2 + i];
    };
    //build vP
    for (int i = 0; i < n - 1; ++i)
    {
        P(i, i) = p[i];
    }
    for (int dis = 1; dis < n - 2; ++dis)
    {
        for (int i = 0; i < n - 1 - dis; ++i)
        {
            int j = i + dis;
            P(i, j) = P(i, j - 1) + p[j];
        }
    }
    std::vector<L4v<std::complex<double>>> vJ(n * (n - 1) / 2);
    //we always have j>=i
    auto J = [&](int i, int j) -> L4v<std::complex<double>>&
    {
        return vJ[j * (j + 1) / 2 + i];
    };
    //first get Jii
    for (int i = 0; i < n - 1; ++i)
    {
        J(i, i) = eps[i];
    }
    for (int dis = 1; dis < n - 2; ++dis)
    {
        for (int i = 0; i < n - 1 - dis; ++i)
        {
            int j = dis + i;
            //calculate three point vertex contribution
            for (int k = i; k <= j - 1; ++k)
            {
                const L4v<std::complex<double>>& Ja = J(i, k);
                const L4v<std::complex<double>>& Jb = J(k + 1, j);
                const L4v<double>& Pa = P(i, k);
                const L4v<double>& Pb = P(k + 1, j);
                J(i, j) += dot(Ja, Jb) * (Pa - Pb) + 2. * dot(Pb, Ja) * Jb - 2. * dot(Pa, Jb) * Ja;
            }
            //calculate four point vertex contribution
            for (int k = 0; k <= j - 2; ++k)
            {
                for (int l = k + 1; l <= j - 1; ++l)
                {
                    J(i, j) += dot(J(i, k), J(k + 1, l)) * J(l + 1, j) - 2. * dot(J(i, k), J(l + 1, j)) * J(k + 1, l) + dot(J(k + 1, l), J(l + 1, j)) * J(i, k);
                }
            }
            J(i, j) /= dot(P(i, j), P(i, j));
        }
    }
    //conclude the amplitude
    int i = 0;
    int j = n - 2;
    std::complex<double> res;
    //three point vertex contribution
    for (int k = i; k <= j - 1; ++k)
    {
        const L4v<std::complex<double>>& Ja = J(i, k);
        const L4v<std::complex<double>>& Jb = J(k + 1, j);
        const L4v<double>& Pa = P(i, k);
        const L4v<double>& Pb = P(k + 1, j);
        res += dot(Ja, Jb) * dot((Pa - Pb), eps[n - 1]) + 2. * dot(Pb, Ja) * dot(Jb, eps[n - 1]) - 2. * dot(Pa, Jb) * dot(Ja, eps[n - 1]);
    }
    //four point vertex contribution
    for (int k = 0; k <= j - 2; ++k)
    {
        for (int l = k + 1; l <= j - 1; ++l)
        {
            res += dot(J(i, k), J(k + 1, l)) * dot(J(l + 1, j), eps[n - 1]) - 2. * dot(J(i, k), J(l + 1, j)) * dot(J(k + 1, l), eps[n - 1]) + dot(J(k + 1, l), J(l + 1, j)) * dot(J(i, k), eps[n - 1]);
        }
    }
    return res;
}
#else
std::complex<double> amp_CO_BG(const std::vector<L4v<double>>& p, std::vector<int> h)
{
    int n = p.size();
    //first calculate all polarization vector
    std::vector<L4v<std::complex<double>>> eps;
    for (size_t i = 0; i < p.size(); ++i)
    {
        if (h[i] == -1)
        {
            eps.push_back(athl::wf_one_neg<std::complex<double>>(p[i]));
        }
        else
        {
            eps.push_back(conj(athl::wf_one_neg<std::complex<double>>(p[i])));
        }
    }
    std::vector<L4v<double>> vP(n * (n - 1) / 2);
    //we always have j>=i
    auto P = [&](int i, int j) -> L4v<double>&
    {
        return vP[j * (j + 1) / 2 + i];
    };
    //build vP
    for (int i = 0; i < n - 1; ++i)
    {
        P(i, i) = p[i];
    }
    for (int dis = 1; dis < n - 2; ++dis)
    {
        for (int i = 0; i < n - 1 - dis; ++i)
        {
            int j = i + dis;
            P(i, j) = P(i, j - 1) + p[j];
        }
    }
    std::vector<L4v<std::complex<double>>> vJ(n * (n - 1) / 2);
    //we always have j>=i
    auto J = [&](int i, int j) -> L4v<std::complex<double>>&
    {
        return vJ[j * (j + 1) / 2 + i];
    };
    //first get Jii
    for (int i = 0; i < n - 1; ++i)
    {
        J(i, i) = eps[i];
    }
    std::vector<std::complex<double>> vS(n * (n - 1) / 2);
    //we always have j>=i
    auto S = [&](int i, int j) -> std::complex<double>&
    {
        return vS[j * (j + 1) / 2 + i];
    };
    std::vector<athl::L4t<std::complex<double>>> vT(n * (n - 1) / 2);
    //we always have j>=i
    auto T = [&](int i, int j) -> athl::L4t<std::complex<double>>&
    {
        return vT[j * (j + 1) / 2 + i];
    };
    for (int dis = 1; dis < n - 2; ++dis)
    {
        for (int i = 0; i < n - 1 - dis; ++i)
        {
            int j = dis + i;
            //calculate three point vertex contribution
            for (int k = i; k <= j - 1; ++k)
            {
                const L4v<std::complex<double>>& Ja = J(i, k);
                const L4v<std::complex<double>>& Jb = J(k + 1, j);
                const L4v<double>& Pa = P(i, k);
                const L4v<double>& Pb = P(k + 1, j);
                J(i, j) += dot(Ja, Jb) * (Pa - Pb) + 2. * dot(Pb, Ja) * Jb - 2. * dot(Pa, Jb) * Ja;
            }
            //calculate Sij and Tij
            for (int k = i; k <= j - 1; ++k)
            {
                S(i, j) += dot(J(i, k), J(k + 1, j));
                for (int l1 = 0; l1 < 4; ++l1)
                {
                    for (int l2 = 0; l2 < 4; ++l2)
                    {
                        T(i, j)[l1][l2] += J(i, k)[l1] * J(k + 1, j)[l2];
                    }
                }
            }
            //calculate four point vertex contribution
            for (int k = i + 1; k <= j - 2; ++k)
            {
                J(i, j) += 2.0 * contract(J(k + 1, j), T(i, k))
                           - S(i, k) * J(k + 1, j) - S(k, j) * J(i, k - 1);
            }
            J(i, j) /= dot(P(i, j), P(i, j));
        }
    }
    //conclude the amplitude
    int i = 0;
    int j = n - 2;
    std::complex<double> res;
    //three point vertex contribution
    for (int k = i; k <= j - 1; ++k)
    {
        const L4v<std::complex<double>>& Ja = J(i, k);
        const L4v<std::complex<double>>& Jb = J(k + 1, j);
        const L4v<double>& Pa = P(i, k);
        const L4v<double>& Pb = P(k + 1, j);
        res += dot(Ja, Jb) * dot((Pa - Pb), eps[n - 1]) + 2. * dot(Pb, Ja) * dot(Jb, eps[n - 1]) - 2. * dot(Pa, Jb) * dot(Ja, eps[n - 1]);
    }
    //four point vertex contribution
    for (int k = 1; k <= j - 1; ++k)
    {
        res +=
            2. * dot(contract(J(k + 1, j), T(i, k)), eps[n - 1])
            - S(i, k) * dot(J(k + 1, j), eps[n - 1])
            - S(k, j) * dot(J(i, k - 1), eps[n - 1]);
    }
    return res;
}
#endif

acca amp_CO_BG(const std::vector<L4v<arra>>& p, std::vector<int> h)
{
    int n = p.size();
    //first calculate all polarization vector
    std::vector<L4v<acca>> eps;
    for (size_t i = 0; i < p.size(); ++i)
    {
        if (h[i] == -1)
        {
            eps.push_back(athl::wf_one_neg<acca>(p[i]));
        }
        else
        {
            eps.push_back(conj(athl::wf_one_neg<acca>(p[i])));
        }
    }
    std::vector<L4v<arra>> vP(n * (n - 1) / 2);
    //we always have j>=i
    auto P = [&](int i, int j) -> L4v<arra>&
    {
        return vP[j * (j + 1) / 2 + i];
    };
    //build vP
    for (int i = 0; i < n - 1; ++i)
    {
        P(i, i) = p[i];
    }
    for (int dis = 1; dis < n - 2; ++dis)
    {
        for (int i = 0; i < n - 1 - dis; ++i)
        {
            int j = i + dis;
            P(i, j) = P(i, j - 1) + p[j];
        }
    }
    std::vector<L4v<acca>> vJ(n * (n - 1) / 2);
    //we always have j>=i
    auto J = [&](int i, int j) -> L4v<acca>&
    {
        return vJ[j * (j + 1) / 2 + i];
    };
    //first get Jii
    for (int i = 0; i < n - 1; ++i)
    {
        J(i, i) = eps[i];
    }
    for (int dis = 1; dis < n - 2; ++dis)
    {
        for (int i = 0; i < n - 1 - dis; ++i)
        {
            int j = dis + i;
            //calculate three point vertex contribution
            for (int k = i; k <= j - 1; ++k)
            {
                const auto& Ja = J(i, k);
                const auto& Jb = J(k + 1, j);
                const auto& Pa = P(i, k);
                const auto& Pb = P(k + 1, j);
                J(i, j) += dot(Ja, Jb) * (Pa - Pb) + 2. * dot(Pb, Ja) * Jb - 2. * dot(Pa, Jb) * Ja;
            }
            //calculate four point vertex contribution
            for (int k = 0; k <= j - 2; ++k)
            {
                for (int l = k + 1; l <= j - 1; ++l)
                {
                    J(i, j) += dot(J(i, k), J(k + 1, l)) * J(l + 1, j) - 2. * dot(J(i, k), J(l + 1, j)) * J(k + 1, l) + dot(J(k + 1, l), J(l + 1, j)) * J(i, k);
                }
            }
            J(i, j) /= dot(P(i, j), P(i, j));
        }
    }
    //conclude the amplitude
    int i = 0;
    int j = n - 2;
    acca res;
    //three point vertex contribution
    for (int k = i; k <= j - 1; ++k)
    {
        const auto& Ja = J(i, k);
        const auto& Jb = J(k + 1, j);
        const auto& Pa = P(i, k);
        const auto& Pb = P(k + 1, j);
        res += dot(Ja, Jb) * dot((Pa - Pb), eps[n - 1]) + 2. * dot(Pb, Ja) * dot(Jb, eps[n - 1]) - 2. * dot(Pa, Jb) * dot(Ja, eps[n - 1]);
    }
    //four point vertex contribution
    for (int k = 0; k <= j - 2; ++k)
    {
        for (int l = k + 1; l <= j - 1; ++l)
        {
            res += dot(J(i, k), J(k + 1, l)) * dot(J(l + 1, j), eps[n - 1]) - 2. * dot(J(i, k), J(l + 1, j)) * dot(J(k + 1, l), eps[n - 1]) + dot(J(k + 1, l), J(l + 1, j)) * dot(J(i, k), eps[n - 1]);
        }
    }
    return res;
}
