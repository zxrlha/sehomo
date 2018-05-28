#ifndef _COIP_RCM_HPP
#define _COIP_RCM_HPP 11

#include <cmath>
#include <cassert>
#include "lorentz.hpp"

template<typename T>
T rcm_pmax(const L4v<T>& Q)
{
    T Esq = contract(Q, Q);
    return sqrt(Esq) / T(2.0);
}

template<typename T>
void rcm_2(const L4v<T>& Q, T x[], L4v<T> p[], T& jkb)
{
    T pi(M_PI);
    T theta = x[0] * pi;
    T phi = x[1] * 2. * pi;
    T pmax = rcm_pmax(Q);
    T E1CM = pmax;
    T E2CM = pmax;
    T ECM = sqrt(contract(Q, Q));
    jkb = pmax * sin(theta) / 8. / ECM;
    L4v<T> p1cm;
    p1cm[0] = E1CM;
    p1cm[1] = pmax * sin(theta) * cos(phi);
    p1cm[2] = pmax * sin(theta) * sin(phi);
    p1cm[3] = pmax * cos(theta);
    L4v<T> QCM{ECM, T(), T(), T()};
    lorentz_transform(QCM, Q, p1cm, p[0]);
    p[1] = Q - p[0];
}

//The RCM method
template<typename T>
void rcm(int n, const L4v<T>& Q, T x[],  L4v<T> p[], T& jkb)
{
    assert(n >= 2);
    if (n == 2)
    {
        return rcm_2(Q, x, p, jkb);
    }
    else
    {
        T pi(M_PI);
        T theta = pi * x[0];
        T phi = 2. * pi * x[1];
        T pmax = rcm_pmax(Q);
        T p1 = pmax * x[2];
        L4v<T> p1cm;
        p1cm[0] = p1;
        p1cm[1] = p1 * sin(theta) * cos(phi);
        p1cm[2] = p1 * sin(theta) * sin(phi);
        p1cm[3] = p1 * cos(theta);
        L4v<T> QCM{sqrt(contract(Q, Q)), T(), T(), T()};
        lorentz_transform(QCM, Q, p1cm, p[0]);
        L4v<T> Qremain = Q - p[0];
        rcm(n - 1, Qremain, x + 3, p + 1, jkb);
        jkb *= pmax * pmax * pmax * x[2] * x[2] * sin(theta) / 8 / pi / p1cm[0];
    }
}

#endif
