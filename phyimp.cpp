#include "phyimp.hpp"
#include <boost/numeric/odeint.hpp>

namespace odeint = boost::numeric::odeint;

using namespace std::literals::complex_literals;

void phyimp::set(std::vector<double>&& vkin)
{
    std::cout<<_vc.size()<<" "<<_v0.size()<<std::endl;
    _H.resize(_n - 1, _n - 3);
    _M.resize(_n - 1);
    _F.resize(_n - 1);
    _z_gf.resize(_n - 3);
    //we try to matching the overall scale
    //using L1 norm
    double ss = 0;
    double ss0 = 0;
    for (size_t i = 0; i < vkin.size(); ++i) {
        ss += std::abs(vkin[i]);
        ss0 += std::abs(_v0[i]);
    }
    std::cout<<ss<<" "<<ss0<<std::endl;
    double r = ss0 / ss;
    _vc = std::move(vkin);
    for (size_t i = 0; i < _vc.size(); ++i) {
        _vc[i] *= r;
        std::cout<<_vc[i]<<" "<<r<<" "<<_v0[i]<<std::endl;
        _vc[i] -= _v0[i];
    }
}

void phyimp::set_z(const std::vector<std::complex<double> >& z)
{
    for (size_t i = 0; i < _n - 3; ++i) {
        _z_gf[i] = z[i];
    }
}

void phyimp::calculate_H(const std::complex<double>& t)
{
    for (size_t i = 0; i <= _n - 2; ++i) {
        std::complex<double> sum = 0;
        for (size_t j = 2; j <= _n - 2; ++j) {
            if (j != i) {
                std::complex<double> tmp = -(s0(i, j) + t * c(i, j)) / ((z(i) - z(j)) * (z(i) - z(j)));
                _H(i, j - 2) = tmp;
                sum += tmp;
            }
        }
        if (i >= 2) {
            _H(i, i - 2) = -sum + (s0(0, i) + t * c(0, i)) / (z(i) * z(i)) + (s0(1, i) + t * c(1, i)) / ((z(i) - 1.) * (z(i) - 1.));
        }
    }
}

void phyimp::calculate_M(const std::complex<double>& t)
{
    for (size_t i = 0; i <= _n - 2; ++i) {
        std::complex<double> sum = 0;
        for (size_t j = 0; j <= _n - 2; ++j) {
            if (j != i) {
                sum += c(j, i) / (z(i) - z(j));
            }
        }
        _M(i) = sum;
    }
}

void phyimp::calculate_F(const std::complex<double>& t)
{
    for (size_t i = 0; i <= _n - 2; ++i) {
        std::complex<double> sum = 0;
        for (size_t j = 0; j <= _n - 2; ++j) {
            if (j != i) {
                sum += (s0(i, j) + t * c(i, j)) / (z(i) - z(j));
            }
        }
        _F(i) = sum;
    }
}

void phyimp::calculate_dzdt(std::vector<std::complex<double> >& dzdt)
{
    //NOTE: using inplace decomposition wouldn't improve performance, since it is a small matrix
    //Eigen::VectorXcd tmp = _H.fullPivHouseholderQr().solve(_M);
    Eigen::VectorXcd tmp = _H.colPivHouseholderQr().solve(_M);
    //auto adjH = _H.adjoint();
    //Eigen::VectorXcd tmp = (adjH * _H).fullPivLu().solve(adjH * _M);
    //Eigen::VectorXcd tmp = (adjH * _H).llt().solve(adjH * _M);
    dzdt.resize(_n - 3);
    for (size_t i = 0; i < _n - 3; ++i) {
        dzdt[i] = tmp(i);
    }
}

void phyimp::solve(int index, std::vector<std::complex<double> >& sol)
{
    sol = _v0_solutions[index];
    solve_homotopy_continuation(sol, 0.0, 1.0);
    solve_newton_raphson(sol, 1.0);
}

void phyimp::solve_newton_raphson(std::vector<std::complex<double> >& z, double t)
{
    double max_error = 1e-15;
    bool flag = true;
    while (flag) {
        this->set_z(z);
        this->calculate_F(t);
        this->calculate_H(t);
        Eigen::VectorXcd delta = _H.fullPivLu().solve(_F);
        flag = false;
        for (size_t i = 0; i < z.size(); ++i) {
            z[i] += delta(i);
            if (std::abs(delta(i)) / (std::abs(z[i]) + std::abs(delta(i))) > max_error) {
                flag = true;
            }
        }
    }
}

void phyimp::solve_homotopy_continuation(std::vector<std::complex<double> >& z, double t0, double t1)
{
    using state_type = std::vector<std::complex<double> >;
    auto difffunc = [this](const state_type& z, state_type& dzdt, double t) {
        this->set_z(z);
        this->calculate_H(t);
        this->calculate_M(t);
        this->calculate_dzdt(dzdt);
    };
    auto integrator = odeint::make_controlled<odeint::runge_kutta_fehlberg78<state_type> >(0.0, 1e-15);
//check initial condition
#if 0
    for (size_t i = 0; i < z.size(); ++i)
    {
        std::cout << z[i] << " ";
    }
    std::cout << std::endl;
    set_z(z);
    calculate_F(t0);
    std::cout << "F:";
    for (size_t i = 0; i < _F.size(); ++i)
    {
        std::cout << _F(i) << " ";
    }
    std::cout << std::endl;
#endif
    //do integration
    size_t steps = odeint::integrate_adaptive(integrator, difffunc,
        z, t0, t1, 0.001);
//check after DE
#if 1
    std::cout << steps << " ";
    for (size_t i = 0; i < z.size(); ++i) {
        std::cout << z[i] << " ";
    }
    std::cout << std::endl;
    set_z(z);
    calculate_F(t1);
    std::cout << "F:";
    for (size_t i = 0; i < _F.size(); ++i) {
        std::cout << _F(i) << " ";
    }
    std::cout << std::endl;
#endif
}
