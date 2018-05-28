#include <boost/numeric/odeint.hpp>
#include "comimp.hpp"
#include "posimp.hpp"

namespace odeint = boost::numeric::odeint;

using namespace std::literals::complex_literals;

void comimp::init(std::vector<double>&& vkin, std::vector<double>&& v0kin)
{
    _vmt.clear();
    _vmt.push_back(0.0);
    _vmt.push_back(std::complex<double>(0.5,0.5));
    _vmt.push_back(1.0);
    _v0 = std::move(v0kin);
    //we try to matching the overall scale
    //using L1 norm
    double ss = 0;
    double ss0 = 0;
    for (size_t i = 0; i < vkin.size(); ++i)
    {
        ss += std::abs(vkin[i]);
        ss0 += std::abs(_v0[i]);
    }
    //std::cout<<ss<<" "<<ss0<<std::endl;
    double r = ss0 / ss;
    _vc = std::move(vkin);
    for (size_t i = 0; i < _vc.size(); ++i)
    {
        _vc[i] *= r;
        //std::cout<<_vc[i]<<" "<<r<<" "<<_v0[i]<<std::endl;
        _vc[i] -= _v0[i];
    }
}

void comimp::set_z(const std::vector<std::complex<double>>& z)
{
    for (size_t i = 0; i < _n - 3; ++i)
    {
        _z_gf[i] = z[i];
    }
}

void comimp::calculate_H(const std::complex<double>& t)
{
    for (size_t i = 0; i <= _n - 2; ++i)
    {
        std::complex<double> sum = 0;
        for (size_t j = 2; j <= _n - 2; ++j)
        {
            if (j != i)
            {
                std::complex<double> tmp = -(s0(i, j) + t * c(i, j)) / ((z(i) - z(j)) * (z(i) - z(j)));
                _H(i, j - 2) = tmp;
                sum += tmp;
            }
        }
        if (i >= 2)
        {
            _H(i, i - 2) = -sum + (s0(0, i) + t * c(0, i)) / (z(i) * z(i)) + (s0(1, i) + t * c(1, i)) / ((z(i) - 1.) * (z(i) - 1.));
        }
    }
}

void comimp::calculate_M(const std::complex<double>& t)
{
    for (size_t i = 0; i <= _n - 2; ++i)
    {
        std::complex<double> sum = 0;
        for (size_t j = 0; j <= _n - 2; ++j)
        {
            if (j != i)
            {
                sum += c(j, i) / (z(i) - z(j));
            }
        }
        _M(i) = sum;
    }
}

void comimp::calculate_F(const std::complex<double>& t)
{
    for (size_t i = 0; i <= _n - 2; ++i)
    {
        std::complex<double> sum = 0;
        for (size_t j = 0; j <= _n - 2; ++j)
        {
            if (j != i)
            {
                sum += (s0(i, j) + t * c(i, j)) / (z(i) - z(j));
            }
        }
        _F(i) = sum;
    }
}

void comimp::calculate_dzdt(std::vector<std::complex<double>>& dzdt)
{
    //NOTE: using inplace decomposition wouldn't improve performance, since it is a small matrix
    Eigen::VectorXcd tmp = _H.fullPivHouseholderQr().solve(_M);
    //Eigen::VectorXcd tmp = _H.colPivHouseholderQr().solve(_M);
    //auto adjH = _H.adjoint();
    //Eigen::VectorXcd tmp = (adjH * _H).fullPivLu().solve(adjH * _M);
    //Eigen::VectorXcd tmp = (adjH * _H).llt().solve(adjH * _M);
    dzdt.resize(_n - 3);
    for (size_t i = 0; i < _n - 3; ++i)
    {
        dzdt[i] = tmp(i);
    }
}

void comimp::solve()
{
    posimp nps(_n);
    nps.init(std::vector<double>(_v0));
    nps.solve();
    //iterate over solutions
    for (size_t i = 0; i < nps.solutions().size(); ++i)
    {
        std::vector<std::complex<double>> sol;
        for (size_t j = 0; j < nps.solutions()[i].size(); ++j)
        {
            sol.push_back(nps.solutions()[i][j]);
        }
        for (size_t mi = 0; mi < _vmt.size() - 1; ++mi)
        {
            solve_homotopy_continuation(sol, _vmt[mi], _vmt[mi + 1]);
            solve_newton_raphson(sol, _vmt[mi + 1]);
        }
        _v_solutions.push_back(sol);
    }
}

void comimp::solve_newton_raphson(std::vector<std::complex<double>>& z, const std::complex<double>& t)
{
    double max_error = 1e-15;
    bool flag = true;
    while (flag)
    {
        this->set_z(z);
        this->calculate_F(t);
        this->calculate_H(t);
        Eigen::VectorXcd delta = _H.fullPivLu().solve(_F);
        flag = false;
        for (size_t i = 0; i < z.size(); ++i)
        {
            z[i] += delta(i);
            if (std::abs(delta(i)) / (std::abs(z[i]) + std::abs(delta(i))) > max_error) { flag = true; }
        }
    }
}

void comimp::solve_homotopy_continuation(std::vector<std::complex<double>>& z, const std::complex<double>& t0, const std::complex<double>& t1)
{
    using state_type = std::vector<std::complex<double>>;
    auto difffunc = [this, t0, t1](const state_type & z, state_type & dzdt, double tpara)
    {
        this->set_z(z);
        std::complex<double> t = t0 + (t1 - t0) * tpara;
        this->calculate_H(t);
        this->calculate_M(t);
        this->calculate_dzdt(dzdt);
        for (size_t i = 0; i < dzdt.size(); ++i)
        {
            dzdt[i] *= t1 - t0;
        }
    };
    auto integrator = odeint::make_controlled<odeint::runge_kutta_fehlberg78<state_type>>(0.0, 1e-15);
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
                   z, 0.0, 1.0, 0.001);
    //check after DE
#if 1
    std::cout << steps << " ";
    for (size_t i = 0; i < z.size(); ++i)
    {
        std::cout << z[i] << " ";
    }
    std::cout << std::endl;
    set_z(z);
    calculate_F(t1);
    std::cout << "F:";
    for (size_t i = 0; i < _F.size(); ++i)
    {
        std::cout << _F(i) << " ";
    }
    std::cout << std::endl;
#endif
}
