#include <boost/numeric/odeint.hpp>
#include <wat/solve.hpp>
#include "posimp.hpp"

namespace odeint = boost::numeric::odeint;

void posimp::init(std::vector<double>&& vkin)
{
    _v = std::move(vkin);
    double sum = 0;
    for (int i = 1; i < (_n - 2) * (_n - 3) / 2; ++i)
    {
        sum += _v[i];
    }
    _sAB0 = -sum;
}

void posimp::set_z(const std::vector<double>& z)
{
    for (size_t i = 0; i < _n - 3; ++i)
    {
        _z_gf(i) = z[i];
    }
}
void posimp::set_z(const Eigen::VectorXd& z)
{
    _z_gf = z;
}

void posimp::calculate_H(double t)
{
    //the first and second line of H
    for (size_t j = 2; j < _n - 2; ++j)
    {
        _H(0, j - 2) = -s(0, j) / (z(j) * z(j));
        _H(1, j - 2) = -s(1, j) / ((1 - z(j)) * (1 - z(j)));
    }
    _H(0, _n - 4) = - t * s(0, _n - 2) / (z(_n - 2) * z(_n - 2));
    _H(1, _n - 4) = - t * s(1, _n - 2) / ((1 - z(_n - 2)) * (1 - z(_n - 2)));
    //the middle lines
    for (size_t i = 2; i < _n - 2; ++i)
    {
        double sum = 0;
        for (size_t j = 2; j < _n - 2; ++j)
        {
            if (j != i)
            {
                double tmp = -s(i, j) / ((z(i) - z(j)) * (z(i) - z(j)));
                _H(i, j - 2) = tmp;
                sum += tmp;
            }
        }
        size_t j = _n - 2;
        double tmp = -t * s(i, j) / ((z(i) - z(j)) * (z(i) - z(j)));
        _H(i, _n - 4) = tmp;
        sum += tmp;
        _H(i, i - 2) = -sum + s(0, i) / (z(i) * z(i)) + s(1, i) / ((z(i) - 1) * (z(i) - 1));
    }
    //the last line
    size_t i = _n - 2;
    double sum = 0;
    for (size_t j = 2; j < _n - 2; ++j)
    {
        double tmp = -s(i, j) / ((z(i) - z(j)) * (z(i) - z(j)));
        _H(i, j - 2) = tmp;
        sum += tmp;
    }
    _H(i, i - 2) = -sum + s(0, i) / (z(i) * z(i)) + s(1, i) / ((z(i) - 1) * (z(i) - 1));
}

void posimp::calculate_M(double t)
{
    for (size_t i = 0; i <= _n - 3; ++i)
    {
        _M(i) = s(_n - 2, i) / (z(i) - z(_n - 2));
    }
    _M(0) += -s(0, 1) + _sAB0;
    _M(1) += s(0, 1) - _sAB0;
}

void posimp::calculate_F(double t)
{
    for (size_t i = 0; i <= _n - 3; ++i)
    {
        double sum = 0;
        for (size_t j = 0; j <= _n - 3; ++j)
        {
            if (j != i)
            {
                sum += s(i, j) / (z(i) - z(j));
            }
        }
        sum += t * s(i, _n - 2) / (z(i) - z(_n - 2));
        _F(i) = sum;
    }
    _F(0) += -(t - 1) * (s(0, 1) - _sAB0);
    _F(1) += (t - 1) * (s(0, 1) - _sAB0);
    double sum = 0;
    for (size_t j = 0; j <= _n - 3; ++j)
    {
        sum += s(_n - 2, j) / (z(_n - 2) - z(j));
    }
    _F(_n - 2) = sum;
}

void posimp::calculate_dzdt(Eigen::VectorXd& dzdt)
{
    //NOTE: using inplace decomposition wouldn't improve performance, since it is a small matrix
    dzdt = _H.colPivHouseholderQr().solve(_M);
}
void posimp::calculate_dzdt(std::vector<double>& dzdt)
{
    //NOTE: using inplace decomposition wouldn't improve performance, since it is a small matrix
    Eigen::VectorXd tmp;
    calculate_dzdt(tmp);
    dzdt.resize(_n - 3);
    for (size_t i = 0; i < _n - 3; ++i)
    {
        dzdt[i] = tmp(i);
    }
}

void posimp::solve()
{
    if (_n == 4)
    {
        _v_solutions.push_back({-s(0, 2) / s(0, 1)});
    }
    else
    {
        //first, get kinematics for _n-1 points
        std::vector<double> nv((_n - 2) * (_n - 3) / 2);
        double sum = 0;
        for (int i = 1; i < (_n - 2) * (_n - 3) / 2; ++i)
        {
            nv[i] = _v[i];
            sum += nv[i];
        }
        nv[0] = -sum;
        //build _n-1 points posimp
        posimp nps(_n - 1);
        nps.init(std::move(nv));
        nps.solve();
        //iterate over solutions
        for (size_t i = 0; i < nps.solutions().size(); ++i)
        {
            solve_through_soft_limit(nps.solutions()[i]);
        }
    }
}

void posimp::solve_through_soft_limit(const std::vector<double>& vns)
{
    //first solve the soft equation to get n-3 solutions
    std::vector<double> vss = solve_soft_equation(vns);
    for (size_t i = 0; i < vss.size(); ++i)
    {
        std::vector<double> vtmp(vns);
        vtmp.push_back(vss[i]);
        size_t nsteps = 1;
        for (size_t i = 0; i < nsteps; ++i)
        {
            solve_homotopy_continuation(vtmp, double(i) / nsteps, double(i + 1) / nsteps);
            //solve_newton_raphson(vtmp, double(i + 1) / nsteps);
        }
        _v_solutions.push_back(std::move(vtmp));
    }
}

void posimp::solve_newton_raphson(std::vector<double>& z, double t)
{
    double max_error = 1e-15;
    bool flag = true;
    while (flag)
    {
        this->set_z(z);
        this->calculate_F(t);
        this->calculate_H(t);
        Eigen::VectorXd delta = _H.fullPivLu().solve(_F);
        flag = false;
        for (size_t i = 0; i < z.size(); ++i)
        {
            z[i] += delta(i);
            if (delta(i) > max_error) { flag = true; }
        }
    }
}

void posimp::solve_homotopy_continuation(std::vector<double>& z, double t0, double t1)
{
    //using state_type = std::vector<double>;
    using state_type = Eigen::VectorXd;
    auto difffunc = [this](const state_type & z, state_type & dzdt, double t)
    {
        this->set_z(z);
        this->calculate_H(t);
        this->calculate_M(t);
        this->calculate_dzdt(dzdt);
    };
    auto integrator = odeint::make_controlled<odeint::runge_kutta_fehlberg78<state_type>>(0.0, 1e-15);
    //check initial condition
    /*
    for (size_t i = 0; i < z.size(); ++i)
    {
        std::cout << z[i] << " ";
    }
    std::cout << std::endl;
    set_z(z);
    calculate_F(0);
    std::cout << "F:";
    for (size_t i = 0; i < _F.size(); ++i)
    {
        std::cout << _F(i) << " ";
    }
    std::cout << std::endl;
    */
    //do integration
    size_t steps = odeint::integrate_adaptive(integrator, difffunc,
                   z, t0, t1, 0.001);
    //check after DE
    std::cout << steps << " ";
    for (size_t i = 0; i < z.size(); ++i)
    {
        std::cout << z[i] << " ";
    }
    std::cout << std::endl;
    /*
    set_z(z);
    calculate_F(1);
    std::cout << "F:";
    for (size_t i = 0; i < _F.size(); ++i)
    {
        std::cout << _F(i) << " ";
    }
    std::cout << std::endl;
    */
}

std::vector<double> posimp::solve_soft_equation(const std::vector<double>& vss)
{
    std::vector<double> res;
    //first, we sort
    std::vector<double> vtmp(vss);
    vtmp.push_back(0);
    vtmp.push_back(1);
    std::sort(vtmp.begin(), vtmp.end());
    for (size_t i = 0; i < vtmp.size() - 1; ++i)
    {
        auto func = [this, &vss](double z)
        {
            double sum = 0;
            sum += this->s(0, _n - 2) / (z);
            sum += this->s(1, _n - 2) / (z - 1);
            for (size_t i = 2; i < _n - 2; ++i)
            {
                sum += this->s(i, _n - 2) / (z - vss[i - 2]);
            }
            return sum;
        };
        double z = wat::solve_brent(func, vtmp[i] + 1e-7, vtmp[i + 1] - 1e-7, 1e-16);
        res.push_back(z);
    }
    return std::move(res);
}
