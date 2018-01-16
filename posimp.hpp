#ifndef _NUMCHY_POSIMP_HPP
#define _NUMCHY_POSIMP_HPP 11

#include <cassert>
#include <vector>
#include <Eigen/Dense>

class posimp
{
public:
    posimp(size_t n)
        : _H(n - 1, n - 3), _M(n - 1), _F(n - 1), _z_gf(n - 3), _v((n - 1) * (n - 2) / 2)
    {
        _n = n;
    }

    void init(std::vector<double>&& vkin);

    double s(size_t i, size_t j) const
    {
        assert(i <= _n - 2);
        assert(j <= _n - 2);
        if (i >= j)
        {
            return _v[(i - 1) * i / 2 + j];
        }
        else
        {
            return _v[(j - 1) * j / 2 + i];
        }
    }
    double z(size_t i) const
    {
        assert(i <= _n - 2);
        if (i == 0) { return 0; }
        if (i == 1) { return 1; }
        return _z_gf[i - 2];
    }

    void solve();

    const std::vector<std::vector<double>>& solutions() const { return _v_solutions; }
protected:
    void set_z(const std::vector<double>& z);
    void calculate_H(double t);
    void calculate_M(double t);
    void calculate_F(double t);

    void calculate_dzdt(std::vector<double>& dzdt);

    void solve_through_soft_limit(const std::vector<double>& vns);
    std::vector<double> solve_soft_equation(const std::vector<double>& vns);
    void solve_homotopy_continuation(std::vector<double>& vres, double t0, double t1);
    void solve_newton_raphson(std::vector<double>& vres, double t);

    size_t _n;
    std::vector<double> _v;
    double _sAB0;

    std::vector<std::vector<double>> _v_solutions;

    //To simplify interface and also avoid memory reallocation, we store the temporary values
    Eigen::Matrix < double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor | Eigen::AutoAlign, 12, 12 > _H;
    //Eigen::MatrixXd _H;
    Eigen::VectorXd _M;
    Eigen::VectorXd _F;
    std::vector<double> _z_gf;//the (n-2) gauge-fixed z
};

#endif
