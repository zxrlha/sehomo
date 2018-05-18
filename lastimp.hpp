#ifndef _NUMCHY_LASTIMP_HPP
#define _NUMCHY_LASTIMP_HPP 11

#include <cassert>
#include <vector>
#include <complex>
#include <Eigen/Dense>

class lastimp
{
public:
    lastimp(size_t n)
        : _H(n - 1, n - 3), _M(n - 1), _F(n - 1), _z_gf(n - 3),  _v0((n - 1) * (n - 2) / 2), _vc((n - 1) * (n - 2) / 2)
    {
        _n = n;
    }

    void init(std::vector<double>&& vkin, std::vector<double>&& v0kin);

    double c(size_t i, size_t j) const
    {
        assert(i <= _n - 2);
        assert(j <= _n - 2);
        if (i >= j)
        {
            return _vc[(i - 1) * i / 2 + j];
        }
        else
        {
            return _vc[(j - 1) * j / 2 + i];
        }
    }
    double s0(size_t i, size_t j) const
    {
        assert(i <= _n - 2);
        assert(j <= _n - 2);
        if (i >= j)
        {
            return _v0[(i - 1) * i / 2 + j];
        }
        else
        {
            return _v0[(j - 1) * j / 2 + i];
        }
    }
    std::complex<double> z(size_t i) const
    {
        assert(i <= _n - 2);
        if (i == 0) { return 0; }
        if (i == 1) { return 1; }
        return _z_gf[i - 2];
    }

    void solve();

    const std::vector<std::vector<std::complex<double>>>& solutions() const { return _v_solutions; }
protected:
    void set_z(const std::vector<std::complex<double>>& z);
    void calculate_H(const std::complex<double>& t);
    void calculate_M(const std::complex<double>& t);
    void calculate_F(const std::complex<double>& t);

    void calculate_dzdt(std::vector<std::complex<double>>& dzdt);

    void solve_homotopy_continuation(std::vector<std::complex<double>>& vres, const std::complex<double>& t0, const std::complex<double>& t1);
    void solve_newton_raphson(std::vector<std::complex<double>>& vres, const std::complex<double>& t);

    size_t _n;
    std::vector<double> _v0;
    std::vector<double> _vc;

    std::vector<std::complex<double>> _vmt;//middle point of p

    std::vector<std::vector<std::complex<double>>> _v_solutions;

    //To simplify interface and also avoid memory reallocation, we store the temporary values
    Eigen::Matrix < std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor | Eigen::AutoAlign, 12, 12 > _H;
    Eigen::VectorXcd _M;
    Eigen::VectorXcd _F;
    std::vector<std::complex<double>> _z_gf;//the (n-3) gauge-fixed z
};

#endif
