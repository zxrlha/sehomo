#ifndef _NUMCHY_PHYIMP_HPP
#define _NUMCHY_PHYIMP_HPP 11

#include <Eigen/Dense>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <cassert>
#include <complex>
#include <vector>

class phyimp {
public:
    friend class boost::serialization::access;

    phyimp()
    {
    }

    void set(std::vector<double>&& vkin);

    void set_init_condition(int n, const std::vector<double>& v0, const std::vector<std::vector<std::complex<double> > >& v0sol)
    {
        _n = n;
        _H.resize(n - 1, n - 3);
        _M.resize(n - 1);
        _F.resize(n - 1);
        _z_gf.resize(n - 3);
        _v0 = v0;
        _v0_solutions = v0sol;
    }

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& _n;
        ar& _v0;
        ar& _v0_solutions;
    }

    double c(size_t i, size_t j) const
    {
        assert(i <= _n - 2);
        assert(j <= _n - 2);
        if (i >= j) {
            return _vc[(i - 1) * i / 2 + j];
        } else {
            return _vc[(j - 1) * j / 2 + i];
        }
    }
    double s0(size_t i, size_t j) const
    {
        assert(i <= _n - 2);
        assert(j <= _n - 2);
        if (i >= j) {
            return _v0[(i - 1) * i / 2 + j];
        } else {
            return _v0[(j - 1) * j / 2 + i];
        }
    }
    std::complex<double> z(size_t i) const
    {
        assert(i <= _n - 2);
        if (i == 0) {
            return 0;
        }
        if (i == 1) {
            return 1;
        }
        return _z_gf[i - 2];
    }

    void solve(int index, std::vector<std::complex<double>>& res);

protected:
    void set_z(const std::vector<std::complex<double> >& z);
    void calculate_H(const std::complex<double>& t);
    void calculate_M(const std::complex<double>& t);
    void calculate_F(const std::complex<double>& t);

    void calculate_dzdt(std::vector<std::complex<double> >& dzdt);

    void solve_homotopy_continuation(std::vector<std::complex<double> >& vres, double t0, double t1);
    void solve_newton_raphson(std::vector<std::complex<double> >& vres, double t);

    size_t _n;
    std::vector<double> _v0;
    std::vector<double> _vc;

    std::vector<std::vector<std::complex<double> > > _v0_solutions;

    //To simplify interface and also avoid memory reallocation, we store the temporary values
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor | Eigen::AutoAlign, 12, 12> _H;
    Eigen::VectorXcd _M;
    Eigen::VectorXcd _F;
    std::vector<std::complex<double> > _z_gf; //the (n-3) gauge-fixed z
};

#endif
