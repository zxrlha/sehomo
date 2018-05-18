#ifndef _WAMAM_MOMENTUM_HPP
#define _WAMAM_MOMENTUM_HPP 1

#include <array>

namespace wamam
{
    template<typename T>
    class momentum: public std::array<T, 4>
    {
    public:
        momentum() = default;
        momentum(const momentum<T>& other) = default;
        momentum(std::initializer_list<T>&& ilist)
        {
            std::copy(ilist.begin(), ilist.end(), this->begin());
        }
        const T& e() const { return (*this)[0]; }
        const T& px() const { return (*this)[1]; }
        const T& py() const { return (*this)[2]; }
        const T& pz() const { return (*this)[3]; }
        T& e() { return (*this)[0]; }
        T& px() { return (*this)[1]; }
        T& py() { return (*this)[2]; }
        T& pz() { return (*this)[3]; }

        T pt() const { return sqrt(px() * px() + py() * py()); }

        const momentum<T>& operator*=(T factor)
        {
            for (int i = 0; i < 4; ++i)
            {
                (*this)[i] *= factor;
            }
            return *this;
        }

        const momentum<T>& operator+=(const momentum<T>& a)
        {
            for (int i = 0; i < 4; ++i)
            {
                (*this)[i] += a[i];
            }
            return *this;
        }

        const momentum<T>& operator-=(const momentum<T>& a)
        {
            for (int i = 0; i < 4; ++i)
            {
                (*this)[i] -= a[i];
            }
            return *this;
        }
    };

    template<typename T>
    T operator*(const momentum<T>& p1, const momentum<T>& p2)
    {
        return p1.e() * p2.e() - p1.px() * p2.px() - p1.py() * p2.py() - p1.pz() * p2.pz();
    }

    template<typename T1, typename T2>
    auto operator*(const momentum<T1>& p1, const momentum<T2>& p2)->decltype(T1()*T2())
    {
        return p1.e() * p2.e() - p1.px() * p2.px() - p1.py() * p2.py() - p1.pz() * p2.pz();
    }

    template<typename T>
    momentum<T> operator+(const momentum<T>& p1, const momentum<T>& p2)
    {
        momentum<T> res(p1);
        res += p2;
        return res;
    }
    template<typename T>
    momentum<T> operator-(const momentum<T>& p1, const momentum<T>& p2)
    {
        momentum<T> res(p1);
        res -= p2;
        return res;
    }

    template<typename T>
    momentum<T> operator*(const momentum<T>& p, T factor)
    {
        momentum<T> res(p);
        res *= factor;
        return res;
    }

    template<typename T>
    momentum<T> operator*(T factor, const momentum<T>& p)
    {
        momentum<T> res(p);
        res *= factor;
        return res;
    }
}

#endif
