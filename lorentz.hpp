#ifndef _ATHL_LORENTZOBJECTS_HPP
#define _ATHL_LORENTZOBJECTS_HPP 1

#include <cassert>
#include <array>
#include <complex>
#include <iostream>

template<typename T>
class L4v: public std::array<T, 4>
{
public:
    L4v() = default;
    L4v(const L4v<T>&) = default;
    L4v(std::initializer_list<T>&& ilist)
    {
        std::move(ilist.begin(), ilist.end(), this->begin());
    }
    //convert from other type
    template<typename U, typename = std::enable_if_t<std::is_convertible<U, T>::value>>
    L4v(const L4v<U>& other)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            (*this)[i] = other[i];
        }
    }
    template<typename U, typename = std::enable_if_t<std::is_convertible<U, T>::value>>
    L4v<T>& operator=(const L4v<U>& other)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            (*this)[i] = other[i];
        }
        return *this;
    }

    const L4v<T>& operator+=(const L4v<T>& other)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            (*this)[i] += other[i];
        }
        return *this;
    }
    const L4v<T>& operator-=(const L4v<T>& other)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            (*this)[i] -= other[i];
        }
        return *this;
    }
    const L4v<T>& operator*=(T factor)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            (*this)[i] *= factor;
        }
        return *this;
    }
    const L4v<T>& operator/=(T factor)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            (*this)[i] /= factor;
        }
        return *this;
    }

    T pt() const
    {
        return sqrt((*this)[1] * (*this)[1] + (*this)[2] * (*this)[2]);
    }
    T phi() const
    {
        return atan2((*this)[2], (*this)[1]);
    }

    T theta() const
    {
        return atan2(pt(), (*this)[3]);
    };
};

template<typename T>
class L4t: public std::array<std::array<T, 4>, 4>
{
public:
    L4t() = default;
    L4t(const L4t<T>&) = default;
    //convert from other type
    template<typename U, typename = std::enable_if_t<std::is_convertible<U, T>::value>>
    L4t(const L4t<U>& other)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            for (size_t j = 0; j < 4; ++j)
            {
                (*this)[i][j] = other[i][j];
            }
        }
    }
    template<typename U, typename = std::enable_if_t<std::is_convertible<U, T>::value>>
    L4t<T>& operator=(const L4t<U>& other)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            for (size_t j = 0; j < 4; ++j)
            {
                (*this)[i][j] = other[i][j];
            }
        }
        return *this;
    }
    const L4t<T>& operator+=(const L4t<T>& other)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            for (size_t j = 0; j < 4; ++j)
            { (*this)[i][j] += other[i][j]; }
        }
        return *this;
    }
    const L4t<T>& operator-=(const L4t<T>& other)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            for (size_t j = 0; j < 4; ++j)
            { (*this)[i][j] -= other[i][j]; }
        }
        return *this;
    }
    const L4t<T>& operator*=(T factor)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            for (size_t j = 0; j < 4; ++j)
            { (*this)[i][j] *= factor; }
        }
        return *this;
    }
    const L4t<T>& operator/=(T factor)
    {
        for (size_t i = 0; i < 4; ++i)
        {
            for (size_t j = 0; j < 4; ++j)
            { (*this)[i][j] /= factor; }
        }
        return *this;
    }
};

template<typename T>
L4t<T> metric_tensor()
{
    L4t<T> r;
    r[0][0] = 1;
    r[1][1] = -1;
    r[2][2] = -1;
    r[3][3] = -1;
    return r;
}

template<typename T>
L4v<T> conj(const L4v<T>& a)
{
    return {conj(a[0]), conj(a[1]), conj(a[2]), conj(a[3])};
}
template<typename T>
L4t<T> conj(const L4t<T>& a)
{
    return {conj(a[0]), conj(a[1]), conj(a[2]), conj(a[3])};
}

//-------------------contract----------------------
template<typename T>
T contract(const L4v<T>& a, const L4v<T>& b)
{
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
}

template<typename T>
L4t<T> contract(const L4t<T>& a, const L4t<T>& b)
{
    L4t<T> r;
    for (size_t i = 0; i < 4; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            r[i][j] = a[i][0] * b[0][j] - a[i][1] * b[1][j] - a[i][2] * b[2][j] - a[i][3] * b[3][j];
        }
    }
    return r;
}
template<typename T>
L4v<T> contract(const L4t<T>& a, const L4v<T>& b)
{
    L4v<T> r;
    for (size_t i = 0; i < 4; ++i)
    {
        r[i] = a[i][0] * b[0] - a[i][1] * b[1] - a[i][2] * b[2] - a[i][3] * b[3];
    }
    return r;
}
template<typename T>
L4v<T> contract(const L4v<T>& a, const L4t<T>& b)
{
    L4v<T> r;
    for (size_t i = 0; i < 4; ++i)
    {
        r[i] = a[0] * b[0][i] - a[1] * b[1][i] - a[2] * b[2][i] - a[3] * b[3][i];
    }
    return r;
}
//-------------------outer product--------------------
template<typename T>
L4t<T> outer_product(const L4v<T>& a, const L4v<T>& b)
{
    L4t<T> r;
    for (size_t i = 0; i < 4; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            r[i][j] = a[i] * b[j];
        }
    }
    return r;
}
//------------------dot product----------------------
template<typename T>
T dot_product(const L4v<T>& a, const L4v<T>& b)
{
    return contract(a, b);
}

//-------------------linear operations------------------
template<typename T>
L4v<T> operator+(const L4v<T>& left, const L4v<T>& right)
{
    L4v<T> res(left);
    res += right;
    return res;
}
template<typename T>
L4v<T> operator-(const L4v<T>& left, const L4v<T>& right)
{
    L4v<T> res(left);
    res -= right;
    return res;
}
//scaling functions:
template<typename T>
L4v<T> operator*(const L4v<T>& v, const T& factor)
{
    L4v<T> res(v);
    res *= factor;
    return res;
}
template<typename T>
L4v<T> operator*(const T& factor, const L4v<T>& v)
{
    L4v<T> res(v);
    res *= factor;
    return res;
}
template<typename T>
L4v<T> operator/(const L4v<T>& v, const T& factor)
{
    L4v<T> res(v);
    res /= factor;
    return res;
}
//converting with scaling
template<typename T1, typename T2>
L4v<decltype(T1()*T2())> operator*(const L4v<T1>& v, const T2& factor)
{
    using T3 = decltype(T1() * T2());
    return L4v<T3>(v) * T3(factor);
}
template<typename T1, typename T2>
L4v<decltype(T1()*T2())> operator*(const T2& factor, const L4v<T1>& v)
{
    using T3 = decltype(T1() * T2());
    return L4v<T3>(v) * T3(factor);
}
template<typename T1, typename T2>
L4v < decltype(T1() / T2()) > operator/(const L4v<T1>& v, const T2& factor)
{
    using T3 = decltype(T1() * T2());
    return L4v<T3>(v) / T3(factor);
}

//-------------------linear operations------------------
template<typename T>
L4t<T> operator+(const L4t<T>& left, const L4t<T>& right)
{
    L4t<T> res(left);
    res += right;
    return res;
}
template<typename T>
L4t<T> operator-(const L4t<T>& left, const L4t<T>& right)
{
    L4t<T> res(left);
    res -= right;
    return res;
}
//scaling functions:
template<typename T>
L4t<T> operator*(const L4t<T>& v, const T& factor)
{
    L4t<T> res(v);
    res *= factor;
    return res;
}
template<typename T>
L4t<T> operator/(const L4t<T>& v, const T& factor)
{
    L4t<T> res(v);
    res /= factor;
    return res;
}
template<typename T>
L4t<T> operator*(const T& factor, const L4t<T>& v)
{
    L4t<T> res(v);
    res *= factor;
    return res;
}
//converting with scaling
template<typename T1, typename T2>
L4t<decltype(T1()*T2())> operator*(const L4t<T1>& v, const T2& factor)
{
    using T3 = decltype(T1() * T2());
    return L4t<T3>(v) * T3(factor);
}
template<typename T1, typename T2>
L4t<decltype(T1()*T2())> operator*(const T2& factor, const L4t<T1>& v)
{
    using T3 = decltype(T1() * T2());
    return L4t<T3>(v) * T3(factor);
}
template<typename T1, typename T2>
L4t < decltype(T1() / T2()) > operator/(const L4t<T1>& v, const T2& factor)
{
    using T3 = decltype(T1() * T2());
    return L4t<T3>(v) / T3(factor);
}

template<typename T>
void lorentz_transform(const L4v<T>& k1, const L4v<T>& k2, const L4v<T>& p1, L4v<T>& p2)
{
    T Ksq = contract(k1, k1);
    L4v<T> kksum = k1 + k2;
    p2 = p1 - T(2.) * contract(kksum, p1) / contract(kksum, kksum) * kksum + T(2.) * contract(k1, p1) / Ksq * k2;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const L4v<T>& v)
{
    for (int i = 0; i < 4; ++i)
    {
        if (i != 0) { os << " "; }
        os << v[i];
    }
    return os;
}

#endif
