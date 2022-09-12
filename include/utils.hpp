/*
 * Some useful utilities functions/templates,
 * reference is given for other people's codes.
 * Author: Rui-Hao Bi
 * Date: Aug 2022.
 */

#include <vector>
#ifndef MYLIB_CONSTANTS_H
#define MYLIB_CONSTANTS_H

namespace MyConst
{
  const double mass = 2000;    
}
#endif

// linspace vector implementation in c++
template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
typename std::vector<T>::iterator x;
    T val;
for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
return xs;
}

// Concatenate vectors
template <typename T>
std::vector<T> operator+(std::vector<T> const &x, std::vector<T> const &y)
{
    std::vector<T> vec;
    vec.reserve(x.size() + y.size());
    vec.insert(vec.end(), x.begin(), x.end());
    vec.insert(vec.end(), y.begin(), y.end());
    return vec;
}

template <typename T>
std::vector<T> &operator+=(std::vector<T> &x, const std::vector<T> &y)
{
    x.reserve(x.size() + y.size());
    x.insert(x.end(), y.begin(), y.end());
    return x;
}
