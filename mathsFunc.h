#ifndef MATHSFUNC_H
#define MATHSFUNC_H
#include <vector>

double integrate(const std::vector<double> & f, double h);

std::vector<double> slice(std::vector<double> const &v, int m, int n);

template<typename T>
std::vector<T> arange(T start, T stop, T step);


#endif