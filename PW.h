#ifndef PW_H
#define PW_H

struct Parameters
{
    public:
        const double mass;
        double E;
        const double h;
        const double x0;
        const int n;
        const double k;
};

double V(double x, double k);//Simple harmonic potential

#endif