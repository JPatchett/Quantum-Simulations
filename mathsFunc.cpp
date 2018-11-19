#include "mathsFunc.h"

double integrate(const std::vector<double> & f, double h){ //Computer simspon's rule

    double sum = 0;
    int n = f.size();

    sum += f[0]+f[n-1];

    for(int i=1; i<n-1; i++){
        
        if(i%2 == 0){
            sum += 2.0*f[i];
        }
        else{
            sum += 4.0*f[i];
        }
    }

    sum /= 3.0;
    sum *= h;

    return sum;
}


std::vector<double> slice(std::vector<double> const &v, int m, int n){//Returns a slice of a vector
    
    auto first = v.cbegin()+m;
    auto last = v.cbegin()+n;

    std::vector<double> vec(first, last);
    
    return vec;
}

template<typename T>
std::vector<T> arange(T start, T stop, T step){
    
    std::vector<T> values;
    for(T value = start; value<stop; value += step){
        values.push_back(value);
    }

    return values;
}