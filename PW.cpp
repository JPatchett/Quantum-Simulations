//Simple code for the matrix method of solving Schrodinger's equation.
// hbar = 1

#include "PW.h"
#include "mathsFunc.h"
#include "csvIO.h"
#include <vector>
#include <cmath>
#include <algorithm> 
#include <memory>
#include <functional>
#include <iostream>


//First off use Numerov's algorithm for solving second order differential equations.

inline double negSqrt(double x){return sqrt(-1.0*x);}

double V(double x, double k){ 
       /* if( x<0){
            return 0.5*k*(x+1/2)*(x+1/2);
        }
        else{
            return 0.5*k*(x-1/2)*(x-1/2);
        }*/
        return 0.5*k*x*x;
    }


void kxh(std::vector<double> & kxh, Parameters& p){

    double x = p.x0;

    //Calculate
    for(int i=0; i<kxh.size(); i++){

            kxh[i] =2*p.mass*(p.E-V(x, p.k))*pow(p.h,2)/(12.0);
            x +=  p.h;
    }
    
} 

void printArray(std::vector<double> & vec){

    for(int i = 0; i<vec.size(); i++){

        std::cout << vec[i] << std::endl;
    }
}

void psi( std::vector<double> & knf, std::vector<double> & psi, Parameters& p){ //Foreward calculate psi(x+h) using Numerov's algroithm
    
    //Create some dummy variables for ease of reading
    double num1 = 0;
    double num2 = 0;
    double denom = 1;

    //Itteratively calculate the result of Numerov's algroithm
    for(int i=1; i< psi.size()-1; i++ ){
        
        num1 = 2*psi[i]*(1-5.0*knf[i]);
        num2 = psi[i-1]*(1+knf[i-1]);
        denom = 1.0+knf[i+1];

        psi[i+1] = (num1-num2)/denom; 
    }
}

std::vector<double> wkb(std::vector<double> & knf, Parameters& p){ //Use the WKB approximation to determine first guess of the wavefunction

    int lb = 0;
    int ub = p.x0+p.h*p.n;

    //Itterate through from the beginning to find the classically forbidden region (CFR)
    for(int i=0; i<p.n; i++){

        if(knf[i] > 0 ){
            lb = i;
            break;
        }
    }
    
    for(int j=1; j<p.n;j++){

        if(knf[p.n-j] > 0){
            ub = p.n-j+1;
            break;
        }
    }

    //Take the sqrt of -kxh in the lower CFR
    std::vector<double> tempArr = slice(knf, 0, lb);
    
    std::vector<double> tempL;
    std::transform(tempArr.begin(), tempArr.end(), std::back_inserter(tempL), negSqrt); //Getting some weird itterator invalidation going on here so have to use a 2nd array
    std::vector<double> tempL2 = slice(tempL, 1, lb-1);


    //Same for the upper CFR
    std::vector<double> tempArr2 = slice(knf, ub, p.n);

    std::vector<double> tempU;
    std::transform(tempArr2.begin(),tempArr2.end(), std::back_inserter(tempU), negSqrt);
    std::vector<double> tempU2 = slice(tempU, 0, p.n-ub-1);

    //Now use the WKB approximation to calculate the intialising parameters (~log(psi)):
    double Lpsilb = p.h*exp(-sqrt(12.0)*integrate(tempL, p.h)/p.h)/sqrt(12.0*tempL[0]);  //Strictly speaking the p.h term isn't needed as the wavefunction is un-normalised anyway 
    double Lpsiub = p.h*exp(-sqrt(12.0)*integrate(tempU, p.h/p.h))/sqrt(12.0*tempU[p.n-ub-1]); //But I've left it in to be consistent with the canonical form of the WKB approximation
    double Lpsilb2 = p.h*exp(-sqrt(12.0)*integrate(tempL2, p.h)/p.h)/sqrt(12.0*tempL2[0]);
    double Lpsiub2 = p.h*exp(-sqrt(12.0)*integrate(tempU2, p.h)/p.h)/sqrt(12.0*tempU2[p.n-ub-2]);

    std::vector<double> result;

    result.push_back(Lpsilb);
    result.push_back(Lpsilb2);
    result.push_back(Lpsiub2);
    result.push_back(Lpsiub);

    return result;
}

int main()
{
    double E1 = 0.5;//Upper and lower bounds for the energy
    double E2 = 1;
    double m0 = 1; //Particle mass

    const double Ec = 0.01;
    
    //X-axis discritisation 
    double x0 = -2;
    double xF = 2;
    double  h = 0.01;
    int n = 201;
    int n2 = 401;
    double k =1;

    //Put all of our parameters in a struct
    Parameters p = {m0, E1, h, x0, n, k};
    Parameters p2 = {m0, E2, -h, xF, n, k};

    //Calculate our k^2(x) array
    std::vector<double> kvec(n2, 0.0);
    kxh(kvec, p);
    
    //Run the algorithm to determine the wavefunction
    std::vector<double> psiF( n, 0.0);
    std::vector<double> psiF2( n, 0.0);
    

    //Use the WKB approximation to determine the intial values for our wavefunctions:
    std::vector<double> wkbApprox = wkb(kvec, p);
    psiF[0] = wkbApprox[0];
    psiF[1] = wkbApprox[1];
    psiF2[n2-2] = wkbApprox[2];
    psiF2[n2-1] = wkbApprox[3];

    //Input the guesses for the wavefunction on the boundary
    psi(kvec, psiF, p);


    //Write the data to a CSV filea to be processed in python
    CSVWriter conditions("numerov.csv");
    conditions.addDataInRow(psiF);

    //Write the x position data
    //std::vector<double>  xPos = arange(p.x0, p.x0+(p.n-1)*p.h, p.h);
   // conditions.addDataInRow(xPos);

    return 0;
}
