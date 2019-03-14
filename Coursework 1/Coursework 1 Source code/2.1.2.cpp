#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

double d1(double S, double X, double r, double D0, double sigma, double T, double t)
{
	return (log(S/X) + (r - D0 + 0.5 *sigma*sigma)*(T - t)) / (sigma * sqrt(T -t));
}

double d2(double S, double X, double r, double D0, double sigma, double T, double t)
{
	return d1(S,X,r,D0,sigma,T,t) - sigma * sqrt(T - t);
}

double normalDistribution_builtin(double x)
{
	return 0.5*erfc(-x / sqrt(2.));
}

double P(double S, double X, double r, double D0, double sigma, double T, double t)
{
	return X*exp(-1*r*(T-t))* normalDistribution_builtin(-d2(S,X,r,D0,sigma,T,t)) - S*exp(-1*D0*(T-t))* normalDistribution_builtin(-d1(S,X,r,D0,sigma,T,t));

}

double C(double S, double X, double r, double D0, double sigma, double T, double t)
{
	return S*exp(-1*D0*(T-t))* normalDistribution_builtin(d1(S,X,r,D0,sigma,T,t)) - X*exp(-1*r*(T-t))* normalDistribution_builtin(d2(S,X,r,D0,sigma,T,t));
}

double BC(double S, double X, double r, double D0, double sigma, double T, double t)
{
	return exp(-1*r*(T-t))* normalDistribution_builtin(d2(S,X,r,D0,sigma,T,t));
}

double Analyticsolution(double S, double X1,double X2, double r, double D0, double sigma, double T, double t)
{
	return P(S,X1,r,D0,sigma,T,t) - C(S,X2,r,D0,sigma,T,t) - 2 * X2 * BC(S,X2,r,D0,sigma,T,t);
}


int main()

{
     cout << "\nThe value of portfolio as time t=0, S0=X1 is: " << Analyticsolution(60000, 60000, 105000, 0.03, 0.03, 0.27, 1.75, 0) << endl;
     cout << "\nThe value of portfolio as time t=0, S0=X2 is: " << Analyticsolution(105000, 60000, 105000, 0.03, 0.03, 0.27, 1.75, 0) << endl;
   
     

}