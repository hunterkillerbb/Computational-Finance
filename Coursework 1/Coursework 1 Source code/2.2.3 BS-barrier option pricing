#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

double d1(double S, double X, double r, double D0, double sigma, double T, double t)
{
	return (log(S / X) + (r - D0 + 0.5 *sigma*sigma)*(T - t)) / (sigma * sqrt(T - t));
}

double d2(double S, double X, double r, double D0, double sigma, double T, double t)
{
	return d1(S, X, r, D0, sigma, T, t) - sigma * sqrt(T - t);
}

double lambda(double r, double D0, double sigma)
{
	return (r - D0 + 0.5 *sigma*sigma) / sigma*sigma;
}

double x1(double S, double B, double r, double D0, double sigma, double T, double t)
{
	return (log(S / B) / sigma * sqrt(T - t)) + lambda(r, D0, sigma) * sigma * sqrt(T - t);
}

double y1(double S, double B, double r, double D0, double sigma, double T, double t)
{
	return (log(B / S) / sigma * sqrt(T - t)) + lambda(r, D0, sigma) * sigma * sqrt(T - t);
}

double y(double S, double X, double B, double r, double D0, double sigma, double T, double t)
{
	return (log((B*B) / (S*X)) / sigma * sqrt(T - t)) + lambda(r, D0, sigma) * sigma * sqrt(T - t);

}



double ND(double x)
{
	return 0.5*erfc(-x / sqrt(2.));
}

double N1(double S, double X, double B, double r, double D0, double sigma, double T, double t)
{
	return ND(d1(S, X, r, D0, sigma, T, t)) - ND(x1(S, B, r, D0, sigma, T, t)) + pow((B / S), 2 * lambda(r, D0, sigma))* (ND(-y(S, X, B, r, D0, sigma, T, t)) - ND(-y1(S, B, r, D0, sigma, T, t)));
}

double N2(double S, double X, double B, double r, double D0, double sigma, double T, double t)
{
	return -ND(d2(S, X, r, D0, sigma, T, t)) + ND(x1(S, B, r, D0, sigma, T, t) - sigma*sqrt(T - r)) - pow((B / S), 2 * lambda(r, D0, sigma) - 2)* (ND(-y(S, X, B, r, D0, sigma, T, t) + sigma * sqrt(T - t)) - ND(-y1(S, B, r, D0, sigma, T, t) + sigma*sqrt(T - t)));
}

double Analyticsolution(double S, double X, double B, double r, double D0, double sigma, double T, double t)
{
	return S*exp(-1 * D0*(T - t))* N1(S,X,B,r,D0,sigma,T,t)
		+ X*exp(-1 * r*(T - t))* N2( S,  X,  B,  r,D0,  sigma,  T,  t);
}


int main()

{
	cout << "\nThe barrier option ends up with the value: " << Analyticsolution(66500, 66500, 113200, 0.06, 0.06, 0.41, 0.75, 0) << endl;


}