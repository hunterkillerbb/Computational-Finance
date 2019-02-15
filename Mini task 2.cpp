#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <conio.h>
using namespace std;

// Return a series of variables with parameters r, t, T, kappa, theta and sigma 
double m(double r, double t, double T, double kappa, double theta)

{
	return exp(-1. * kappa * (T - t)) * r + (1. - exp(-1. * kappa *(T - t)))*theta;
}

double n(double r, double t, double T, double kappa, double theta)
{
	return (T - t)*theta + (r - theta)* (1 - exp(-1.* kappa*(T - t))) / kappa;
}

double kSquared(double t, double T, double kappa, double sigma)
{
	return sigma*sigma / 2. / kappa / kappa / kappa
		*(4.*exp(-kappa*(T - t)) - exp(-2.*kappa*(T - t))
		+ 2.*kappa*(T - t) - 3.);
}

double q(double t, double T, double kappa, double sigma)
{
	return ((sigma*sigma) / (2 * kappa*kappa)) * (1 - exp(-1.*kappa* (T - t)))*(1 - exp(-1.*kappa* (T - t)));
}

double vSquared(double t, double T,
	double kappa, double sigma)
{
	return sigma*sigma / 2. / kappa*(1 - exp(-2.*kappa*(T - t)));
}

double f(double r, double t, double T, double kappa, double theta, double sigma)

{
	return m(r, t, T, kappa, theta) - q(t, T, kappa, sigma);
}


double P(double r, double t, double T, double kappa, double theta, double sigma)
{
	return exp(0.5 * kSquared(t, T, kappa, theta) *kSquared(t, T, kappa, theta) - n(r, t, T, kappa, theta));
}
double normalDistribution_builtin(double x)
{
	return 0.5*erfc(-x / sqrt(2.));
}
double N(double Xr, double r, double t, double T, double kappa, double theta, double sigma) 
// Standardize the current normal distribution and convert Xr to h
{

	return normalDistribution_builtin((Xr - f(r, t, T, kappa, theta, sigma)) / sqrt(vSquared(t, T, kappa, sigma)));
}


double V(double Xr, double r, double t, double T, double kappa, double theta, double sigma)
{
	return P(r, t, T, kappa, theta, sigma)* (1 - N(Xr, r, t, T, kappa, theta, sigma));
}

// Main body
int main()
{

	double t = 0., T = 9.;
	double Xr = 0.1, kappa = 0.3528, theta = 0.088, sigma = 0.0284, r = 0.0302;

	cout << "\n\nStep1--The value of the option at time t = 0 is equal to: " << V(Xr, r, t, T, kappa, theta, sigma) << endl;
	cout << "\npress any key to continue: " << endl;
	_getch();

// Generate 3 columns of data
	cout << "\n\nStep2--To output the bond and option prices for different interest rates, first we generate the value of all r: " << endl;

	for (double r = 0; r < 0.201; r += 0.002)
	{
		cout << r << endl;
	}

	cout << "\n\nThen, the values of P(r,t=0,T) are equal to:" << endl;
	cout << "\npress and key to continue: " << endl;
	_getch();

	for (double r = 0; r < 0.201; r += 0.002)
	{
		cout << P(r, t, T, kappa, theta, sigma) << endl;
	}

	cout << "\n\nFinally, the values of V(r,t=0) are equal to: " << endl;
	cout << "\npress and key to continue: " << endl;
	_getch();

	for (double r = 0.; r < 0.201; r += 0.002)
	{
		cout << V(Xr, r, t, T, kappa, theta, sigma) << endl;
	}

	cout << "\nOutput is complete\n\n" << endl;
	cout << "\npress and key to continue: " << endl;
	_getch();

// open up a file stream to write data
	ofstream output;
	
	output.open("P:/Data.csv");

	if (!output.is_open())
	{
		cout << "File not opened \n";
		throw; // stop the program here
	}
	
	for (double r = 0; r <0.201; r += 0.002)
	{
	output << r << " , " << P(r, t, T, kappa, theta, sigma) << " , " << V(Xr, r, t, T, kappa, theta, sigma) << endl;
	}
	// file write successful then close file
	cout << "File write successful \n";
	output.close();

	return 0;
}
