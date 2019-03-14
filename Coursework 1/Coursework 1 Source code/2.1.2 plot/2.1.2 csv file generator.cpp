#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
#include <vector>
#include <fstream>
#include <conio.h>
using namespace std;

double putoption(double S, double X)
{
	return max(X - S, 0.);
}

double calloption(double S, double X)

{
	return max(S - X, 0.);
}

double binarycalloption(double S, double X)
{
	if (S>X)
	{
		return 1;
	}

	else
	{
		return 0;
	}
}


double payoff(double S, double X1, double X2)
{
	return putoption(S, X1) - calloption(S, X2) - 2 * X2 * binarycalloption(S, X2);
}

//the option value at time t = 0
double monteCarlo(double S, double X1, double X2, double r, double D0, double sigma, double T, int N)
{
	// declare the random number generator
	static mt19937 rng;
	// declare the distribution
	normal_distribution<> ND(0., 1.);
	ND(rng);
	// initialise sum
	double sum = 0.;
	for (int i = 0; i<N; i++)
	{
		double phi = ND(rng);
		// calculate stock price at T
		double ST = S * exp((r - D0 - 0.5*sigma*sigma)*T + phi*sigma*sqrt(T));
		// add in payoff
		sum = sum + payoff(ST, X1, X2);
	}
	// return discounted value
	return sum / N*exp(-r*T);
}

int main()
{
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "\nPress any key to output the CSV file..." << endl;
	_getch();


	ofstream output;

	output.open("P:/Data_value.csv");

	if (!output.is_open())
	{
		cout << "File not opened \n";
		throw; // stop the program here
	}

	for (int N = 1000; N <= 100000; N += 1000)
	{
		output << N << " , " << monteCarlo(60000, 60000, 105000, 0.03, 0.03, 0.27, 1.75, N) << " , " << monteCarlo(105000, 60000, 105000, 0.03, 0.03, 0.27, 1.75, N) << endl;
	}
	// file write successful then close file
	cout << "File write successful \n";
	output.close();

	return 0;
}