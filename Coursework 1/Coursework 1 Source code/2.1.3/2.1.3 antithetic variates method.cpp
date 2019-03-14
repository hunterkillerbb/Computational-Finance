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
		double ST_plus = S * exp((r - D0 - 0.5*sigma*sigma)*T + phi*sigma*sqrt(T));
		// add in payoff
		sum = sum + payoff(ST_plus, X1, X2);
		double ST_minus = S * exp((r - D0 - 0.5*sigma*sigma)*T - phi*sigma*sqrt(T));
		sum = sum + payoff(ST_minus,X1,X2);
	}
	// return discounted value
	return (sum / (2*N) )*exp(-r*T);
}

int main()
{
	
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "\nPress any key to generate the confidence intervals with antithetic variates method..." << endl;
	_getch();


	ofstream output1;

	output1.open("P:/Data_CI2.csv");

	if (!output1.is_open())
	{
		cout << "File not opened \n";
		throw; // stop the program here
	}

    auto start = std::chrono::steady_clock::now(); //start counting time

	for (int N = 1000; N <= 100000; N += 1000)
	{
		int M = 100;
		
		vector<double> samples(M);
		


		cout << "\nRun results with M=" << M << " samples from V_N, where N=" << N << "." << endl;
		

		// run some calculations
		for (int i = 0; i<M; i++)
		{
			samples[i] = monteCarlo(60000, 60000, 105000, 0.03, 0.03, 0.27, 1.75, N);
		}
		// estimate the mean from the sample
		double sum = 0.;
		for (int i = 0; i<M; i++)
		{
			sum += samples[i];
		}
		double mean = sum / M;
		cout << " mean = " << mean << endl;

		// estimate the variance from the sample
		double sumvar = 0.;
		for (int i = 0; i<M; i++)
		{
			sumvar += (samples[i] - mean)*(samples[i] - mean);
		}
		double variance = sumvar / (M - 1);
		cout << " variance = " << variance << endl;
		// get the standard deviation of the sample mean
		cout << " variance of the sample mean = " << variance / M << endl;
		double sd = sqrt(variance / M);
		cout << " 95% confident result is in [" << mean - 2.*sd << "," << mean + 2.*sd << "] with " << N*M << " total paths.\n" << endl;

		output1 << N << "," << mean - 2 * sd << "," << mean + 2 * sd << endl;
	}
     
    auto finish = std::chrono::steady_clock::now();// convert into real time in seconds
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
  
    cout << " Total time elapsed for the process is: "<<  elapsed.count() << endl;




	cout << "File write successful \n";
	output1.close();



	return 0;

}

