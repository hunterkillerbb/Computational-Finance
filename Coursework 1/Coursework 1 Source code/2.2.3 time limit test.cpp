#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <conio.h>
#include <time.h>
using namespace std;

double barriercallOption(double S, double  sigma, double  r, double D0, double  T, double  X, double B, int K, int N)
{
	static mt19937 rng;
	normal_distribution<> ND(0., 1.);
	double sum = 0.;
	double total = 0;
	for (int n = 0; n < N; n++)
	{
		// now create a path
		double dt = T / K;
		vector<double> stockPath(K + 1);
		vector<double> stockPath1(K + 1);
		vector<double> stockPath2(K + 1);
		stockPath[0] = S;
		stockPath1[0] = S;
		stockPath2[0] = S;


		for (int i = 1; i <= K; i++)
		{
			double phi = ND(rng);
			stockPath[i] = stockPath[i - 1] * exp((r - D0 - 0.5*sigma*sigma)*dt + phi*sigma*sqrt(dt));
			stockPath1[i] = stockPath1[i - 1] * exp((r - D0 - 0.5*sigma*sigma)*dt - phi*sigma*sqrt(dt));
			stockPath2[i] = (stockPath[i] + stockPath1[i])*0.5;
		}
		// and calculate A
		double A = S;
		for (int i = 1; i <= K; i++)
		{
			A = A + (stockPath2[i] - stockPath2[i - 1]);

			// add in the payoff to the sum

			if (A > B)

			{
				sum = 0;
				break;
			}

			sum = max(A - X, 0.);

		}
		total += sum;
	}
	return total / N*exp(-r*T);


}





int main()
{

	double S0 = 66500, sigma = 0.41, r = 0.06, T = 0.75, X = 66500, B = 113200, D0 = 0.06;
	int K = 35;
	int N = 1000;
	double sum = 0;
	int n;

	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Press any key to generate as many option prices as possible..." << endl;
	_getch();

	time_t start;
	start = time(NULL);

	for (int i = 0; i <= 1000000000; i++)
	{

		if ((time(NULL) - start) > 10)
			break;
		sum += barriercallOption(S0, sigma, r, D0, T, X, B, K, N);
		n = i;


	}

	cout << "\nThe average barrier option prices is: " << sum / n <<" with "<< n << " times simulations"<< endl;

	return 0;












}