#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <conio.h>
using namespace std;

double barriercallOption(double S, double  sigma, double  r, double D0, double  T, double  X, double B, int K, int N)
{
	static mt19937 rng;
	normal_distribution<> ND(0., 1.);
	double sum = 0.;
	double total = 0;
	int barriernum = 0;
	for (int n = 0; n < N; n++)
	{
		// now create a path
		double dt = T / K;
		vector<double> stockPath(K + 1);
		stockPath[0] = S;
		for (int i = 1; i <= K; i++)
		{
			double phi = ND(rng);
			stockPath[i] = stockPath[i - 1] * exp((r - D0 - 0.5*sigma*sigma)*dt + phi*sigma*sqrt(dt));
		}
		// and calculate A
		double A = S;
		for (int i = 1; i <= K; i++)
		{
			A = A + (stockPath[i] - stockPath[i - 1]);

			// add in the payoff to the sum

			if (A > B)

			{
				barriernum += 1;
				cout << "\nHit the barrier " << barriernum << " time(s)..." << endl;
				sum = 0;
				break;
			}

			sum = max(A - X, 0.);

		}

		total += sum;
		cout << "\nIteration for path number: " << N << " is finished..." << endl;
	}
	return total / N*exp(-r*T);


}

int main()
{

	double S0 = 66500, sigma = 0.41, r = 0.06, T = 0.75, X = 66500, B = 113200, D0 = 0.06;
	int K = 35;
	int N = 100000;
	cout << "\nThe up-and-out barrier Call Option value is: " << barriercallOption(S0, sigma, r, D0, T, X, B, K, N) << endl;

	return 0;
}