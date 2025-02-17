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
	
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Press any key to generate the first CSV file..." << endl;
	_getch();

	ofstream output;


	output.open("P:/Data_K&N.csv");

	if (!output.is_open())
	{
		cout << "File not opened \n";
		throw; // stop the program here
	}

	for (int K = 35; K <= 350; K += 35)
	{
		for (int N = 1000; N <= 20000; N += 1000)
		{
			output << K << "," << N << "," << barriercallOption(66500, 0.42, 0.06, 0.06, 0.75, 66500, 113200, K, N) << endl;
		}

	}

	cout << "File write successful \n";
	output.close();
    
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Press any key to generate the first CSV file..." << endl;
	_getch();

	ofstream output1;
		output.open("P:/Data_K&N 2.csv");

	if (!output1.is_open())
	{
		cout << "File not opened \n";
		throw; // stop the program here
	}

	for (int K = 35; K <= 35000; K *= 35)
	{
		for (int N = 100; N <= 100000 ; N *= 10)
		{
			output << K << "," << N << "," << barriercallOption(66500, 0.42, 0.06, 0.06, 0.75, 66500, 113200, K, N) << endl;
		}

	}

	cout << "File write successful \n";
	output1.close();
	
	
	return 0;

}