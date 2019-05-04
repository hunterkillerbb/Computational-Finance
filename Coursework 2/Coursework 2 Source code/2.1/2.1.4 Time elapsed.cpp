#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <conio.h>
#include <chrono>
using namespace std;

// A generic lagrange interpolation function
double lagrangeInterpolation(const vector<double>& y, const vector<double>& x, double x0, unsigned int n)
{
	if (x.size()<n)return lagrangeInterpolation(y, x, x0, x.size());
	if (n == 0)throw;
	int nHalf = n / 2;
	int jStar;
	double dx = x[1] - x[0];
	if (n % 2 == 0)
		jStar = int((x0 - x[0]) / dx) - (nHalf - 1);
	else
		jStar = int((x0 - x[0]) / dx + 0.5) - (nHalf);
	jStar = std::max(0, jStar);
	jStar = std::min(int(x.size() - n), jStar);
	if (n == 1)return y[jStar];
	double temp = 0.;
	for (unsigned int i = jStar; i<jStar + n; i++){
		double  int_temp;
		int_temp = y[i];
		for (unsigned int j = jStar; j<jStar + n; j++){
			if (j == i){ continue; }
			int_temp *= (x0 - x[j]) / (x[i] - x[j]);
		}
		temp += int_temp;
	}
	// end of interpolate
	return temp;
}

void thomasSolve(const std::vector<double> &a, const std::vector<double> &b_, const std::vector<double> &c, std::vector<double> &rhs)
{
	int n = a.size();
	std::vector<double> b(n);
	// initial first value of b
	b[0] = b_[0];
	for (int j = 1; j<n; j++)
	{
		b[j] = b_[j] - c[j - 1] * a[j] / b[j - 1];
		rhs[j] = rhs[j] - rhs[j - 1] * a[j] / b[j - 1];
	}
	// calculate solution
	rhs[n - 1] = rhs[n - 1] / b[n - 1];
	for (int j = n - 2; j >= 0; j--)
		rhs[j] = (rhs[j] - c[j] * rhs[j + 1]) / b[j];
}



double ConvertibleBond_CN(double C, double alpha, double beta, double mu, double kappa, double R, double S0, double r, double X, double sigma, double T,
	double F, int iMax, int jMax, double SMax, double omega, double tol, double iterMax)


{
	double dS = SMax / jMax;
	double dt = T / iMax;
	// create storage for the stock price and option price (old and new)
	std::vector<double> S(jMax + 1), vOld(jMax + 1), vNew(jMax + 1);
	// setup and initialise the stock price 
	for (int j = 0; j <= jMax; j++)
	{
		S[j] = j*dS;
	}
	// setup and initialize the final conditions on option price
	for (int j = 0; j <= jMax; j++)
	{
		vOld[j] = std::max(F, R*S[j]);
		vNew[j] = std::max(F, R*S[j]);
	}

	// run through timesteps
	for (int i = iMax - 1; i >= 0; i--)

	{

		// declare vectors for matrix equations
		std::vector<double> a(jMax + 1), b(jMax + 1), c(jMax + 1), d(jMax + 1);

		// set up matrix equations,define the initial conditions;


		a[0] = 0;
		b[0] = -1. / dt - 0.5 * r - 0.5 * 1 / dS * kappa *  (1 + mu)*X*exp(mu*(i + 0.5)*dt);
		c[0] = 0.5 * 1 / dS * kappa * (1 + mu)*X*exp(mu*(i + 0.5)*dt);
		d[0] = -C*exp(-alpha * (i + 0.5) *dt) + (-b[0] - 2 / dt) * vOld[0] - c[0] * vOld[1];


		for (int j = 1; j<jMax; j++)
		{


			a[j] = 0.25 * sigma * sigma*pow(j, 2 * beta) * pow(dS, 2 * beta - 2) - 0.25 *kappa * 1 / dS *  (1 + mu)*X*exp(mu*(i + 0.5)*dt) + 0.25 * kappa * j;

			b[j] = -1. / dt - 0.5* sigma*sigma*pow(j, 2 * beta) * pow(dS, 2 * beta - 2) - 0.5 * r;


			c[j] = 0.25 * sigma * sigma*pow(j, 2 * beta) * pow(dS, 2 * beta - 2) + 0.25 *kappa * 1 / dS *  (1 + mu)*X*exp(mu*(i + 0.5)*dt) - 0.25 * kappa * j;



			d[j] = -C*exp(-alpha * (i + 0.5) *dt) - a[j] * vOld[j - 1] + (-b[j] - 2. / dt)*vOld[j] - c[j] * vOld[j + 1];


		}


		a[jMax] = -0.5 * 1 / dS;
		b[jMax] = 0.5 * 1 / dS;
		c[jMax] = 0;
		d[jMax] = -0.5 * 1 / dS * vOld[jMax] + 0.5 * 1 / dS * vOld[jMax - 1] + R*exp(-(kappa + r)*(T - (i + 0.5)*dt));


		// solve with thomas solver + interpolation 
		vNew = d;
		thomasSolve(a, b, c, vNew);

		vOld = vNew;
	}
	return lagrangeInterpolation(vNew, S, S0, 4);
}





int main()
{
	// declare parameters
	double C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F;
	// declare grid paramaters 
	int iMax, jMax;
	double SMax;

	// initialise parameters
	C = 1.52; alpha = 0.02; beta = 0.198; mu = 0.0218; kappa = 0.05; R = 4; S0 = 44.12; r = 0.0169; X = 44.12; sigma = 4.98; T = 5;  F = 180;
	// initialise grid paramaters 
	iMax = 400; jMax = 40; SMax = 5 * F;


	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Press any key to run the timer..." << endl;
	_getch();


	double valueOld = 0.;
	iMax = 10000;

	cout << "\n1 on the go..." << endl;

	for (int n = 10; n <= 5000; n *= 2)
	{
		int jMax = n;

		auto start = std::chrono::steady_clock::now();
		double value = ConvertibleBond_CN(C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, iMax, jMax, SMax, 1, 1.e-8, 10000);
		valueOld = value;

		auto finish = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);

		cout << "The time process 1 trial " << n << " used is: " << elapsed.count() << endl;
	}





	valueOld = 0.;
	jMax = 10000;

	cout << "\n2 on the go..." << endl;

	for (int n = 10; n <= 5000; n *= 2)
	{
		int iMax = n;

		auto start = std::chrono::steady_clock::now();
		double value = ConvertibleBond_CN(C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, iMax, jMax, SMax, 1, 1.e-8, 10000);
		valueOld = value;

		auto finish = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
		cout << "The time process 2 trial " << n << " used is: " << elapsed.count() << endl;
	}





	valueOld = 0.;

	cout << "\n3 on the go..." << endl;

	for (int n = 10; n <= 5000; n *= 2)
	{
		int iMax = n;
		int jMax = n;

		auto start = std::chrono::steady_clock::now();
		double value = ConvertibleBond_CN(C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, iMax, jMax, SMax, 1, 1.e-8, 10000);
		valueOld = value;

		auto finish = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
		cout << "The time process 3 trial " << n << " used is: " << elapsed.count() << endl;
	}





}