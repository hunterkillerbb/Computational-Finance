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



double ConvertibleBond_CN_AO(double Cp, double C, double alpha, double beta, double mu, double kappa, double R, double S0, double r, double X, double sigma, double T,
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


		// As S->infty, define the final conditions 

		a[jMax] = -0.5 * 1 / dS;
		b[jMax] = 0.5 * 1 / dS;
		c[jMax] = 0;
		d[jMax] = -0.5 * 1 / dS * vOld[jMax] + 0.5 * 1 / dS * vOld[jMax - 1] + R;


		// solve with PSOR method

		int sor;
		for (sor = 0; sor < iterMax; sor++)
		{
			double error = 0.;
			// implement sor in here
			{
				double y = (d[0] - c[0] * vNew[1]) / b[0];
				// track changes in value
				y = max(0., vNew[0] + omega*(y - vNew[0]));
				if ((i+0.5)* dt <= 2.4712)
				{
					y = min(y, max(Cp, 0.));
				}
				error += fabs(vNew[0] - y);
				vNew[0] = y;
			}
			for (int j = 1; j < jMax; j++)
			{
				double y = (d[j] - a[j] * vNew[j - 1] - c[j] * vNew[j + 1]) / b[j];
				y = max(R*S[j], vNew[j] + omega*(y - vNew[j]));
				if ((i+0.5) * dt <= 2.4712)
				{
					y = min(y, max(Cp, R*S[j]));
				}
				error += fabs(vNew[j] - y);
				vNew[j] = y;
			}

			{
				double y = (d[jMax] - a[jMax] * vNew[jMax - 1]) / b[jMax];
				y = max(R*S[jMax], vNew[jMax] + omega*(y - vNew[jMax]));
				if ((i+0.5) * dt <= 2.4712)
				{
					y = min(y, max(Cp, R*S[jMax]));
				}
				error += fabs(vNew[jMax] - y);
				vNew[jMax] = y;
			}
			if (error<tol)
				break;
		}
		if (sor >= iterMax)
		{
			std::cout << " Error NOT converging within required iterations\n";
			std::cout.flush();
			throw;
		}

		//set vOld = vNew
		vOld = vNew;
	}


	return lagrangeInterpolation(vNew, S, S0, 4);

}




int main()
{
	// declare parameters
	double Cp, C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F;
	// declare grid paramaters 
	int iMax, jMax;
	double SMax;

	// initialise parameters
	Cp = 220;  C = 1.52; alpha = 0.02; beta = 0.198; mu = 0.0218; kappa = 0.05; R = 4; S0 = 44.12; r = 0.0169; X = 44.12; sigma = 4.98; T = 5;  F = 180;
	// initialise grid paramaters 
	iMax = 400; jMax = 40; SMax = (1103 / 225) * F;


	ofstream output;


	double  valueOld1 = 0.;double valueOld2 = 0; 

	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Press any key to iterate our grid parameters..." << endl;
	_getch();



	output.open("P:/2.2.3 best estimate.csv");
	if (!output.is_open())
	{
		cout << "File not opened \n";
		throw; // stop the program here
	}


    auto start = std::chrono::steady_clock::now();
    
	for (int n = 10; n <= 1280; n*= 2)
	{

		int iMax = n;
		int jMax = n;
		
		double value = ConvertibleBond_CN_AO(Cp, C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, iMax, jMax, SMax, 1, 1.e-8, 10000);


		auto finish = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);

		if (elapsed.count() > 1.)
		{
			break;
		}
		cout <<  n << " , "<< elapsed.count() << "  ,  " << value << " , " << (4 * value - valueOld2) / 3 << " , "<<(16 * value - 8*valueOld2 + valueOld1) / 9 <<  endl;
		output <<  n << " , "<< elapsed.count() << "  ,  " << value << " , " << (4 * value - valueOld2) / 3 << " , "<<(16 * value - 8*valueOld2 + valueOld1) / 9 <<  endl;
        valueOld1 = valueOld2;
		valueOld2= value ;


	}



	cout << "\nFile write successful \n";
	output.close();




}