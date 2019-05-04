#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <conio.h>

using namespace std;


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
				if ((i+0.5) * dt < 2.4712)
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
				if ((i+0.5) * dt < 2.4712)
				{
					y = min(y, max(Cp, R*S[j]));
				}
				error += fabs(vNew[j] - y);
				vNew[j] = y;
			}

			{
				double y = (d[jMax] - a[jMax] * vNew[jMax - 1]) / b[jMax];
				y = max(R*S[jMax], vNew[jMax] + omega*(y - vNew[jMax]));
				if ((i+0.5) * dt < 2.4712)
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
	int jStar = S0 / dS;
	double sum = 0.;
	sum += (S0 - S[jStar]) / dS *vNew[jStar + 1];
	sum += (S[jStar + 1] - S0) / dS *vNew[jStar];
	return sum;
}



int main()
{
	// declare parameters
	double C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, Cp;
	// declare grid paramaters 
	int iMax, jMax;
	double SMax;

	// initialise parameters
	Cp = 220; C = 1.52; alpha = 0.02; beta = 0.198; mu = 0.0218; kappa = 0.05; R = 4; S0 = 44.12; r = 0.0169; X = 44.12; sigma = 4.98; T = 5;  F = 180;
	// initialise grid paramaters 
	iMax = 400; jMax = 40; SMax = 5 * F;

	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Press any key to check American option ..." << endl;
	_getch();

	ofstream output;

	output.open("P:/American 1.csv");

	if (!output.is_open())
	{
		cout << "File not opened \n";
		throw; // stop the program here
	}

	for (S0 = 0; S0 <= 105; S0 += 5)
	{
		cout << ConvertibleBond_CN_AO(Cp, C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, iMax, jMax, SMax, 1, 1.e-8, 10000) << endl;
		output << ConvertibleBond_CN_AO(Cp, C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, iMax, jMax, SMax, 1, 1.e-8, 10000) << endl;

	}


	cout << "\nFile write successful \n";
	output.close();


	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Press any key to check American option ..." << endl;
	_getch();



	output.open("P:/American 2.csv");

	if (!output.is_open())
	{
		cout << "File not opened \n";
		throw; // stop the program here
	}

	for (C = 0.76; C < 2.29; C += 0.76)
	{
		for (S0 = 0; S0 <= 105; S0 += 5)

		{
			cout << C << "  ,  " << ConvertibleBond_CN_AO(Cp, C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, iMax, jMax, SMax, 1, 1.e-8, 10000) << endl;
			output << C << "," << ConvertibleBond_CN_AO(Cp, C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, iMax, jMax, SMax, 1, 1.e-8, 10000) << endl;

		}

	}
	cout << "\nFile write successful \n";
	output.close();


}