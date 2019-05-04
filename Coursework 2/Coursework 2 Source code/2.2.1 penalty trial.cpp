#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <conio.h>
#include <chrono>
using namespace std;

//Introduce Thomas Solver
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

//American option embedded with call option
double ConvertibleBond_AO_CN(double Cp, double C, double alpha, double beta, double mu, double kappa, double R, double S0, double r, double X, double sigma, double T,
	double F, int iMax, int jMax, double SMax, double penalty, double tol, int iterMax)

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


		a[jMax] = //0 ;
			-0.5 * 1 / dS;
		b[jMax] = //1;
			0.5 * 1 / dS;
		c[jMax] = 0;
		d[jMax] = //S[jMax] * R * exp(-(kappa + r)*(T -i *dt)) + C/(alpha + r) * (exp(-alpha * t) - exp(-alpha +r*(i*dt - T))) + X*R * ( exp(r*(T -i*dt)) - exp(-(kappa + r)*(T- i*dt)));
			-0.5 * 1 / dS * vOld[jMax] + 0.5 * 1 / dS * vOld[jMax - 1] + R; //*exp(-(kappa + r)*(T - i*dt));


		int iter;

		for (iter = 0; iter<iterMax; iter++)
		{
			// create new vectors containing the FD approximations
			std::vector<double> aHat(a), bHat(b), cHat(c), dHat(d);
			// check on whether to apply the penalty
			for (int j = 0; j <= jMax; j++)
			{

				// if current value suggests apply penalty, adjust matrix equations
				if (vNew[j] < R*S[j])
				{


					if (i * dt < 2.4712)

					{
						bHat[j] = b[j] - penalty; dHat[j] = d[j] - penalty* max(Cp, R*S[j]);
					}
					else
					{
						bHat[j] = b[j] - penalty; dHat[j] = d[j] - penalty* R*S[j];
					}



				}

				else
				{

					bHat[j] = b[j]; dHat[j] = d[j];
				}


			}
			// now solve *Hat matrix problem with thomas
			thomasSolve(aHat, bHat, cHat, dHat);
			// dHat now contains next guess at solution v_j^{q+1}
			// check for differences between dHat and vNew
			double error = 0.;
			for (int j = 0; j <= jMax; j++)
			{
				error += (dHat[j] - vNew[j])*(dHat[j] - vNew[j]);
			}
			// update vNew
			vNew = dHat;
			// make an exit condition when solution is converged
			if (error<tol*tol)
			{
				// cout << " # solved after "<< penaltyIt << " iterations\n\n" << endl;
				break;
			}
		}
		if (iter >= iterMax)
		{
			std::cout << " Error NOT converging within required iterations\n";
			std::cout.flush();
			throw;
		}
		// set old=new
		vOld = vNew;
	}
	int jStar = S0 / dS;
	double sum = 0.;
	sum += (S0 - S[jStar]) / dS * vNew[jStar + 1];
	sum += (S[jStar + 1] - S0) / dS * vNew[jStar];
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
	iMax = 400; jMax = 40; SMax = 6 * F;

	for (S0 = 10; S0 <= 105; S0 += 5)
	{

		cout << ConvertibleBond_AO_CN(Cp, C, alpha, beta, mu, kappa, R, S0, r, X, sigma, T, F, iMax, jMax, SMax, 1e8, 1.e-8, 10000) << endl;

	}



}








