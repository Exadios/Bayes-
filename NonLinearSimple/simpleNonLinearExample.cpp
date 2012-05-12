/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens
 * See accompanying Bayes++.htm for terms and conditions of use.
 *
 * $Id$
 */

/*
 * Example of using Bayesian Filter Class to solve a simple problem.
 *  A Non-linear filter with one state and constant noises
 */

// Use the Unscented filter scheme from Bayes++
#include "BayesFilter/unsFlt.hpp"
#include <iostream>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace Bayesian_filter;
using namespace Bayesian_filter_matrix;


/*
 * Simple Prediction model
 */
class Simple_predict : public Linear_predict_model
{
public:
	Simple_predict() : Linear_predict_model(1,1)
	// Construct a constant model
	{
				// Stationary Prediction model (Identity)
		Fx(0,0) = 1.;
				// Constant Noise model with a large variance
		q[0] = 2.;
		G(0,0) = 1.;
	}
};

/*
 * Nonlinear Observation model
 * Observation is distance from origin squared
 */
class Nonlinear_observe : public Uncorrelated_additive_observe_model
{
public:
	Nonlinear_observe () : Uncorrelated_additive_observe_model(1), hx(1)
	// Construct a noise model
	{
				// Constant Observation Noise model with variance of one
		Zv[0] = 1.;
	}
	const FM::Vec& h(const FM::Vec& x) const
	{
		hx[0] = x[0]*x[0];
		return hx;
	}
	private:
		mutable FM::Vec hx;
};


/*
 * Filter a simple example
 */
int main()
{
	// Global setup for test output
	cout.flags(ios::fixed); cout.precision(4);

	// Construct simple Prediction and Observation models
	Simple_predict my_predict;
	Nonlinear_observe my_observe;

	// Use an 'Unscented' filter scheme with one state
	Unscented_scheme my_filter(1);

	// Setup the initial state and covariance
	Vec x_init (1); SymMatrix X_init (1, 1);
	x_init[0] = 10.;		// Start at 10 with no uncertainty
	X_init(0,0) = 0.;
	my_filter.init_kalman (x_init, X_init);

	cout << "Initial  " << my_filter.x << my_filter.X << endl;

	// Predict the filter forward
	my_filter.predict (my_predict);
	my_filter.update();		// Update the filter, so state and covariance are available

	cout << "Predict  " << my_filter.x << my_filter.X << endl;

	// Make an observation
	Vec z(1);
	z[0] = 11.;				// Observe that we should be at 11
	my_filter.observe (my_observe, z);
	my_filter.update();		// Update the filter to state and covariance are available

	cout << "Filtered " << my_filter.x << my_filter.X << endl;
	return 0;
}
