/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 */

/*
 * Example of using Bayesian Filter Class to solve a simple problem.
 *
 * The example implements a simple quadratic observer.
 *  This trys to estimate the state of system while also trying to
 *  calibrate a simple linear model of the system which includes
 *  a scale factor and a bias.
 *  Estimating both the system state and a scale factor results in a
 *  quadtratic (product of two states and therefore non-linear) observation.
 *  The system model is a 1D brownian motion with a known pertubation.
 */

#include "BayesFilter/allFilters.hpp"
#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>

namespace
{
	namespace FM = Bayesian_filter_matrix;
	using namespace FM;

	// Choose Filtering Scheme to use
	typedef Bayesian_filter::Information_scheme FilterScheme;

	// Square 
	template <class scalar>
	inline scalar sqr(scalar x)
	{
		return x*x;
	}

	// Random numbers from Boost
	boost::mt19937 localRng;
	inline void randNormal(FM::Vec& v)
	{
		boost::normal_distribution<boost::mt19937> gen(localRng);
		std::generate (v.begin(), v.end(), gen);
	}
	inline void randNormal(FM::Vec& v, const Float mean, const Float sigma)
	{
		boost::normal_distribution<boost::mt19937> gen(localRng, mean, sigma);
		std::generate (v.begin(), v.end(), gen);
	}
	inline void randUniform01(FM::Vec& v)
	{
		boost::uniform_01<boost::mt19937> gen(localRng);
		std::generate (v.begin(), v.end(), gen);
	}


	// Constant Dimensions
	const unsigned NX = 3;			// Filter State dimension 	(SystemState, Scale, Bias)

	// Filter Parameters
	// Noise on observing system state
	const Float OBS_NOISE = 0.01;
	// Prediction Noise: no Prediction noise as pertubation is known
	const Float X_NOISE = 0.0;	// System State	
	const Float S_NOISE = 0.0;	// Scale
	const Float B_NOISE = 0.0;	// Bias
	// Filter's Initial state uncertainty: System state is unknown
	const Float i_X_NOISE = 1000.;
	const Float i_S_NOISE = 0.1;
	const Float i_B_NOISE = 0.1;

}//namespace


/*
 * Prediction model
 * Linear state predict model with additive control input
 */
class QCpredict : public Bayesian_filter::Linrz_predict_model
{
	Float motion;
	mutable FM::Vec fx;
public:
	QCpredict();
		;
	void predict(const FM::Vec& u)
	{
		motion = u[0];
	}
	const FM::Vec& f(const FM::Vec& x) const
	{
		// Constant scale and bias, system state pertubed by control input
		fx = x;
		fx[0] += motion;
		return fx;
	};
};

QCpredict::QCpredict() : Bayesian_filter::Linrz_predict_model(NX, NX), fx(NX)
{
	FM::identity (Fx);

	// Setup constant noise model: G is identity
	q[0] = sqr(X_NOISE);
	q[1] = sqr(S_NOISE);
	q[2] = sqr(B_NOISE);
	FM::identity (G);
}


/*
 * Quadratic observation model
 */
class QCobserve : public Bayesian_filter::Linrz_uncorrelated_observe_model
{
	mutable FM::Vec z_pred;
public:
	QCobserve ();
	const FM::Vec& h(const FM::Vec& x) const
	{	// Quadratic Observation model
		z_pred[0] = x[0] * x[1] + x[2];
		return z_pred;
	};
	void state (const FM::Vec& x)
	// Linearised model, Jacobian of h at x
	{
		Hx(0,0) = x[1];
		Hx(0,1) = x[0];
		Hx(0,2) = 1.;
	}
};

QCobserve::QCobserve () :
	Bayesian_filter::Linrz_uncorrelated_observe_model(NX,1), z_pred(1)
{
	// Observation Noise variance
	Zv[0] = OBS_NOISE*OBS_NOISE;
}


/*
 * Quadratic Filter
 *  Derived for FilterScheme class
 */
class QCfilter : public FilterScheme {
public:
	QCfilter (const Float SystemStateGuess);
private:
};

QCfilter::QCfilter (const Float SystemStateGuess)
	: FilterScheme (NX)
/*
 * Construct the filter
 *  Numerical alogrithm requires construction with fixed dimensions
 * Parameters:
 *  An initial guess for the system state
 */
{
	// Initialise filter state and covaince

	// Initial State
	x[0] = SystemStateGuess;
	x[1] = 1.;		// Assumed initial Scale
	x[2] = 0.;		// Assumed initial Bias
	X.clear();
	X(0,0) = sqr(i_X_NOISE);
	X(1,1) = sqr(i_S_NOISE);
	X(2,2) = sqr(i_B_NOISE);

	init ();
}


int main()
{
	// Global setup for test output
	std::cout.flags(std::ios::scientific); std::cout.precision(6);

	// Setup the test filters
	FM::Vec x_true (NX);

	// True State to be observed
	x_true[0] = 10.;	// System State
	x_true[1] = 1.0;	// Scale
	x_true[2] = 0.0;	// Bias

	std::cout << "Quadratic Calibration" << std::endl;
	std::cout << x_true << std::endl;


	// Construct Prediction and Observation model and Calibration filter
	// Give the filter an initial guess of the system state
	QCpredict linearPredict;
	QCobserve nonlinObserve;
	QCfilter obsAndCalib (x_true[0]);

	// Iterate the filter with test observations
	FM::Vec u(1), z_true(1), z(1);
	for (unsigned i = 0; i < 100; i++ )
	{
		// Predict true state using control input brownian
		randNormal (u);				// normally distributed
		x_true[0] += u[0];
		linearPredict.predict (u);

		// Predict filter with known pertubation
		obsAndCalib.predict (linearPredict);

		// True Observation: Quadratic observation model
		z_true[0] = x_true[0] * x_true[1] + x_true[2];

		// Observation with addative noise
		randNormal (z, z_true[0], OBS_NOISE);	// normally distributed mean z_true[0], stdDev OBS_NOISE.

		// Filter observation using model linearised at state estimate x
		nonlinObserve.state (obsAndCalib.x);
		obsAndCalib.observe (nonlinObserve, z);
	}

	// Update the filter to state and covariance are available
	obsAndCalib.update ();

	// Print everything: True, filter, covariance
	std::cout << x_true <<  std::endl;
	std::cout << obsAndCalib.x << std::endl;
	std::cout << obsAndCalib.X << std::endl;
	return 0;
}
