/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Covariance Filter.
 */
#include "bayesFlt.hpp"
#include "matSup.hpp"
#include "covFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


Covariance_scheme::Covariance_scheme (size_t x_size, size_t z_initialsize) :
	Extended_filter(x_size),
	S(Empty), SI(Empty), W(Empty)
/*
 * Initialise filter and set the size of things we know about
 */
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

Covariance_scheme& Covariance_scheme::operator= (const Covariance_scheme& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Extended_filter::operator=(a);
	return *this;
}


void Covariance_scheme::init ()
{
						// Postconditions
	if (!isPSD (X))
		error (Numeric_exception("Initial X not PSD"));
}

void Covariance_scheme::update ()
{
	// Nothing to do, implicit in observation
}

Bayes_base::Float
 Covariance_scheme::predict (Linrz_predict_model& f)
{
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict state covariance
	RowMatrix temp(f.Fx.size1(), X.size2());
	X = prod_SPD(f.Fx,X, temp) + prod_SPD(f.G, f.q);

	assert_isPSD (X);
	return 1;
}

void Covariance_scheme::observe_size (size_t z_size)
/*
 * Optimised dynamic observation sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		S.resize(z_size,z_size);
		SI.resize(z_size,z_size);
		W.resize(x.size(),z_size);
	}
}

Bayes_base::Float
 Covariance_scheme::observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s)
/* correlated innovation observe
 */
{
						// Size consistency, z to model
	if (s.size() != h.Z.size1())
		error (Logic_exception("observation and model size inconsistent"));
	observe_size (s.size());// Dynamic sizing

						// Innovation covariance
	Matrix temp_XZ (prod(X, trans(h.Hx)));
	S.assign (prod(h.Hx, temp_XZ) + h.Z);

						// Inverse innovation covariance
	Float rcond = UdUinversePD (SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

						// Kalman gain, X*Hx'*SI
	W.assign (prod(temp_XZ, SI));

						// State update
	x.plus_assign (prod(W, s));
	X.minus_assign (prod_SPD(W, S, temp_XZ));

	assert_isPSD (X);
	return rcond;
}


Bayes_base::Float
 Covariance_scheme::observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s)
/* uncorrelated innovation observe
 */
{
						// Size consistency, z to model
	if (s.size() != h.Zv.size())
		error (Logic_exception("observation and model size inconsistent"));
	observe_size (s.size());// Dynamic sizing

						// Innovation covariance
	Matrix temp_XZ (prod(X, trans(h.Hx)));
	S.assign (prod(h.Hx, temp_XZ));
	for (size_t i = 0; i < h.Zv.size(); ++i)
		S(i,i) += h.Zv[i];

						// Inverse innovation covariance
	Float rcond = UdUinversePD (SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

						// Kalman gain, X*Hx'*SI
	W.assign (prod(temp_XZ, SI));

						// State update
	x.plus_assign (prod(W, s));
	X.minus_assign (prod_SPD(W, S, temp_XZ));

	assert_isPSD (X);
	return rcond;
}

}//namespace
