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
#include "matSup.h"
#include "covFlt.h"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


Covariance_filter::Covariance_filter (size_t x_size, size_t z_initialsize) :
	Extended_filter(x_size),
	S(Empty), SI(Empty)
/*
 * Initialise filter and set the size of things we know about
 */
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

Covariance_filter& Covariance_filter::operator= (const Covariance_filter& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Extended_filter::operator=(a);
	return *this;
}


void Covariance_filter::init ()
{
						// Preconditions
	if (!isPSD (X))
		filter_error ("Xi not PSD");
}

void Covariance_filter::update ()
{
	// Nothing to do, implicit in observation
}

Bayes_base::Float
 Covariance_filter::predict (Linrz_predict_model& f)
{
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict state covariance
	X = mult_SPD(f.Fx, X) + mult_SPD(f.G, f.q);

	// TODO: Remove check once proved
	assert_isPSD (X);
	return 1.;
}

void Covariance_filter::observe_size (size_t z_size)
/*
 * Optimised dynamic observation sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		S.resize(z_size,z_size);
		SI.resize(z_size,z_size);
	}
}

Bayes_base::Float
 Covariance_filter::observe_innovation (Linrz_correlated_observe_model& h, const Vec& s)
/* correlated innovation observe
 */
{
						// Size consistency, z to model
	if (s.size() != h.Z.size1())
		filter_error("observation and model size inconsistent");
	observe_size (s.size());// Dynamic sizing

						// Innovation covariance
	S = mult_SPD(h.Hx, X) + h.Z;

						// Inverse innovation covariance
	Float rcond = UdUinversePD (SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

						// Kalman gain
	Matrix W = prod(prod(X,trans(h.Hx)), SI);

						// State update
	x += prod(W, s);
	X -= mult_SPD(W, S);
	// TODO: Remove check once proved
	assert_isPSD (X);
	return rcond;
}


Bayes_base::Float
 Covariance_filter::observe_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s)
/* uncorrelated innovation observe
 */
{
						// Size consistency, z to model
	if (s.size() != h.Zv.size())
		filter_error("observation and model size inconsistent");
	observe_size (s.size());// Dynamic sizing

						// Innovation covariance
	S = mult_SPD(h.Hx, X);
	for (size_t i = 0; i < h.Zv.size(); ++i)
		S(i,i) += h.Zv[i];

						// Inverse innovation covariance
	Float rcond = UdUinversePD (SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

						// Kalman gain
	Matrix W = prod(prod(X,trans(h.Hx)), SI);

						// State update
	x += prod(W, s);
	X -= mult_SPD(W, S);
	// TODO: Remove check once proved
	assert_isPSD (X);
	return rcond;
}

}//namespace
