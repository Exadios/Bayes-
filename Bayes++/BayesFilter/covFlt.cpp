/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Covariance Filter.
 */
#include "covFlt.hpp"
#include "matSup.hpp"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


Covariance_scheme::Covariance_scheme (size_t x_size) :
	Kalman_state_filter(x_size),
	tempX(x_size,x_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
}

Covariance_scheme& Covariance_scheme::operator= (const Covariance_scheme& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Kalman_state_filter::operator=(a);
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
/* Standard Linrz predict
 */
{
	x = f.f(x);		// Extended Kalman state predict is f(x) directly
						// Predict state covariance
	noalias(X) = prod_SPD(f.Fx,X, tempX);
	noalias(X) += prod_SPD(f.G, f.q, tempX);

	return 1;
}

Bayes_base::Float
 Covariance_scheme::predict (Gaussian_predict_model& f)
/* Specialised 'stationary' predict, only addative noise
 */
{
						// Predict state covariance, simply add in noise
	noalias(X) += prod_SPD(f.G, f.q, tempX);
  
	return 1;
}


Bayes_base::Float
 Covariance_scheme::observe_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s)
/* Extended_kalman_filter observe, unused byproduct
 */
{
	const size_t z_size = h.Hx.size1();
	Covariance_byproduct S(z_size, z_size);
	Kalman_gain_byproduct b(h.Hx.size2(), z_size);
	return eobserve_innovation (h, s, S, b);
}	

Bayes_base::Float
 Covariance_scheme::observe_innovation (Linrz_correlated_observe_model& h, const Vec& s)
/* Extended_kalman_filter observe, unused byproduct
 */
{
	const size_t z_size = h.Hx.size1();
	Covariance_byproduct S(z_size, z_size);
	Kalman_gain_byproduct b(h.Hx.size2(), z_size);
	return eobserve_innovation (h, s, S, b);
}	

Bayes_base::Float
 Covariance_scheme::eobserve_innovation (Linrz_correlated_observe_model& h, const Vec& s,
 				Covariance_byproduct& S, Kalman_gain_byproduct& b)
 /* Correlated innovation observe with explict byproduct
 */
{
						// Size consistency, z to model
	if (s.size() != h.Z.size1())
		error (Logic_exception("observation and model size inconsistent"));

						// Innovation covariance
	Matrix temp_XZ (prod(X, trans(h.Hx)));
	noalias(S) = prod(h.Hx, temp_XZ) + h.Z;

						// Inverse innovation covariance
	Float rcond = UdUinversePD (b.SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

						// Kalman gain, X*Hx'*SI
	noalias(b.W) = prod(temp_XZ, b.SI);

						// State update
	noalias(x) += prod(b.W, s);
	noalias(X) -= prod_SPD(b.W, S, temp_XZ);

	return rcond;
}


Bayes_base::Float
 Covariance_scheme::eobserve_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s,
 				Covariance_byproduct& S, Kalman_gain_byproduct& b)
/* Uncorrelated innovation observe with explict byproduct
 */
{
						// Size consistency, z to model
	if (s.size() != h.Zv.size())
		error (Logic_exception("observation and model size inconsistent"));

						// Innovation covariance
	Matrix temp_XZ (prod(X, trans(h.Hx)));
	noalias(S) = prod(h.Hx, temp_XZ);
	for (size_t i = 0; i < h.Zv.size(); ++i)
		S(i,i) += Float(h.Zv[i]);	// ISSUE mixed type proxy assignment

						// Inverse innovation covariance
	Float rcond = UdUinversePD (b.SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

						// Kalman gain, X*Hx'*SI
	noalias(b.W) = prod(temp_XZ, b.SI);

						// State update
	noalias(x) += prod(b.W, s);
	noalias(X) -= prod_SPD(b.W, S, temp_XZ);

	return rcond;
}

}//namespace
