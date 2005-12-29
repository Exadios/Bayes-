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


Covariance_bscheme::Covariance_bscheme (std::size_t x_size) :
	Kalman_state_filter(x_size),
	tempX(x_size,x_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
}

Covariance_scheme::Covariance_scheme (std::size_t x_size, std::size_t z_initialsize) :
	Kalman_state_filter(x_size),
	Covariance_bscheme(x_size),
	S(Empty), SI(Empty), W(Empty), XtHx(Empty)
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

Covariance_bscheme& Covariance_bscheme::operator= (const Covariance_bscheme& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Kalman_state_filter::operator=(a);
	return *this;
}

void Covariance_scheme::observe_size (std::size_t z_size)
/*
 * Optimised dynamic observe sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		S.resize(z_size,z_size, false);
		SI.resize(z_size,z_size, false);
		W.resize(x.size(),z_size, false);
	}
}


void Covariance_bscheme::init ()
{
						// Postconditions
	if (!isPSD (X))
		error (Numeric_exception("Initial X not PSD"));
}

void Covariance_bscheme::update ()
{
	// Nothing to do, implicit in observation
}

Bayes_base::Float
 Covariance_bscheme::predict (Linrz_predict_model& f)
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
 Covariance_bscheme::predict (Gaussian_predict_model& f)
/* Specialised 'stationary' predict, only addative noise
 */
{
						// Predict state covariance, simply add in noise
	noalias(X) += prod_SPD(f.G, f.q, tempX);
  
	return 1;
}


Bayes_base::Float
 Covariance_bscheme::byobserve_innovation (Linrz_correlated_observe_model& h, const Vec& s,
 				FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W, FM::Matrix& XtHx)
/* Correlated innovation observe with explict byproduct
*/
{
						// Size consistency, z to model
	if (s.size() != h.Z.size1())
		error (Logic_exception("observation and model size inconsistent"));

						// Innovation covariance
	noalias(XtHx) = prod(X, trans(h.Hx));
	noalias(S) = prod(h.Hx, XtHx) + h.Z;

						// Inverse innovation covariance
	Float rcond = UdUinversePD (SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

						// Kalman gain, X*Hx'*SI
	noalias(W) = prod(XtHx, SI);

						// State update
	noalias(x) += prod(W, s);
	noalias(X) -= prod_SPD(W, S, XtHx);

	return rcond;
}

Bayes_base::Float
 Covariance_bscheme::byobserve_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s,
 				FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W, FM::Matrix& XtHx)
/* Uncorrelated innovation observe with explict byproduct
 */
{
						// Size consistency, z to model
	if (s.size() != h.Zv.size())
		error (Logic_exception("observation and model size inconsistent"));

						// Innovation covariance
	noalias(XtHx) = prod(X, trans(h.Hx));
	noalias(S) = prod(h.Hx, XtHx);
	for (std::size_t i = 0; i < h.Zv.size(); ++i)
		S(i,i) += Float(h.Zv[i]);	// ISSUE mixed type proxy assignment

						// Inverse innovation covariance
	Float rcond = UdUinversePD (SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

						// Kalman gain, X*Hx'*SI
	noalias(W) = prod(XtHx, SI);

						// State update
	noalias(x) += prod(W, s);
	noalias(X) -= prod_SPD(W, S, XtHx);

	return rcond;
}


Bayes_base::Float
 Covariance_scheme::observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s)
{
					// Size consistency, z to model
	if (s.size() != h.Zv.size())
		error (Logic_exception("observation and model size inconsistent"));
	observe_size (s.size());// Dynamic sizing
	return byobserve_innovation (h, s, S, SI, W, XtHx);
}

Bayes_base::Float
 Covariance_scheme::observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s)
{
					// Size consistency, z to model
	if (s.size() != h.Z.size1())
		error (Logic_exception("observation and model size inconsistent"));
	observe_size (s.size());// Dynamic sizing
	return byobserve_innovation (h, s, S, SI, W, XtHx);
}

}//namespace
