/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Covariance Intersection Filter.
 * TODO: Implement useful Omega based on iterative optimization algorithm from the authors of reference [1]
 */
#include "CIFlt.hpp"
#include "matSup.hpp"
#include "models.hpp"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


CI_scheme::CI_scheme (size_t x_size) :
	Kalman_state_filter(x_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
}

CI_scheme& CI_scheme::operator= (const CI_scheme& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Kalman_state_filter::operator=(a);
	return *this;
}


void CI_scheme::init ()
{
						// Postconditions
	if (!isPSD (X))
		error (Numeric_exception("Initial X not PSD"));
}

void CI_scheme::update ()
{
	// Nothing to do, implicit in observation
}

Bayes_base::Float
 CI_scheme::predict (Linrz_predict_model& f)
{
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict state covariance
	RowMatrix tempX(f.Fx.size1(), X.size2());
	noalias(X) = prod_SPD(f.Fx,X, tempX);
	noalias(X) += prod_SPD(f.G, f.q, tempX);

	return 1;
}


Bayes_base::Float
 CI_scheme::eobserve_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s,
				Covariance_byproduct& S, Kalman_gain_byproduct& b)
/*
 * Iterated Extended Kalman Filter
 * Bar-Shalom and Fortmann p.119 (full scheme)
 * A hard limit is placed on the iterations whatever the
 * the normal terminal condition is to guarantee termination
 * Uncorrelated noise
 */
{
						// ISSUE: Implement simplified uncorrelated noise equations
	size_t z_size = s.size();
	SymMatrix Z(z_size,z_size);

	Adapted_Linrz_correlated_observe_model hh(h);
	return eobserve_innovation (hh, s, S, b);
}


Bayes_base::Float
 CI_scheme::eobserve_innovation (Linrz_correlated_observe_model& h, const Vec& s,
				Covariance_byproduct& S, Kalman_gain_byproduct& b)
/* Correlated innovation observe
 */
{
	const Float one = 1;
						// size consistency, z to model
	if (s.size() != h.Z.size1())
		error (Logic_exception("observation and model size inconsistent"));

						// Linear conditioning for omega
	SymMatrix invZ(h.Z.size1(),h.Z.size2());
	Float rcond = UdUinversePD (invZ, h.Z);
	rclimit.check_PSD(rcond, "Z not PSD in observe");

	Matrix HTinvZ (prod(trans(h.Hx), invZ));
	SymMatrix HTinvZH (prod(HTinvZ, h.Hx));

	SymMatrix invX(X.size1(),X.size2());
	rcond = UdUinversePD (invX, X);
	rclimit.check_PD(rcond, "X not PD in observe");


						// find omeage
	Float omega = Omega(invX, HTinvZH, X);

						// calculate predicted innovation
	Matrix XHT (prod(X, trans(h.Hx)));
	Matrix HXHT (prod(h.Hx, XHT));
	S = HXHT * (one-omega) + h.Z * omega;

						// inverse innovation covariance
	rcond = UdUinversePD (b.SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

	Matrix K (prod(XHT*(one-omega), b.SI));

						// state update
	noalias(x) += prod(K, s);
						// inverse covariance
	invX *= omega;						
	noalias(invX) += HTinvZH*(one-omega);
						// covariance
	rcond = UdUinversePD (X, invX);
	rclimit.check_PD(rcond, "inverse covariance not PD in observe");
	return rcond;
}


}//namespace
