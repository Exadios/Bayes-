/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Covariance Intersection Filter.
 * TODO: Implement useful Omega based on iterative optimization algorithm from the authors of reference [1]
 */
#include "bayesFlt.hpp"
#include "CIFlt.hpp"
#include "matSup.hpp"
#include "models.hpp"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


CI_filter::CI_filter (size_t x_size, size_t z_initialsize) :
	Extended_filter(x_size),
	S(Empty), SI(Empty)
/*
 * Initialise filter and set the size of things we know about
 */
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

CI_filter& CI_filter::operator= (const CI_filter& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Extended_filter::operator=(a);
	return *this;
}


void CI_filter::init ()
{
						// Postconditions
	if (!isPSD (X))
		filter_error ("Initial X not PSD");
}

void CI_filter::update ()
{
	// Nothing to do, implicit in observation
}

Bayes_base::Float
 CI_filter::predict (Linrz_predict_model& f)
{
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict state covariance
	RowMatrix temp(f.Fx.size1(), X.size2());
	X = prod_SPD(f.Fx,X, temp) + prod_SPD(f.G, f.q);

	assert_isPSD (X);
	return 1;
}

void CI_filter::observe_size (size_t z_size)
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
 CI_filter::observe_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s)
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
	return observe_innovation (hh, s);
}


Bayes_base::Float
 CI_filter::observe_innovation (Linrz_correlated_observe_model& h, const Vec& s)
/* Correlated innovation observe
 */
{
	const Float one = 1;
						// size consistency, z to model
	if (s.size() != h.Z.size1())
		filter_error("observation and model size inconsistent");
	observe_size (s.size());	// dynamic sizing

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
	rcond = UdUinversePD (SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");

	Matrix K (prod(XHT*(one-omega), SI));

						// state update
	x += prod(K, s);
						// inverse covariance
	invX = invX * omega + HTinvZH*(one-omega);
						// covariance
	rcond = UdUinversePD (X, invX);
	rclimit.check_PD(rcond, "inverse covariance not PD in observe");
	return rcond;
}


}//namespace
