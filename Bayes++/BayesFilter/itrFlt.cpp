/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Iterated Covariance Filter.
 *	A non-linear Covariance (Kalman) filter as an Abstract class
 */
#include "matSup.h"
#include "itrFlt.h"
#include "models.h"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


Iterated_covariance_filter::Iterated_covariance_filter(size_t x_size, size_t z_initialsize) :
		Linrz_filter(x_size),
		S(Empty), SI(Empty),
		s(Empty), HxT(Empty)
/*`
 * Initialise filter and set the size of things we know about
 */
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

Iterated_covariance_filter&
 Iterated_covariance_filter::operator= (const Iterated_covariance_filter& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Linrz_filter::operator=(a);
	return *this;
}

void Iterated_covariance_filter::init ()
{
						// Preconditions
	if (!isPSD (X))
		filter_error ("Xi not PSD");
}

void Iterated_covariance_filter::update ()
{
	// Nothing to do, implicit in observation
}

Bayes_base::Float
 Iterated_covariance_filter::predict (Linrz_predict_model& f)
{
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict state covariance
	X = mult_SPD(f.Fx, X) + mult_SPD(f.G, f.q);

	// TODO: Remove check once proved
	assert_isPSD (X);
	return 1.;
}


void Iterated_covariance_filter::observe_size (size_t z_size)
/*
 * Optimised dynamic observation sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		S.resize(z_size,z_size);
		SI.resize(z_size,z_size);
		HxT.resize(x.size(),z_size);
	}
}

Bayes_base::Float
 Iterated_covariance_filter::observe (Linrz_uncorrelated_observe_model& h, const Vec& z)
/*
 * Iterated Extended Kalman Filter
 * Bar-Shalom and Fortmann p.119 (full scheme)
 * A hard limit is placed on the iterations whatever the
 * the normal terminal condition is to guarantee termination
 * Uncorrelated noise
 */
{
						// ISSUE: Implement simplied uncorrelated noise equations
	size_t z_size = z.size();
	SymMatrix Z(z_size,z_size);
	
	Adapted_Linrz_correlated_observe_model hh(h);
	return observe (hh, z);
}

Bayes_base::Float
 Iterated_covariance_filter::observe (Linrz_correlated_observe_model& h, const Vec& z)
/*
 * Iterated Extended Kalman Filter
 * Bar-Shalom and Fortmann p.119 (full scheme)
 * A hard limit is placed on the iterations whatever the
 * the normal terminal condition is to guarantee termination
 * returned rcond is of S (or 1 if no iterations are performed)
 */
{
	size_t x_size = x.size();
	size_t z_size = z.size();
	Vec xpred(x_size);
	SymMatrix Xpred(x_size,x_size);
	SymMatrix XpredI(x_size,x_size);
	SymMatrix ZI(z_size,z_size);

						// Dynamic sizing
	observe_size (z.size());

	xpred = x;				// Initialise iteration
	Xpred = X;
							// Inverse predicted covariance
	Float rcond = UdUinversePD (XpredI, Xpred);
	rclimit.check_PD(rcond, "Xpred not PD in observe");

							// Inverse observation covariance
	rcond = UdUinversePD (ZI, h.Z);
	rclimit.check_PD(rcond, "Z not PD in observe");
				
	do {
		Vec zp = h.h(x);		// Observation model
		h.normalise(zp, z);
		HxT.assign (trans(h.Hx));
		s = z - zp;			// Innovation and covariance
		S = mult_SPD(h.Hx, Xpred) + h.Z;

							// Inverse innovation covariance
		rcond = UdUinversePD (SI, S);
		rclimit.check_PD(rcond, "S not PD in observe");

							// Iterative observe
		X = Xpred - mult_SPD (Xpred, mult_SPD(HxT,SI));
		// TODO: Remove check once proved
		assert_isPSD (X);
		x += prod(prod(X,HxT),prod(ZI,s)) - prod(prod(X,XpredI), (x - xpred));
	} while (!observe_iteration_end());
	return rcond;
}

}//namespace
