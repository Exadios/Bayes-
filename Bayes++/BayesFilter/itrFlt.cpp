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
 */
#include "bayesFlt.hpp"
#include "matSup.hpp"
#include "itrFlt.hpp"
#include "models.hpp"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


bool Counted_iterated_terminator::term_or_relinearize (const Iterated_covariance_scheme& f)
{
	--i;
	if (i != 0)
		m.relinearise (f.x);
	return (--i) == 0;
}


Iterated_covariance_scheme::Iterated_covariance_scheme(size_t x_size, size_t z_initialsize) :
		Linrz_kalman_filter(x_size),
		S(Empty), SI(Empty),
		s(Empty), HxT(Empty)
/*`
 * Initialise filter and set the size of things we know about
 */
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

Iterated_covariance_scheme&
 Iterated_covariance_scheme::operator= (const Iterated_covariance_scheme& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Kalman_state_filter::operator=(a);
	return *this;
}

void Iterated_covariance_scheme::init ()
{
						// Postconditions
	if (!isPSD (X))
		filter_error ("Initial X not PSD");
}

void Iterated_covariance_scheme::update ()
{
	// Nothing to do, implicit in observation
}

Bayes_base::Float
 Iterated_covariance_scheme::predict (Linrz_predict_model& f)
{
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict state covariance
	RowMatrix temp_FxX(f.Fx.size1(),X.size2());
	X = prod_SPD(f.Fx,X, temp_FxX) + prod_SPD(f.G,f.q);

	assert_isPSD (X);
	return 1;
}


void Iterated_covariance_scheme::observe_size (size_t z_size)
/*
 * Optimised dynamic observation sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		s.resize(z_size);
		S.resize(z_size,z_size);
		SI.resize(z_size,z_size);
		HxT.resize(x.size(),z_size);
	}
}

Bayes_base::Float
 Iterated_covariance_scheme::observe (Linrz_uncorrelated_observe_model& h, Iterated_terminator& term, const FM::Vec& z)
/*
 * Iterated Extended Kalman Filter
 * Bar-Shalom and Fortmann p.119 (full scheme)
 * A hard limit is placed on the iterations whatever the
 * the normal terminal condition is to guarantee termination
 * Uncorrelated noise
 */
{
						// ISSUE: Implement simplified uncorrelated noise equations
	Adapted_Linrz_correlated_observe_model hh(h);
	return observe (hh, term, z);
}

Bayes_base::Float
 Iterated_covariance_scheme::observe (Linrz_correlated_observe_model& h, Iterated_terminator& term, const FM::Vec& z)
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
	SymMatrix ZI(z_size,z_size);

							// Dynamic sizing
	observe_size (z.size());

	Vec xpred = x;			// Initialise iteration
	SymMatrix Xpred = X;
							// Inverse predicted covariance
	SymMatrix XpredI(x_size,x_size);
	Float rcond = UdUinversePD (XpredI, Xpred);
	rclimit.check_PD(rcond, "Xpred not PD in observe");

							// Inverse observation covariance
	rcond = UdUinversePD (ZI, h.Z);
	rclimit.check_PD(rcond, "Z not PD in observe");
				
	RowMatrix HxXtemp(h.Hx.size1(),X.size2());
	RowMatrix temp1(x_size,x_size), temp2(x_size,z_size);
	SymMatrix temp3(x_size,x_size);

	do {
							// Observation model, linearize about new x
		const Vec& zp = h.h(x);
		
		HxT.assign (trans(h.Hx));
							// Innovation
		h.normalise(s = z, zp);
		s.minus_assign (zp);
							// Innovation covariance
		S = prod_SPD(h.Hx, Xpred, HxXtemp) + h.Z;
							// Inverse innovation covariance
		rcond = UdUinversePD (SI, S);
		rclimit.check_PD(rcond, "S not PD in observe");

							// Iterative observe
		temp3.assign (prod_SPD(HxT,SI, temp2));
		temp1.assign (prod(Xpred,temp3));
		X.assign (Xpred - prod(temp1,Xpred));
		assert_isPSD (X);
							// New state iteration
		temp2.assign (prod(X,HxT));
		temp1.assign (prod(X,XpredI));
		x += prod(temp2,prod(ZI,s)) - prod(temp1, (x - xpred));
	} while (!term.term_or_relinearize(*this));
	return rcond;
}

}//namespace
