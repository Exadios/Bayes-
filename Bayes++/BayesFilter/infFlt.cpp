/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Information Filter.
 */
#include "bayesFlt.hpp"
#include "matSup.hpp"
#include "infFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


Information_filter::Information_filter (size_t x_size, size_t z_initialsize) :
		Information_form_filter(x_size), Extended_filter(x_size),
		i(x_size), I(x_size,x_size),
		ZI(Empty)
/*
 * Initialise filter and set the size of things we know about
 */
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
	update_required = true;	// Not a valid state, init is required before update can be used
}

Information_filter::Linear_predict_byproducts::Linear_predict_byproducts (size_t x_size, size_t q_size) :
/* Set size of by-products for linear predict
 */
		 tempA(x_size,x_size), A(x_size,x_size), tempG(q_size,x_size),
		 B(q_size,q_size), tempY(x_size,x_size),
		 invq(q_size)
{}

Information_filter& Information_filter::operator= (const Information_filter& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Extended_filter::operator=(a);
	y = a.y;
	Y = a.Y;
	return *this;
}


void Information_filter::init ()
/*
 * Initialise the filter from x,X
 * Precondition:
 *		x, X
 * Postcondition:
 *		x, X is PD
 *		y, Y is PSD
 */
{
						// Information
	Float rcond = UdUinversePD (Y, X);
	rclimit.check_PD(rcond, "Initial X not PD");
						// Information state
	y.assign (prod(Y,x));
	update_required = false;
}

void Information_filter::init_yY ()
/*
 * Initialisation directly from Information
 * Precondition:
 *		y, Y
 * Postcondition:
 *		y, Y is PSD
 */
{
						// Postconditions
	if (!isPSD (Y))
		filter_error ("Initial Y not PSD");
	update_required = true;
}

void Information_filter::update_yY ()
/*
 * Postcondition:
 *		y, Y is PSD
 */
{
}

void Information_filter::update ()
/*
 * Recompute x,X from y,Y
 *  Optimised using update_required (postcondition met iff update_required false)
 * Precondition:
 *		y, Y is PD
 * Postcondition:
 *		x=X*y, X=inv(Y) is PSD
 *		y, Y is PD
 */
{
	if (update_required)
	{		// Covariance
		Float rcond = UdUinversePD (X, Y);
		rclimit.check_PD(rcond, "Y not PD in update");

		x.assign (prod(X,y));
		update_required = false;
	}
}

Bayes_base::Float
 Information_filter::predict (Linrz_predict_model& f)
/*
 * Extented linrz information prediction
 *  Computation is through state to accommodate linearied model
 */
{
	update ();			// x,X required
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict information matrix, and state covariance
	RowMatrix temp(f.Fx.size1(), X.size2());
	X = prod_SPD(f.Fx,X, temp) + prod_SPD(f.G, f.q);

						// Information
	Float rcond = UdUinversePD (Y, X);
	rclimit.check_PD(rcond, "X not PD in predict");
						// Predict information state
	y.assign (prod(Y,x));
	return rcond;
}

Float Information_filter::predict (Linear_invertable_predict_model& f, Linear_predict_byproducts& b)
/*
 * Linear information prediction
 *  Computation is through information state only
 *  Uses x(k+1|k) = Fx * x(k|k) instead of extended x(k+1|k) = f(x(k|k))
 * Prediction is done completely on y,Y
 * Requires y(k|k), Y(k|k)
 * Predicts y(k+1|k), Y(k+1|k)
 * 
 * The numerical solution used is particularly flexible. It takes
 * particular care to avoid invertibilty requirements for the noise and noise coupling g,Q
 * Therefore both zero noises and zeros in the couplings can be used
 */
{
						// A = invFx'*Y*invFx ,Inverse Predict covariance
	b.tempA.assign (prod(Y, f.inv.Fx));
	b.A.assign (prod(trans(f.inv.Fx), b.tempA));
						// B = G'*A*G+invQ , A in coupled additive noise space
	b.tempG.assign (prod(trans(f.G), b.A));
	b.B.assign (prod(b.tempG, f.G));
	for (size_t i = 0; i < f.q.size(); ++i)
	{
		if (f.q[i] < 0)	// allow PSD q, let infinity propogate into B
			filter_error ("Predict q Not PSD");
		b.invq[i] = Float(1) / f.q[i];
	}
	diag(b.B).plus_assign (b.invq);
	
						// invert B ,Addative noise
	Float rcond = UdUinversePDignoreInfinity (b.B);
	rclimit.check_PD(rcond, "(G'invFx'.Y.invFx.G + invQ) not PD in predict");

						// G*invB*G' ,in state space
	b.tempG.assign (prod(b.B,trans(f.G)));
	Y.assign (prod(f.G,b.tempG));
						// I - A* G*invB*G', information gain
	FM::identity(b.tempY);
	b.tempY.minus_assign (prod(b.A,Y));
						// Information
	Y.assign (prod(b.tempY,b.A));
						// Information state
	y = prod(prod(b.tempY,trans(f.inv.Fx)), y);
	
	update_required = true;
	assert_isPSD (Y);
	return rcond;
}


inline void Information_filter::observe_size (size_t z_size)
/*
 * Optimised dynamic observation sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		ZI.resize(z_size,z_size);
	}
}

Bayes_base::Float
 Information_filter::observe_innovation (Linrz_correlated_observe_model& h, const Vec& s)
/* correlated innovation observe
 */
{
						// Size consistency, z to model
	if (s.size() != h.Z.size1())
		filter_error("observation and model size inconsistent");
	observe_size (s.size());// Dynamic sizing

	Vec zz(s + prod(h.Hx,x));		// Strange EIF obsevation object object

						// Observation Information
	Float rcond = UdUinversePD (ZI, h.Z);
	rclimit.check_PD(rcond, "Z not PD in observe");

	RowMatrix HxT (trans(h.Hx));
	RowMatrix HxTZI (prod(HxT, ZI));
												// Calculate EIF i = Hx'*ZI*zz
	i.assign (prod(HxTZI, zz));
												// Calculate EIF I = Hx'*ZI*Hx
	I.assign (prod(HxTZI, trans(HxT)));				// use column matrix trans(HxT)

	y.plus_assign (i);
	Y.plus_assign (I);
	update_required = true;

	assert_isPSD (Y);
	return rcond;
}

Bayes_base::Float
 Information_filter::observe_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s)
/* Extended linrz uncorrelated observe
 */
{
						// Size consistency, z to model
	if (s.size() != h.Zv.size())
		filter_error("observation and model size inconsistent");
	observe_size (s.size());// Dynamic sizing

	Vec zz(s + prod(h.Hx,x));		// Strange EIF obsevation object object

						// Observation Information
	Float rcond = UdUrcond(h.Zv);
	rclimit.check_PD(rcond, "Zv not PD in observe");

	RowMatrix HxT (trans(h.Hx));      			// HxT = Hx'*inverse(Z)
	for (size_t w = 0; w < h.Zv.size(); ++w)
		column(HxT, w) *= 1 / h.Zv[w];
												// Calculate EIF i = Hx'*ZI*zz
	i.assign (prod(HxT, zz));
												// Calculate EIF I = Hx'*ZI*Hx
	I.assign (prod(HxT, h.Hx));

	y.plus_assign (i);
	Y.plus_assign (I);
	update_required = true;

	assert_isPSD (Y);
	return rcond;
}

}//namespace
