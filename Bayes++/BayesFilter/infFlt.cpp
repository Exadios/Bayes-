/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Information Filter.
 */
#include "infFlt.hpp"
#include "matSup.hpp"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


Information_scheme::Information_scheme (size_t x_size) :
		Kalman_state_filter(x_size), Information_state_filter(x_size),
		tempX(x_size,x_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
	update_required = true;	// Not a valid state, init is required before update can be used
}

Information_scheme::Predict_linear_byproduct::Predict_linear_byproduct (size_t x_size, size_t q_size) :
/* Set size of by-products for linear predict
 */
		 A(x_size,x_size), tempG(x_size,q_size),
		 B(q_size,q_size),
		 y(x_size)
{}

Information_scheme& Information_scheme::operator= (const Information_scheme& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Information_state_filter::operator=(a);
	Kalman_state_filter::operator=(a);
	return *this;
}


void Information_scheme::init ()
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
	noalias(y) = prod(Y,x);
	update_required = false;
}

void Information_scheme::init_yY ()
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
		error (Numeric_exception("Initial Y not PSD"));
	update_required = true;
}

void Information_scheme::update_yY ()
/*
 * Postcondition:
 *		y, Y is PSD
 */
{
}

void Information_scheme::update ()
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

		noalias(x) = prod(X,y);
		update_required = false;
	}
}

Bayes_base::Float
 Information_scheme::predict (Linrz_predict_model& f)
/*
 * Extended_kalman_filter predict via state
 *  Computation is through state to accommodate linearied model
 */
{
	update ();			// x,X required
	x = f.f(x);			// Extended Kalman state predict is f(x) directly
						// Predict information matrix, and state covariance
	noalias(X) = prod_SPD(f.Fx,X,tempX) + prod_SPD(f.G, f.q);

						// Information
	Float rcond = UdUinversePD (Y, X);
	rclimit.check_PD(rcond, "X not PD in predict");
						// Predict information state
	noalias(y) = prod(Y,x);
	return rcond;
}

Bayes_base::Float
 Information_scheme::predict (Linear_invertable_predict_model& f)
/* Linear information predict, without byproduct
 */
{
	Predict_linear_byproduct b(f.Fx.size1(),f.q.size());
	return epredict (f, b);	
}

Bayes_base::Float
 Information_scheme::epredict (Linear_invertable_predict_model& f, Predict_linear_byproduct& b)
/*
 * Linear information predict
 *  Computation is through information state y,Y only
 *  Uses x(k+1|k) = Fx * x(k|k) instead of extended x(k+1|k) = f(x(k|k))
 * Requires y(k|k), Y(k|k)
 * Predicts y(k+1|k), Y(k+1|k)
 *
 * The numerical solution used is particularly flexible. It takes
 * particular care to avoid invertibilty requirements for the noise and noise coupling g,Q
 * Therefore both zero noises and zeros in the couplings can be used
 */
{
						// A = invFx'*Y*invFx ,Inverse Predict covariance
	noalias(b.A) = prod_SPDT(f.inv.Fx, Y, tempX);
						// B = G'*A*G+invQ , A in coupled additive noise space
	noalias(b.B) = prod_SPDT(f.G, b.A, b.tempG);
	for (size_t i = 0; i < f.q.size(); ++i)
	{
		if (f.q[i] < 0)	// allow PSD q, let infinity propogate into B
			error (Numeric_exception("Predict q Not PSD"));
		b.B(i,i) += Float(1) / f.q[i];
	}

						// invert B ,Addative noise
	Float rcond = UdUinversePDignoreInfinity (b.B);
	rclimit.check_PD(rcond, "(G'invFx'.Y.invFx.G + invQ) not PD in predict");

						// G*invB*G' ,in state space
	noalias(Y) = prod_SPD(f.G,b.B, b.tempG);
						// I - A* G*invB*G', information gain
	FM::identity(tempX);
	noalias(tempX) -= prod(b.A,Y);
						// Information
	noalias(Y) = prod(tempX,b.A);
						// Information state
	noalias(b.y) = prod(trans(f.inv.Fx), y);
	noalias(y) = prod(tempX, b.y);

	update_required = true;
	return rcond;
}


Bayes_base::Float
 Information_scheme::observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s)
/* Extended_kalman_filter observe, unused byproduct
 */
{
	const size_t x_size = h.Hx.size2();
	State_byproduct i(x_size);
	Covariance_byproduct I(x_size,x_size);
	return eobserve_innovation (h, s, i,I);
}	

Bayes_base::Float
 Information_scheme::observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s)
/* Extended_kalman_filter observe, unused byproduct
 */
{
	const size_t x_size = h.Hx.size2();
	State_byproduct i(x_size);
	Covariance_byproduct I(x_size,x_size);
	return eobserve_innovation (h, s, i,I);
}	

Bayes_base::Float
 Information_scheme::eobserve_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s, State_byproduct& i, Covariance_byproduct& I)
/* Correlated innovation observe with explict byproduct
 */
{
						// Size consistency, z to model
	if (s.size() != h.Z.size1())
		error (Logic_exception("observation and model size inconsistent"));

	Vec zz(s + prod(h.Hx,x));		// Strange EIF observation object

						// Observation Information
	SymMatrix ZI(s.size(),s.size());
	Float rcond = UdUinversePD (ZI, h.Z);
	rclimit.check_PD(rcond, "Z not PD in observe");

	RowMatrix HxT (trans(h.Hx));
	RowMatrix HxTZI (prod(HxT, ZI));
												// Calculate EIF i = Hx'*ZI*zz
	noalias(i) = prod(HxTZI, zz);
												// Calculate EIF I = Hx'*ZI*Hx
	noalias(I) = prod(HxTZI, trans(HxT));				// use column matrix trans(HxT)

	noalias(y) += i;
	noalias(Y) += I;
	update_required = true;

	return rcond;
}

Bayes_base::Float
 Information_scheme::eobserve_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s, State_byproduct& i, Covariance_byproduct& I)
/* Uncorrelated innovation observe with explict byproduct
 */
{
						// Size consistency, z to model
	if (s.size() != h.Zv.size())
		error (Logic_exception("observation and model size inconsistent"));

	Vec zz(s + prod(h.Hx,x));		// Strange EIF observation object

						// Observation Information
	Float rcond = UdUrcond(h.Zv);
	rclimit.check_PD(rcond, "Zv not PD in observe");

	RowMatrix HxT (trans(h.Hx));      			// HxT = Hx'*inverse(Z)
	for (size_t w = 0; w < h.Zv.size(); ++w)
		column(HxT, w) *= 1 / h.Zv[w];
												// Calculate EIF i = Hx'*ZI*zz
	noalias(i) = prod(HxT, zz);
												// Calculate EIF I = Hx'*ZI*Hx
	noalias(I) = prod(HxT, h.Hx);

	noalias(y) += i;
	noalias(Y) += I;
	update_required = true;

	return rcond;
}

}//namespace
