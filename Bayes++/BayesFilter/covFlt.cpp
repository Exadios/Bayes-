/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Id$
 */

/*
 * Covariance Filter.
 */
#include "covFlt.hpp"
#include "matSup.hpp"

extern void mdebug(char*, Bayesian_filter_matrix::RowMatrix m);
extern void mdebug(char*, Bayesian_filter_matrix::SymMatrix m);
extern void mdebug(char*, Bayesian_filter_matrix::UTriMatrix m);
extern void mdebug(char*, Bayesian_filter_matrix::Vec v);

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


Covariance_bscheme::Covariance_bscheme (std::size_t x_size) :
	Kalman_state(x_size),
	tempX(x_size,x_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
}

Covariance_scheme::Covariance_scheme (std::size_t x_size, std::size_t z_initialsize) :
	Kalman_state(x_size),
	Covariance_bscheme(x_size),
	S(Empty), Sci(Empty), W(Empty), Wc(Empty), XtHx(Empty)
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

Covariance_bscheme& Covariance_bscheme::operator= (const Covariance_bscheme& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Kalman_state::operator=(a);
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
		Sci.resize(z_size,z_size, false);
		W.resize(x.size(),z_size, false);
		Wc.resize(x.size(),z_size,false);
		XtHx.resize(x.size(),z_size, false);
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
 				FM::SymMatrix& S, FM::UTriMatrix& Sci, FM::Matrix& W, FM::Matrix& Wc, FM::Matrix& XtHx)
/* Correlated innovation observe with explict byproduct
*/
{
						// Size consistency, z to model
	if (s.size() != h.Z.size1())
		error (Logic_exception("observation and model size inconsistent"));

						// Innovation covariance
	noalias(XtHx) = prod(X, trans(h.Hx));
	noalias(S) = prod(h.Hx, XtHx) + h.Z;

						// Inverse Cholesky factorisation
	Float rcond = UCfactor (Sci, S);
	rclimit.check_PD(rcond, "S not PD in observe");
	UTinverse (Sci);	// NOTE cannot be singular as PD
						
						// Kalman gain, X*Hx'*SI
	noalias(Wc) = prod(XtHx, trans(Sci));
	noalias(W) = prod(Wc, Sci);

						// State and Covariance update
	noalias(x) += prod(W, s);
	noalias(X) -= prod_SPD(Wc);

	return rcond;
}

Bayes_base::Float
 Covariance_bscheme::byobserve_innovation (Linrz_uncorrelated_observe_model& h, const Vec& s,
 				FM::SymMatrix& S, FM::UTriMatrix& Sci, FM::Matrix& W, FM::Matrix& Wc, FM::Matrix& XtHx)
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

						// Inverse Cholesky factorisation
	Float rcond = UCfactor (Sci, S);
	rclimit.check_PD(rcond, "S not PD in observe");
	UTinverse (Sci);	// NOTE cannot be singular as PD
						
						// Kalman gain, X*Hx'*SI
	noalias(Wc) = prod(XtHx, trans(Sci));
	noalias(W) = prod(Wc, Sci);

						// State and Covariance update
	noalias(x) += prod(W, s);
	noalias(X) -= prod_SPD(Wc);

	return rcond;
}


Bayes_base::Float
 Covariance_scheme::observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s)
{
					// Size consistency, z to model
	if (s.size() != h.Zv.size())
		error (Logic_exception("observation and model size inconsistent"));
	observe_size (s.size());// Dynamic sizing
	return byobserve_innovation (h, s, S, Sci, W, Wc, XtHx);
}

Bayes_base::Float
 Covariance_scheme::observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s)
{
					// Size consistency, z to model
	if (s.size() != h.Z.size1())
		error (Logic_exception("observation and model size inconsistent"));
	observe_size (s.size());// Dynamic sizing
	return byobserve_innovation (h, s, S, Sci, W, Wc, XtHx);
}

}//namespace
