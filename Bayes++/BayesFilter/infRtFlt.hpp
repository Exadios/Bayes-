#ifndef _BAYES_FILTER_INFORMATION_ROOT
#define _BAYES_FILTER_INFORMATION_ROOT

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Information Root Filter Scheme.
 *  A extended 'Square-root' Information filter as an Abstract class
 *
 * Algorithm: Square-root information propogation using QR factorisation
 * Ref:	P. Dyer and S. McReynolds, "Extension of Square-Root Filtering to Include Process Noise",
 * [1] Journal of Optimization Theory and Applications, Vol.3 No.6 1969
 * Filter maintains r,R where
 *   inv(R)*inv(R)' = X
 *   r = R*x
 *   R is upper triangular but not strictly a Cholesky factor as diagonal may be negative
 * Observe algorithm has been extended to include linearised models
 * Discontinous observe models require that state is normailised with respect to the observation.
 *
 * The filter is operated by performing a
 *  predict, observe
 * cycle defined by the base class
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

class Information_root_bscheme : public Extended_kalman_filter
{
public:
	FM::Vec r;			// Information Root state
	FM::UTriMatrix R;	// Information Root

	Information_root_bscheme (std::size_t x_size);

	void init ();
	void update ();
	// Covariance form state interface

	Float bypredict (Linear_invertible_predict_model& f)
	// Use linear form for r, and use inv.Fx from invertible model
	{
		return bypredict(f, f.inv.Fx, true);
	}
	Float bypredict (Linrz_predict_model& f, const FM::ColMatrix& invFx, bool linear_r);
	/* Explict form, using precomputed inverse of f.Fx */

	Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s);
	Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s);
	// Extended_kalman_filter observe

	static void inverse_Fx (FM::DenseColMatrix& invFx, const FM::Matrix& Fx);
	/* Numerical Inversion of Fx using LU factorisation */
};


class Information_root_scheme : public Information_root_bscheme
{
public:
	Information_root_scheme (std::size_t x_size);

	Float predict (Linrz_predict_model& f);
	// Extended_kalman_filter predict - use linrz form for r, computes inverse model using inverse_Fx
	Float predict (Linear_predict_model& f);
	/* Use linear form for r, computes inverse model using inverse_Fx */
};


/*
 * Information Root Filter Scheme with exposed information state
 * Augments Information_root_filter with y,Y in the interface
 */

class Information_root_info_scheme : public Information_root_scheme, virtual public Information_state
{
public:
	Information_root_info_scheme (std::size_t x_size);

	void init_yY ();
	void update_yY ();
	// Information form state interface
};


}//namespace
#endif
