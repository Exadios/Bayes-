#ifndef _BAYES_FILTER_CI
#define _BAYES_FILTER_CI

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Covariance Intersection Filter Scheme.
 *
 * References
 *  [1] "A Non divergent Estimation Algorithm in the Presence of Unknown Correlations"
 *   Simon J Julier, Jeffrey K Uhlmann
 *
 * CI provides a generalised consistent method to combine mean and covariances
 * of two estimates. The combination can be optimised by chosing a norm of the
 * combined correlations. The norm (omega) is restrict to 0..1 inclusive an effectively
 * scales the combination
 * Here is CI with a predict and obeserve model to form a filter.
 *
 * The Omega norm chosen here is the fixed value of 0.5
 * The Omega function should be overloaded to produce more useful results
 *
 * The filter is operated by performing a
 *  predict, observe
 * cycle derived from the Extended_filter
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

class CI_bscheme : public Extended_kalman_filter
{
public:
	CI_bscheme (std::size_t x_size);
	CI_bscheme& operator= (const CI_bscheme&);
	// Optimise copy assignment to only copy filter state

	void init ();
	void update ();
	Float predict (Linrz_predict_model& f);

	Float byobserve_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s,
				FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W);
	Float byobserve_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s,
				FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W);
	// Observe with explict byproduct

	virtual Float Omega(const FM::SymMatrix& Ai, const FM::SymMatrix& Bi, const FM::SymMatrix& A)
	// Determine norm Omega 0..1 for the CI combination
	// Default norm is the fixed value 0.5
	{
		return 0.5;
	}

};


class CI_scheme : public CI_bscheme
{
public:
	CI_scheme (std::size_t x_size);

	Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s)
	{	// Extended_kalman_filter observe
		const std::size_t z_size = h.Hx.size1();
		FM::SymMatrix S(z_size,z_size), SI(z_size,z_size);
		FM::Matrix W(h.Hx.size2(), z_size);
		return byobserve_innovation (h, s, S,SI,W);
	}
	Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s)
	{	// Extended_kalman_filter observe
		const std::size_t z_size = h.Hx.size1();
		FM::SymMatrix S(z_size,z_size), SI(z_size,z_size);
		FM::Matrix W(h.Hx.size2(), z_size);
		return byobserve_innovation (h, s, S,SI,W);
	}
};

}//namespace
#endif
