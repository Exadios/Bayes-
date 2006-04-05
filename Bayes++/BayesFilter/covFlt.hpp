#ifndef _BAYES_FILTER_COVARIANCE
#define _BAYES_FILTER_COVARIANCE

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Id$
 */

/*
 * Covariance Filter Scheme.
 *  Implemention of extended Kalman filter
 * 
 * To work with with Linear and Linrz models
 *  a) a state seperate from covariance predict is used.
 *  b) a EKF innovation update algorithm is used.
 * Discontinous observe models require that predict is normailised with
 * respect to the observation.
 *
 * A initial observation size may also be specified for efficiency.
 * 
 * The filter is operated by performing a
 * predict, observe
 * cycle defined by the base class
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

class Covariance_bscheme : virtual public Kalman_state
{
public:
	Covariance_bscheme (std::size_t x_size);
	Covariance_bscheme& operator= (const Covariance_bscheme&);
	// Optimise copy assignment to only copy filter state

	void init ();
	void update ();

	Float predict (Gaussian_predict_model& f);
	// Specialised 'stationary' predict, only additive noise

	Float byobserve_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s,
				FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W, FM::Matrix& XtHx);
	Float byobserve_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s,
				FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W, FM::Matrix& XtHx);
	// Observe with explict byproduct

protected:			   		// Permenantly allocated temps
	FM::RowMatrix tempX;
};


class Covariance_scheme : public Extended_kalman_filter, public Covariance_bscheme
{
public:
	Covariance_scheme (std::size_t x_size, std::size_t z_initialsize = 0);

	Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s);
	Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s);
	// Extended_kalman_filter observe

	Float predict (Linrz_predict_model& f);
	// Extended_kalman_filter predict

protected:					// Allow fast operation if z_size remains constant
	std::size_t last_z_size;
	void observe_size (std::size_t z_size);
							// Numerical byproducts
	FM::SymMatrix S, SI;		// Innovation Covariance and Inverse
	FM::Matrix W;				// Kalman Gain
	FM::Matrix XtHx;			// X * Hx'
};

}//namespace
#endif
