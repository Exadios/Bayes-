#ifndef _BAYES_FILTER_COVARIANCE
#define _BAYES_FILTER_COVARIANCE

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.htm for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
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

class Covariance_scheme : public Extended_kalman_filter
{
public:
	Covariance_scheme (size_t x_size);
	Covariance_scheme& operator= (const Covariance_scheme&);
	// Optimise copy assignment to only copy filter state

	void init ();
	void update ();

	Float predict (Linrz_predict_model& f);
	// Extended_kalman_filter predict
	Float predict (Gaussian_predict_model& f);
	// Specialised 'stationary' predict, only addative noise

	Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s);
	Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s);
	// Extended_kalman_filter observe
	Float eobserve_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s,
				Covariance_byproduct& S, Kalman_gain_byproduct& b);
	Float eobserve_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s,
				Covariance_byproduct& S, Kalman_gain_byproduct& b);
	// Observe with explict byproduct

protected:			   		// Permenantly allocated temps
	FM::RowMatrix tempX;
};


}//namespace
#endif
