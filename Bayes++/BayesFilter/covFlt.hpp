#ifndef _BAYES_FILTER_COVARIANCE
#define _BAYES_FILTER_COVARIANCE

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Covariance Filter.
 *	Implemention of extended Kalman filter
 * 
 * To work with with Linear and Linrz models
 *  a) a state seperate from covariance prediction is used.
 *  b) a EKF innovation update algorithm is used.
 * Discontinous observe models require that prediction is normailised with
 * respect to the observation.
 *
 * A initial observation size may also be specified for efficiency.
 * 
 * The filter is operated by performing a
 * 	predict, observe
 * cycle defined by the base class
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

class Covariance_filter : public Extended_filter
{
public:
	Covariance_filter (size_t x_size, size_t z_initialsize = 0);
	Covariance_filter& operator= (const Covariance_filter&);
	// Optimise copy assignment to only copy filter state

	void init ();
	void update ();
	Float predict (Linrz_predict_model& f);
	Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s);
	Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s);

public:						// Exposed Numerical Results
	FM::SymMatrix S, SI;		// Innovation Covariance and Inverse

protected:					// allow fast operation if z_size remains constant
	size_t last_z_size;
	void observe_size (size_t z_size);
};


}//namespace
#endif
