#ifndef _BAYES_FILTER_ITERATED_COVARIANCE
#define _BAYES_FILTER_ITERATED_COVARIANCE

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
 *
 * The observe algorithm uses the iterated non-linear formulation 
 * from Bar-Shalom and Fortmann p.119 (full scheme)
 * Discontinous observe models require that prediction is normailised with
 * respect to the observation.
 *
 * The filter is operated by performing a
 * 	predict, observe
 * cycle defined by the base class
 */
#include "covFlt.h"

/* Filter namespace */
namespace Bayesian_filter
{

class Iterated_covariance_filter : public Linrz_filter
{
public:
	Iterated_covariance_filter (size_t x_size, size_t z_initialsize = 0);
	/* Initialised filter requries an addition iteration limit for the
	   observe algorithm */
	Iterated_covariance_filter& operator= (const Iterated_covariance_filter&);
	// Optimise copy assignment to only copy filter state

	using Linrz_filter::init;
	void init ();
	void update ();
	Float predict (Linrz_predict_model& f);

	Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z);
	Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z);

public:						// Exposed Numerical Results
	FM::SymMatrix S, SI;		// Innovation Covariance and Inverse

protected:
	virtual bool observe_iteration_end () = 0;
	/* An algorithm must be supplied to signal end of observe iteration
	 * Precond: x,X,s,S computed for the last iteration
	 * Iteration count should be reset on construction and when true is returned
	 */

protected:					// allow fast operation if z_size remains constant
	size_t last_z_size;
	void observe_size (size_t z_size);
							// Permenantly allocated temps
	FM::Vec s;
	FM::Matrix HxT;
};


}//namespace
#endif
