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
 *  A non-linear Covariance (Kalman) filter as an Abstract class
 *
 * The observe algorithm uses the iterated non-linear formulation 
 * from Bar-Shalom and Fortmann p.119 (full scheme)
 * Discontinous observe models require that prediction is normailised with
 * respect to the observation.
 *
 * The filter is operated by performing a
 *  predict, observe
 * cycle defined by the base class
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

class Iterated_terminator : public Bayes_base
/*
 * Termination condition for filter Iteration
 *  Used by iterated observe to parameterise termination condition
 *
 * Defaults to immidately terminating the iteration
 *
 * A more useful terminator can built by derivation.
 * For example terminator constructed with a reference to the filter can
 * detect convergence of x and/or X
 */
{
public:
	virtual void relinearize (const FM::Vec& x)
	{}	// linearize observation model about new x
	virtual bool term_or_relinearize ()
	{
		return true;
	}
};

class Counted_iterated_terminator : public Bayes_base
/*
 * Termination condition for filter Iteration
 *  Used by iterated observe to parameterise termination condition
 *
 */
{
public:
	virtual void relinearize (const FM::Vec& x)
	{}	// linearize observation model about new x
	virtual bool term_or_relinearize ()
	{
		return true;
	}
};



class Iterated_covariance_filter : public Linrz_filter
{
public:
	Iterated_covariance_filter (size_t x_size, size_t z_initialsize = 0);
	/* Initialised filter requries an addition iteration limit for the
	   observe algorithm */
	Iterated_covariance_filter& operator= (const Iterated_covariance_filter&);
	// Optimise copy assignment to only copy filter state

	void init ();
	void update ();
	Float predict (Linrz_predict_model& f);

	Float observe (Linrz_uncorrelated_observe_model& h, Iterated_terminator& term, const FM::Vec& z);
	Float observe (Linrz_correlated_observe_model& h, Iterated_terminator& term, const FM::Vec& z);
	// Observe with iteration
	Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z)
	{	// Observe with default termination
		Iterated_terminator term;
		return observe (h, term, z);
	}
	Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z)
	{	// Observe with default termination
		Iterated_terminator term;
		return observe (h, term, z);
	}

public:						// Exposed Numerical Results
	FM::SymMatrix S, SI;		// Innovation Covariance and Inverse

protected:					// allow fast operation if z_size remains constant
	size_t last_z_size;
	void observe_size (size_t z_size);
							// Permenantly allocated temps
	FM::Vec s;
	FM::Matrix HxT;
};


}//namespace
#endif
