#ifndef _BAYES_FILTER_ITERATED_COVARIANCE
#define _BAYES_FILTER_ITERATED_COVARIANCE

/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
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
 * Derived filters must supply predict and observe model functions.
 * State and control input sizes should remain constant.
 * The filter is operated by performing a
 * 	predict, observe
 * cycle derived from the covariance_filter
 */
#include "covFlt.h"

/* Filter namespace */
namespace Bayesian_filter
{

class Iterated_covariance_filter : public Linrz_filter
{
public:
	Iterated_covariance_filter (FM::Subscript x_size, FM::Subscript z_initialsize = 0);
	/* Initialised filter requries an addition iteration limit for the
	   observe algorithm */
	Iterated_covariance_filter& operator= (const Iterated_covariance_filter&);
	// Optimise copy assignment to only copy filter state

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
	FM::Subscript last_z_size;
	void observe_size (FM::Subscript z_size);
							// Permenantly allocated temps
	FM::Vec s;
	FM::Matrix HxT;
};


}//namespace
#endif
