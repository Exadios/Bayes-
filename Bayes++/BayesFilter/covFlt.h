#ifndef _BAYES_FILTER_COVARIANCE
#define _BAYES_FILTER_COVARIANCE

/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
 * $Header$
 * $NoKeywords: $
 */

/*
 * Covariance Filter.
 *	A Covariance (Extended Kalman) filter as an Abstract class
 *
 * To work with with Linear and Linrz models
 *  a) a state seperate from covariance prediction is used.
 *  b) a EKF innovation update algorithm is used.
 * Discontinous observe models require that prediction is normailised with
 * respect to the observation.
 *
 * Derived filters must supply predict and observe model functions.
 * State and control input sizes should remain constant.
 * A initial observation size may also be specified for efficiency.
 * 
 * The filter is operated by performing a
 * 	predict, observe
 * cycle derived from the Linrz_filter
 */
#include "bayesFlt.h"

/* Filter namespace */
namespace Bayesian_filter
{

class Covariance_filter : public Extended_filter
{
public:
	Covariance_filter (FM::Subscript x_size, FM::Subscript z_initialsize = 0);
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
	FM::Subscript last_z_size;
	void observe_size (FM::Subscript z_size);
};


}//namespace
#endif
