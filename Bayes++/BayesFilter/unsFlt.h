#ifndef _BAYES_FILTER_UNSCENTED
#define _BAYES_FILTER_UNSCENTED

/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
 * $Header$
 * $NoKeywords: $
 */

/*							
 * Unscented Filter.
 *	A Julier-Uhlmann unscented non-linear Kalman filter as an Abstract class
 *
 * The innovation update algorithm is used.
 * Predictions of state and state covariance (and observation) use
 * unscented transformations to interpolate the non-linear predict and observe
 * models. unscented transforms can be further optimised by vary the Kappa
 * parameter from its usual value of 1.
 * Discontinous observe models require that a normailisation function.
 *
 * The predict model is represtented by the state prediction function and a 
 * seperate prediction noise matrix.
 * The observe model is represtented by the observation prediction function and
 * a function to normalise observeations.
 * Derived filters must supply the predict functions and observe function. The
 * observer normailise function is only required for discontinous functions.
 *
 * The filter is operated by performing a
 * 	predict, observe
 * cycle derived from the bayes_filter
 */
#include "bayesFlt.h"

/* Filter namespace */
namespace Bayesian_filter
{

class Unscented_predict_model : public Predict_model_base
/* Specific Unscented prediction model for Addative noise
 *  x(k|k-1) = f(x(k-1|k-1)) + w(x(k))
 *
 * Unscented filter requires
 *  f the function part of the non-linear model
 *  Q the covariance of the addative w(x(k)), w is specificly allow to be a function of state
 */
{
public:
	Unscented_predict_model (FM::Subscript q_size)
	{
		q_unscented = q_size;
	}

	virtual const FM::Vec& f(const FM::Vec& x) const = 0;
	// Functional part of addative model
	// Note: Reference return value as a speed optimisation, MUST be copied by caller.

	virtual const FM::SymMatrix& Q(const FM::Vec& x) const = 0;
	// Covariance of addative noise
	// Note: Reference return value as a speed optimisation, MUST be copied by caller.
private:
	friend class Unscented_filter;	// Filter implementation need to know noise size
	FM::Subscript q_unscented;
};


class Unscented_filter : public Extended_filter, public Functional_filter
{
private:	// TODO: make XX public, requires specification of postcond on XX
	FM::Subscript q_max;	// Maxiumum size allocated for noise model, constructed before XX
	FM::ColMatrix XX;		// Unscented form of state
public:

	Unscented_filter (FM::Subscript x_size, FM::Subscript z_initialsize = 0);
	Unscented_filter& operator= (const Unscented_filter&);
	// Optimise copy assignment to only copy filter state

	void init ();
	void update ();

	void predict (Unscented_predict_model& f);
	// Efficient Unscented prediction 
	
	void predict (Functional_predict_model& f);
	void predict (Addative_predict_model& f);
	Float predict (Linrz_predict_model& f)
	{	// Adapt to use the more general addative model
		predict(static_cast<Addative_predict_model&>(f));
		return 1.;		// Always well condition for addative predict
	}
	
	Float observe (Uncorrelated_addative_observe_model& h, const FM::Vec& z);
	Float observe (Correlated_addative_observe_model& h, const FM::Vec& z);
	Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z)
	{	// Adapt to use the more general addative model
		return observe (static_cast<Uncorrelated_addative_observe_model&>(h),z);
	}
	Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z)
	{	// Adapt to use the more general addative model
		return observe (static_cast<Correlated_addative_observe_model&>(h),z);
	}
	Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s)
	{	// Adapt to use the more general addative model
		return observe (static_cast<Uncorrelated_addative_observe_model&>(h),s);
	}
	Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s)
	{	// Adapt to use the more general addative model
		return observe (static_cast<Correlated_addative_observe_model&>(h),s);
	}

public:						// Exposed Numerical Results
	FM::Vec s;					// Innovation
	FM::SymMatrix S, SI;		// Innovation Covariance and Inverse

protected:
	virtual Float predict_Kappa (unsigned size) const;
	virtual Float observe_Kappa (unsigned size) const;
	/* unscented Kappa values
	   default uses the rule to minimise mean squared error of 4th order term
	*/

protected:					// allow fast operation if z_size remains constant
	FM::Subscript last_z_size;
	void observe_size (FM::Subscript z_size);

private:
	void unscented (FM::ColMatrix& XX, const FM::Vec& x, const FM::SymMatrix& X, Float Scale);
	/* Determine Unscented points for a distribution */
	FM::Subscript x_size;
	FM::Subscript XX_size;	// 2*x_size+1

protected:			   		// Permenantly allocated temps
	FM::ColMatrix fXX;
};


}//namespace
#endif
