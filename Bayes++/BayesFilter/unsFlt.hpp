#ifndef _BAYES_FILTER_UNSCENTED
#define _BAYES_FILTER_UNSCENTED

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens
 * See accompanying Bayes++.htm for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/** Unscented Filter Scheme.
 *
 *  A Julier-Uhlmann Unscented non-linear Kalman filter
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

/** Specific Unscented prediction model for Addative noise.
 * 
 * Prediction equation x(k|k-1) = f(x(k-1|k-1)) + w(x(k))
 *
 * Unscented filter requires
 *  f the function part of the non-linear model
 *  Q the covariance of the addative w(x(k)), w is specificly allow to be a function of state
 */
class Unscented_predict_model : public Predict_model_base
{
public:
	Unscented_predict_model (size_t q_size)
	/// Construct with fixed noise size
	{
		q_unscented = q_size;
	}

	virtual const FM::Vec& f(const FM::Vec& x) const = 0;
	/// Functional part of addative model
	//  Note: Reference return value as a speed optimisation, MUST be copied by caller.

	virtual const FM::SymMatrix& Q(const FM::Vec& x) const = 0;
	/// Covariance of addative noise
	//  Note: Reference return value as a speed optimisation, MUST be copied by caller.
private:
	friend class Unscented_filter;	// Filter implementation need to know noise size
	size_t q_unscented;
};


/** A Julier-Uhlmann Unscented non-linear Kalman filter.
 * 
 * The Unscented transform is used for non-linear state and observation predictions.
 * Uses the original Duplex Unscented transform implementation
 *
 * Observations are fused using innovation gain equations from a Covariance filter.
 *
 * Predictions of state and state covariance (and observation) use
 * unscented transformations to interpolate the non-linear predict and observe
 * models. The Unscented transforms can be optimised by varying the Kappa
 * parameter.
 * Discontinous observe models require a normailisation function.
 *
 * The predict model is represtented by the state prediction function and a 
 * seperate prediction noise matrix.
 * The observe model is represtented by the observation prediction function and
 * a function to normalise observations.
 */
class Unscented_scheme : public Linrz_kalman_filter, public Functional_filter
{
private:
	size_t q_max;			// Maxiumum size allocated for noise model, constructed before XX
public:
	FM::ColMatrix XX;		///< Unscented form of state, with associated kappa
	Float kappa;

	Unscented_scheme (size_t x_size, size_t z_initialsize = 0);
	Unscented_scheme& operator= (const Unscented_scheme&);

 	void init ();
	void init_XX ();
	void update ();
	void update_XX (Float kappa);

	void predict (Unscented_predict_model& f);
	///< Efficient Unscented prediction 
	void predict (Functional_predict_model& f);
	void predict (Addative_predict_model& f);
	Float predict (Linrz_predict_model& f)
	{	// Adapt to use the more general addative model
		predict(static_cast<Addative_predict_model&>(f));
		return 1.;		// Always well condition for addative predict
	}
	
	Float observe (Uncorrelated_addative_observe_model& h, const FM::Vec& z);
	Float observe (Correlated_addative_observe_model& h, const FM::Vec& z);
	///< Unscented filter implements general addative observe models 
	
	Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z)
	{	///< Adapt to use the more general addative model
		return observe (static_cast<Uncorrelated_addative_observe_model&>(h),z);
	}
	Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z)
	{	///< Adapt to use the more general addative model
		return observe (static_cast<Correlated_addative_observe_model&>(h),z);
	}

public:						// Exposed Numerical Results
	FM::Vec s;					///< Last observation innovation
	FM::SymMatrix S, SI;		///< Innovation Covariance and its Inverse

protected:
	//@{
	/** Unscented transform Kappa values.
	    Defaults use the rule which minimise mean squared error of 4th order term
	*/
	virtual Float predict_Kappa (size_t size) const;
	virtual Float observe_Kappa (size_t size) const;
	//@}

protected:					// allow fast operation if z_size remains constant
	size_t last_z_size;
	void observe_size (size_t z_size);

private:
	void unscented (FM::ColMatrix& XX, const FM::Vec& x, const FM::SymMatrix& X, Float scale);
	// Determine Unscented points for a distribution
	size_t x_size;
	size_t XX_size;	// 2*x_size+1

protected:			   		// Permenantly allocated temps
	FM::ColMatrix fXX;
};


}//namespace
#endif
