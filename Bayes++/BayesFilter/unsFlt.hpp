#ifndef _BAYES_FILTER_UNSCENTED
#define _BAYES_FILTER_UNSCENTED

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Id$
 */

/*
 * Unscented Filter Scheme.
 *
 *  A Julier-Uhlmann Unscented non-linear Kalman filter
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

/** Specific Unscented predict model for Additive noise.
 * 
 * predict equation x(k|k-1) = f(x(k-1|k-1)) + w(x(k))
 *
 * Unscented filter requires
 *  f the function part of the non-linear model
 *  Q the covariance of the additive w(x(k)), w is specificly allow to be a function of state
 */
class Unscented_predict_model : public Predict_model_base
{
public:
	Unscented_predict_model ()
	{}

	virtual const FM::Vec& f(const FM::Vec& x) const = 0;
	/// Functional part of additive model
	//  Note: Reference return value as a speed optimisation, MUST be copied by caller.

	virtual const FM::SymMatrix& Q(const FM::Vec& x) const = 0;
	/// Covariance of additive noise
	//  Note: Reference return value as a speed optimisation, MUST be copied by caller.
};


/** A Julier-Uhlmann Unscented non-linear Kalman filter.
 * 
 * The Unscented transform is used for non-linear state and observation predicts.
 * Uses the original Duplex Unscented transform implementation
 *
 * Observations are fused using innovation gain equations from a Covariance filter.
 *
 * Predicts of state and state covariance (and observation) use
 * unscented transformations to interpolate the non-linear predict and observe
 * models. The Unscented transforms can be optimised by varying the Kappa
 * parameter.
 * Discontinous observe models require a normailisation function.
 *
 * The predict model is represtented by the state predict function and a 
 * seperate predict noise matrix.
 * The observe model is represtented by the observation predict function and
 * a function to normalise observations.
 */
class Unscented_scheme : public Linrz_kalman_filter, public Functional_filter
{
private:
	std::size_t q_max;			// Maxiumum size allocated for noise model, constructed before XX
public:
	FM::ColMatrix XX;		///< Unscented form of state, with associated kappa
	Float kappa;

	Unscented_scheme (std::size_t x_size);
	Unscented_scheme& operator= (const Unscented_scheme&);

 	void init ();
	void init_XX ();
	void update ();
	void update_XX (Float kappa);

	void predict (Unscented_predict_model& f);
	///< Efficient Unscented predict 
	void predict (Functional_predict_model& f);
	void predict (Additive_predict_model& f);
	Float predict (Linrz_predict_model& f)
	{	///< Linrz_kalman_filter predict
		predict(static_cast<Additive_predict_model&>(f));
		return 1.;		// Always well condition for additive predict
	}
	
	Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z)
	{	///< Linrz_kalman_filter observe
		const std::size_t z_size = h.Hx.size1();
		FM::Vec s(z_size);
		FM::SymMatrix S(z_size,z_size), SI(z_size,z_size);
		FM::Matrix W(h.Hx.size2(), z_size);
		return byobserve (static_cast<Uncorrelated_additive_observe_model&>(h), z, s,S,SI,W);
	}
	Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z)
	{	///< Linrz_kalman_filter observe
		const std::size_t z_size = h.Hx.size1();
		FM::Vec s(z_size);
		FM::SymMatrix S(z_size,z_size), SI(z_size,z_size);
		FM::Matrix W(h.Hx.size2(), z_size);
		return byobserve (static_cast<Correlated_additive_observe_model&>(h), z, s,S,SI,W);
	}
	
	Float byobserve (Uncorrelated_additive_observe_model& h, const FM::Vec& z,
				FM::Vec& s, FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W);
	Float byobserve (Correlated_additive_observe_model& h, const FM::Vec& z,
				FM::Vec& s, FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W);
	///< General additive observe models form, with explict byproduct
	

protected:
	//@{
	/** Unscented transform Kappa values.
	    Defaults use the rule which minimise mean squared error of 4th order term
	*/
	virtual Float predict_Kappa (std::size_t size) const;
	virtual Float observe_Kappa (std::size_t size) const;
	//@}

private:
	void unscented (FM::ColMatrix& XX, const FM::Vec& x, const FM::SymMatrix& X, Float scale);
	// Determine Unscented points for a distribution
	std::size_t x_size;
	std::size_t XX_size;	// 2*x_size+1

protected:			   		// Permenantly allocated temps
	FM::ColMatrix fXX;
};


}//namespace
#endif
