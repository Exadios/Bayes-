#ifndef _BAYES_FILTER_UD
#define _BAYES_FILTER_UD
/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Id$
 */

/*
 * UdU' Factorisation of covariance Filter Scheme.
 *  Implementation of a 'Square-root' linearised kalman filter
 * 
 * Bierman's UD factorisatised update algorithm using Agee-Turner UdU' factorisation rank 1 update
 * Thornton's MWG-S factorisation predict algorithm
 * References
 * [1] "Factorisation Methods for Discrete Sequential Estimation" Gerald J. Bierman ISBN 0-12-097350-2
 * [2] "Kalman Filtering, Theory and Practice", Mohinder S. Grewal, Angus P. Andrews ISBN 0-13-211335-X
 *
 * The filter is operated by performing a
 *  predict, observe
 * cycle defined by the base class
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{


class UD_bscheme : virtual public Kalman_state
{
private:
	std::size_t q_max;	// Maxiumum size allocated for noise model, constructed before UD
public:
	FM::Matrix UD;	// UDU factorisation of X with D on diagonal. Lower triangle used as workspace

	struct Predict_byproduct
	{
		Predict_byproduct (std::size_t x_size, std::size_t q_size);
		FM::Vec d, dv, v;
	};

	class Sequential_observe_model : public Parametised_observe_model
	{
	public:
		Sequential_observe_model (std::size_t x_size, std::size_t z_size) :
			Parametised_observe_model(z_size)
		{}
		virtual const FM::Vec& ho (const FM::Vec& x, const std::size_t o, Float& Zv_o, FM::Vec& Hx_o) = 0;
		/* Supplied model (h) for observation using state x, z allows normalisation and model variation
		   Fast model of a single element (o) in observation model
		   Precondition: Hx_o is conformantly dimensioned
		   Postcondition:
			z(k|k-1) = h(x(k|k-1)
			Zv_o(x(k|k-1) = observe noise variance (element o)
			Hx_o(x(k|k-1) = jacobian of h with respect to state x (row o)
		*/
	};

	struct Observe_innovation_byproduct
	// Kalman gain and associated innovation and variance from sequential observe
	{
		Observe_innovation_byproduct (std::size_t x_size, std::size_t z_size);
		FM::Vec s;
		FM::Vec Sv;
		FM::Matrix W;
	};

	struct Observe_byproduct
	{	// Numerical byproducts of observe using observeUD
		Observe_byproduct (std::size_t x_size, std::size_t z_size);
		FM::Vec a;				// observe UD temporary
		FM::Vec w;
		FM::Vec znorm;		// normalised innovation
	};

	struct Observe_linear_byproduct : public Observe_byproduct
	{	// Numerical byproducts of observe with Linear_correlated_observe_model
		Observe_linear_byproduct (std::size_t x_size, std::size_t z_size);
		FM::Vec zpdecol;	// decorrelated zp
		FM::Matrix Gz;			// Z coupling
		FM::Matrix GIHx;		// modified model for linear decorrelation
	};

	UD_bscheme (std::size_t x_size, std::size_t q_maxsize);
	UD_bscheme& operator= (const UD_bscheme&);
	// Optimise copy assignment to only copy filter state

	void init ();
	void update ();

	Float bypredict (Linrz_predict_model& f, Predict_byproduct& b);
	// Predict with explict byproduct

	Float byobserve (Linrz_uncorrelated_observe_model& h, const FM::Vec& z,
			Observe_innovation_byproduct& g, Observe_byproduct& b);
	// Observe with explict byproduct
	Float byobserve (Linear_correlated_observe_model& h, const FM::Vec& z,
			Observe_innovation_byproduct& g, Observe_linear_byproduct& b);
	// Special Linear observe for correlated Z, fast Z decorrelation
	Float byobserve (Sequential_observe_model& h, const FM::Vec& z,
			Observe_innovation_byproduct& g, Observe_byproduct& b);
	// Special Linrz observe using fast sequential model

protected:
	Float predictGq (const FM::Matrix& Fx, const FM::Matrix& G, const FM::Vec& q, Predict_byproduct& b);
	Float observeUD (const Float r, FM::Vec& a, FM::Vec& b, Float& alpha);
};


class UD_scheme : public Linrz_kalman_filter, public UD_bscheme
{
public:
	UD_scheme (std::size_t x_size, std::size_t q_maxsize);
	
	Float predict (Linrz_predict_model& f);
	// Linrz_kalman_filter predict

	Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z);
	// Linrz_kalman_filter observe
	Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z);
	// NO solution for Correlated noise and Linearised model
};

}//namespace
#endif
