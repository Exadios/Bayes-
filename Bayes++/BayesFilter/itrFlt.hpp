#ifndef _BAYES_FILTER_ITERATED_COVARIANCE
#define _BAYES_FILTER_ITERATED_COVARIANCE

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Iterated Covariance Filter Scheme.
 *  A non-linear Covariance (Kalman) filter with relinearisation and iteration
 *
 * The observe algorithm uses the iterated non-linear formulation 
 * from Bar-Shalom and Fortmann p.119 (full scheme)
 * Discontinous observe models require that state is normailised with
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

class Iterated_covariance_scheme;

class Iterated_observe_model : virtual public Observe_model_base
/* Linrz observation model which can be iterated
    Hx can be relinearised
 */
{
protected: // model is not sufficient, it is used to build observe model's
	Iterated_observe_model ()
	{}
public:
	virtual void relinearise (const FM::Vec& x) = 0;
	// Relinearised about state x
};


class Iterated_terminator : public Bayes_base
/*
 * Termination condition for filter Iteration
 *  Used by iterated observe to parameterise termination condition
 *  If iteration continues, the terminator must also relinearise the model about the filters new state
 *
 * Defaults to immediately terminating the iteration
 *
 * A more useful terminator can built by derivation.
 * For example terminator constructed with a reference to the filter and model can
 * detect convergence of x and/or X
 */
{
public:
	virtual bool term_or_relinearize (const Iterated_covariance_scheme& f)
	{
		return true;
	}
};

class Counted_iterated_terminator : public Iterated_terminator
/*
 * Termination condition with a simple fixed number of iterations 
 */
{
public:
	Counted_iterated_terminator (Iterated_observe_model& model, unsigned iterations) :
		m(model), i(iterations)
	{}
	bool term_or_relinearize (const Iterated_covariance_scheme& f);
	Iterated_observe_model& m;
	unsigned i;
};



class Iterated_covariance_scheme : public Linrz_kalman_filter
{
public:
	Iterated_covariance_scheme (std::size_t x_size);
	/* Initialised filter requries an addition iteration limit for the
	   observe algorithm */
	Iterated_covariance_scheme& operator= (const Iterated_covariance_scheme&);
	// Optimise copy assignment to only copy filter state

	void init ();
	void update ();

	Float predict (Linrz_predict_model& f);
	// Linrz_kalman_filter predict
	Float predict (Gaussian_predict_model& f);
	// Specialised 'stationary' predict, only addative noise

	Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z)
	{	// Linrz_kalman_filter observe with default termination
		Iterated_terminator term;
		const std::size_t z_size = h.Hx.size1();
		FM::Vec s(z_size);
		FM::SymMatrix S(z_size,z_size), SI(z_size, z_size);
		FM::Matrix W(h.Hx.size2(), z_size);
		return byobserve (h, term, z, s,S,SI,W);
	}
	Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z)
	{	// Linrz_kalman_filter observe with default termination
		Iterated_terminator term;
		const std::size_t z_size = h.Hx.size1();
		FM::Vec s(z_size);
		FM::SymMatrix S(z_size,z_size), SI(z_size, z_size);
		FM::Matrix W(h.Hx.size2(), z_size);
		return byobserve (h, term, z, s,S,SI,W);
	}
	Float byobserve (Linrz_uncorrelated_observe_model& h, Iterated_terminator& term, const FM::Vec& z,
				FM::Vec& s, FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W);
	Float byobserve (Linrz_correlated_observe_model& h, Iterated_terminator& term, const FM::Vec& z,
				FM::Vec& s, FM::SymMatrix& S, FM::SymMatrix& SI, FM::Matrix& W);
	// Observe with iteration and explict byproduct

protected:			   		// Permenantly allocated temps
	FM::RowMatrix tempX;
};


}//namespace
#endif
