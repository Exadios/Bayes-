#ifndef _BAYES_FILTER__MODELS
#define _BAYES_FILTER__MODELS

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Predict and Observe models
 *  These models extend, adapt and simpilify the fundamental Bayesian filter models
 *  Simple : Simplify model construction and use
 *  General: Combine model types to form more general models
 *  Adapted: Adapt one model type into another
 */
#include <boost/function.hpp>

/* Filter namespace */
namespace Bayesian_filter
{

typedef boost::function1<const FM::Vec&, const FM::Vec&> State_function;
// A generalised function of state. Compatible with predict and observe models


class Simple_addative_predict_model : public Addative_predict_model
// Addative predict model initialised from function and model matricies
{
	State_function ff;
public:
	Simple_addative_predict_model (State_function f_init, const FM::Matrix& G_init, const FM::Vec& q_init);
	// Precondition: G, q are conformantly dimensioned (not checked)

	// No default assignment operator

	virtual const FM::Vec& f(const FM::Vec& x) const
	{	return ff(x);
	}
};

class Simple_linrz_predict_model : public Linrz_predict_model
// Linrz predict model initialised from function and model matricies
{
	State_function ff;
public:
	Simple_linrz_predict_model (State_function f_init, const FM::Matrix& Fx_init, const FM::Matrix& G_init, const FM::Vec& q_init);
	// Precondition: Fx, G, q are conformantly dimensioned (not checked)

	// No default assignment operator

	virtual const FM::Vec& f(const FM::Vec& x) const
	{	return ff(x);
	}
};

class Simple_linear_predict_model : public Linear_predict_model
// Linear predict model initialised from model matricies
{
public:
	Simple_linear_predict_model (const FM::Matrix& Fx_init, const FM::Matrix& G_init, const FM::Vec& q_init);
	// Precondition: Fx, q and G are conformantly dimensioned (not checked)
};


class Simple_linrz_correlated_observe_model : public Linrz_correlated_observe_model
// Linrz observe model initialised from function and model matricies
{
	State_function ff;
public:
	Simple_linrz_correlated_observe_model (State_function f_init, const FM::Matrix& Hx_init, const FM::SymMatrix& Z_init);
	// Precondition: Hx, Z are conformantly dimensioned (not checked)
	// No default assignment operator

	virtual const FM::Vec& h(const FM::Vec& x) const
	{	return ff(x);
	}
};

class Simple_linrz_uncorrelated_observe_model : public Linrz_uncorrelated_observe_model
// Linrz observe model initialised from function and model matricies
{
	State_function ff;
public:
	Simple_linrz_uncorrelated_observe_model (State_function f_init, const FM::Matrix& Hx_init, const FM::Vec& Zv_init);
	// Precondition: Hx, Zv are conformantly dimensioned (not checked)
	// No default assignment operator

	virtual const FM::Vec& h(const FM::Vec& x) const
	{	return ff(x);
	}
};

class Simple_linear_correlated_observe_model : public Linear_correlated_observe_model
// Linear observe model initialised from model matricies
{
public:
	Simple_linear_correlated_observe_model (const FM::Matrix& Hx_init, const FM::SymMatrix& Z_init);
	// Precondition: Hx, Z are conformantly dimensioned (not checked)
};

class Simple_linear_uncorrelated_observe_model : public Linear_uncorrelated_observe_model
// Linear observe model initialised from model matricies
{
public:
	Simple_linear_uncorrelated_observe_model (const FM::Matrix& Hx_init, const FM::Vec& Zv_init);
	// Precondition: Hx, Zv are conformantly dimensioned (not checked)
};



/*
 * Generalise addative observe model to likelihood observe model
 */

// General Linearised Uncorrelated Addative
class General_LzUnAd_observe_model: public Linrz_uncorrelated_observe_model, public Likelihood_observe_model
{
public:
	General_LzUnAd_observe_model (size_t x_size, size_t z_size) :
		Linrz_uncorrelated_observe_model(x_size, z_size),
		Likelihood_observe_model(z_size),
		Zv_inv(z_size)
	{
		zset = false;
	}
	virtual Float L(const FM::Vec& x) const;
	// Definition of likelihood for addative noise model given zz
	virtual void Lz (const FM::Vec& zz);
	// Fix the observation zz about which to evaluate the Likelihood function
	// Zv is also fixed

private:
	FM::Vec Zv_inv;			// Inverse Noise Covariance given zz
	Float logdetZ;			// log(det(Z))
	mutable bool zset;	
};

// General Linear Uncorrelated Addative
class General_LiUnAd_observe_model: public General_LzUnAd_observe_model
{
public:
	General_LiUnAd_observe_model (size_t x_size, size_t z_size) :
		General_LzUnAd_observe_model(x_size, z_size),
		hx(z_size)
	{}
	const FM::Vec& h(const FM::Vec& x) const
	{	// Provide a linear implementation of functional h assumes model is already Linrz for Hx
		hx.assign (FM::prod(Hx,x));
		return hx;
	}
private:
	mutable FM::Vec hx;
};

// General Linearised Correlated Addative
class General_LzCoAd_observe_model: public Linrz_correlated_observe_model, public Likelihood_observe_model
{
public:
	General_LzCoAd_observe_model (size_t x_size, size_t z_size) :
		Linrz_correlated_observe_model(x_size, z_size),
		Likelihood_observe_model(z_size),
		Z_inv(z_size,z_size)
	{
		zset = false;
	}
	virtual Float L(const FM::Vec& x) const;
	// Definition of likelihood for addative noise model given zz
	virtual void Lz (const FM::Vec& zz);
	// Fix the observation zz about which to evaluate the Likelihood function
	// Z is also fixed

private:
	FM::SymMatrix Z_inv;	// Inverse Noise Covariance
	Float logdetZ;			// log(det(Z)
	static Float scaled_vector_square(const FM::Vec& v, const FM::SymMatrix& V);
	mutable bool zset;	
};

// General Linear Correlated Addative
class General_LiCoAd_observe_model: public General_LzCoAd_observe_model
{
public:
	General_LiCoAd_observe_model (size_t x_size, size_t z_size) :
		General_LzCoAd_observe_model(x_size, z_size),
		hx(z_size)
	{}
	const FM::Vec& h(const FM::Vec& x) const
	{	// Provide a linear implementation of functional h assumes model is already Linrz for Hx
		hx.assign (FM::prod(Hx,x));
		return hx;
	}
private:
	mutable FM::Vec hx;
};


/*
 * Model Adaptors:
 */

class Adapted_Correlated_addative_observe_model : public Correlated_addative_observe_model
/*
 * Adapt Uncorrelated_addative_observe_model to an equivilent
 * Correlated_addative_observe_model_adaptor
 */
{
public:
	Adapted_Correlated_addative_observe_model (Uncorrelated_addative_observe_model& adapt);
	const FM::Vec& h(const FM::Vec& x) const
	{
		return unc.h(x);
	}
	inline void normalise (FM::Vec& z_denorm, const FM::Vec& z_from) const
	{
		unc.normalise (z_denorm, z_from);
	};
private:
	Uncorrelated_addative_observe_model& unc;
};

class Adapted_Linrz_correlated_observe_model : public Linrz_correlated_observe_model
/*
 * Adapt Linrz_uncorrelated_observe_model to an equivilent
 * Linrz_correlated_observe_model
 */
{
public:
	Adapted_Linrz_correlated_observe_model (Linrz_uncorrelated_observe_model& adapt);
	const FM::Vec& h(const FM::Vec& x) const
	{
		return unc.h(x);
	}
	inline void normalise (FM::Vec& z_denorm, const FM::Vec& z_from) const
	{
		unc.normalise (z_denorm, z_from);
	};
protected:
	Linrz_uncorrelated_observe_model& unc;
};

}// namespace

#endif
