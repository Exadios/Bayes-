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
 * Predict and Observe model
 *  These model extend, adapt and simpilify the fundematale Bayes filter models
 *  Simple : Simplify model construction and use
 *  General: Combine model types to form more general models
 *  Adapted: Adapt one model type into another
 */


/* Filter namespace */
namespace Bayesian_filter
{

class Simple_addative_predict_model : public Addative_predict_model
// Addative predict model initialised from function and model matricies
{
	Predict_function& ff;
public:
	Simple_addative_predict_model (Predict_function& f_init, const FM::Matrix& G_init, const FM::Vec& q_init);
	// Precondition: G, q are conformantly dimensioned (not checked)
	// No default assignment operator

	virtual const FM::Vec& f(const FM::Vec& x) const
	{	return ff.fx(x);
	}
};

class Simple_linrz_predict_model : public Linrz_predict_model
// Linrz predict model initialised from function and model matricies
{
	Predict_function& ff;
public:
	Simple_linrz_predict_model (Predict_function& f_init, const FM::Matrix& Fx_init, const FM::Matrix& G_init, const FM::Vec& q_init);
	// Precondition: Fx, G, q are conformantly dimensioned (not checked)
	// No default assignment operator

	virtual const FM::Vec& f(const FM::Vec& x) const
	{	return ff.fx(x);
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
	Observe_function &ff;
public:
	Simple_linrz_correlated_observe_model (Observe_function& f_init, const FM::Matrix& Hx_init, const FM::SymMatrix& Z_init);
	// Precondition: Hx, Z are conformantly dimensioned (not checked)
	// No default assignment operator

	virtual const FM::Vec& h(const FM::Vec& x) const
	{	return ff.h(x);
	}
};

class Simple_linrz_uncorrelated_observe_model : public Linrz_uncorrelated_observe_model
// Linrz observe model initialised from function and model matricies
{
	Observe_function& ff;
public:
	Simple_linrz_uncorrelated_observe_model (Observe_function& f_init, const FM::Matrix& Hx_init, const FM::Vec& Zv_init);
	// Precondition: Hx, Zv are conformantly dimensioned (not checked)
	// No default assignment operator

	virtual const FM::Vec& h(const FM::Vec& x) const
	{	return ff.h(x);
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
 * Generalise a addative predict model to sampled predict model
 *  To instantiate template sqrt is required
 */
struct General_sampled_predict_random
/* Random number generators interface
 *  Helper to allow polymorthic use of random number generators 
 */
{
	virtual void normal(FM::DenseVec& v) = 0;
};

template <class Predict_model>
class General_sampled_predict_model: public Predict_model, public Sampled_predict_model
{
public:
	typedef General_sampled_predict_random Random;

	General_sampled_predict_model (size_t x_size, size_t q_size, Random& random_helper) :
		Predict_model(x_size, q_size),
		Sampled_predict_model(),
		genn(random_helper),
		xp(x_size),
		n(q_size), rootq(q_size)
	{
		first_init = true;
	}

	virtual const FM::Vec& fx(const FM::Vec& x) const
	/*
	 * Definition of sampler for addative noise model given state x
	 *  Generate Gaussian correlated samples
	 * Precond: init_GqG, automatic on first use
	 */
	{
		if (first_init)
			init_GqG();
							// Predict state using supplied functional predict model
		xp = f(x);
							// Additive random noise
		genn.normal(n);				// independant zero mean normal
									// multiply elements by std dev
		for (FM::DenseVec::iterator ni = n.begin(); ni != n.end(); ++ni) {
			*ni *= rootq[ni.index()];
		}
		xp += FM::prod(G,n);			// add correlated noise
		return xp;
	}

	void init_GqG() const
	/* initialise predict given a change to q,G
	 *  Implementation: Update rootq
	 */
	{
		first_init = false;
		for (FM::Vec::const_iterator qi = q.begin(); qi != q.end(); ++qi) {
			if (*qi < 0.)
				throw Bayesian_filter::Bayes_filter_exception ("Negative q in init_GqG");
			rootq[qi.index()] = std::sqrt(*qi);
		}
	}
private:
	Random& genn;
	mutable FM::Vec xp;
	mutable FM::DenseVec n;
	mutable FM::Vec rootq;		// Optimisation of sqrt(q) calculation, automatic on first use
	mutable bool first_init;	
};

/*
 * General Model: Combine properties of mutliple models
 *  Names a shortened to first two letters of their model properties
 */
typedef General_sampled_predict_model<Linear_predict_model> General_LiAd_predict_model;
typedef General_sampled_predict_model<Linear_invertable_predict_model> General_LiInAd_predict_model;



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
