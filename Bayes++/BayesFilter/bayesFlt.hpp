#ifndef _BAYES_FILTER
#define _BAYES_FILTER

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

// Common headers required for declerations
#include <exception>
#include "matSupSub.hpp"	// Expect to find the actual matrix support headers elsewhere


/* Filter namespace */
namespace Bayesian_filter
{
	// Allow use of a short name for matrix namespace
	namespace FM = Bayesian_filter_matrix;

/*
 * Abstraction support classes
 */

class Bayes_base {
/*
 * A very abstract base representation!
 * Interface provides: type, internal error handing, and destruction
 */
public:
	typedef Bayesian_filter_matrix::Float Float;
	// Type used thoughout as a number representation for state etc

	virtual ~Bayes_base() = 0;
	// Provide abstract destruction
protected:
	void filter_error (const char* error_description) const;
	// Report a filter error, throw a Bayes_filter_exception
	// No exception saftey rules are specified, assume the object is invalid
};


class Bayes_filter_exception : public std::exception {
/*
 * Bayesian Filter Exception
 *	Base class for all exception produced by filter heirachy
 */
public:
	Bayes_filter_exception (const char* description)
	{	error_description = description;
	};
	const char *what() const throw()
	{	return error_description;
	}
private:
	const char* error_description;
};


class Numerical_rcond : private Bayes_base {
/*
 * Numerical comparison of reciprocal condition numbers
 *  Required for all linear algebra in models and filters
 *  Implements minimum allowable reciprocal condition number for PD Matrix factorisations
 */
public:
	Numerical_rcond()
	{	limit_PD = limit_PD_init;
	}
	void set_limit_PD(Float nl)
	{	limit_PD = nl;
	}
	inline void check_PSD (Float rcond, const char* error_description) const
	/* Checks a the reciprocal condition number
	 * Generates a Bayes_filter_exception if value represents a NON PSD matrix
	 * The rcond != rcond provides a test for IEC 559 NaN values
	 */
	{	if (rcond < 0 || rcond != rcond)
			throw Bayes_filter_exception (error_description);
	}

	inline void check_PD (Float rcond, const char* error_description) const
	/* Checks a reciprocal condition number
	 * Generates a Bayes_filter_exception if value represents a NON PD matrix
	 * I.e. rcond is bellow given conditioning limit
	 * The rcond != rcond provides a test for IEC 559 NaN values
	 */
	{	if (rcond < limit_PD || rcond != rcond)
			throw Bayes_filter_exception (error_description);
	}
private:
	Float limit_PD;		
	const static Float limit_PD_init;	// Initial common value for limit_PD
};


/*
 * Abstract model
 *  A generic function object
 */
class Function_model
{
public:
	virtual const FM::Vec& fx(const FM::Vec& x) const = 0;
	// Note: Reference return value as a speed optimisation, MUST be copied by caller.
};


/*
 * Abstract Prediction models
 *  Prediction model are used to parameterise predict functions if filters
 */
class Predict_model_base : public Bayes_base
{
	// Empty
};

class Sampled_predict_model : virtual public Predict_model_base, public Function_model
/* Sampled stochastic prediction model
    x*(k) = fx(x(k-1), w(k))
   fx should generate samples from the stochastic variable w(k)
   The fundamental model that is used instead of the prediction likelihood function L(x*|x)
   Since drawing samples from an arbitary L is non-trivial (see MCMC theory)
   the burden is place on the model to generate these samples.
   Defines an Interface without data members
 */
{
	// only fx
};

class Functional_predict_model : public Sampled_predict_model
/* Functional (non-stochastic) prediction model f
    x*(k) = fx(x(k-1))
   The fundamental model that is used instead of the prediction likelihood function L(x*|x)
   L is a delta function which isn't much use for numerical systems.
   Conviently the model is a Sampled_predict_model but without any 'w'
   Defines an Interface without data members
 */
{
public:
	const FM::Vec& operator()(const FM::Vec& x) const
	// Operator form of functional model
	{	return fx(x);
	}
};

class Addative_predict_model : virtual public Predict_model_base
/* Addative Gaussian noise prediction model
   This fundamental model for linear/linearised filtering
	x(k|k-1) = f(x(k-1|k-1)) + G(k)w(k)
	q(k) = state noise covariance, q(k) is covariance of w(k)
	G(k) = state noise coupling
*/
{
public:
	Addative_predict_model (size_t x_size, size_t q_size) :
		q(q_size), G(x_size, q_size)
	{}

	FM::Vec q;			// Noise Covariance
	FM::Matrix G;		// Noise Coupling

	virtual const FM::Vec& f(const FM::Vec& x) const = 0;
	// Functional part of addative model
	// Note: Reference return value as a speed optimisation, MUST be copied by caller.

	Numerical_rcond rclimit;
};

class Linrz_predict_model : public Addative_predict_model
/* Linrz prediction model Fx, of functional part of addative model f about state x
	x(k|k-1) = f(x(k-1|k-1)
	Fx(x(k-1|k-1) = Jacobian of f with respect to state x
 */
{
public:
	Linrz_predict_model (size_t x_size, size_t q_size);
	FM::Matrix Fx;		// Model
};

class Linear_predict_model : public Linrz_predict_model
/* Linear prediction model Fx about state x (fixed size)
	x(k|k-1) = Fx(k-1|k-1) * x(k-1|k-1)
 */
{
public:
	Linear_predict_model (size_t x_size, size_t q_size);
	/* Set constant sizes for
		x_size of the state vector
		q_size of the noise vector
	   Postcondition:
		Fx, q and G are conformantly dimensioned
	 */
	const FM::Vec& f(const FM::Vec& x) const
	{	// Provide a linear implementation of functional f assumes model is already Linrz for Fx
		xp.assign (FM::prod(Fx,x));
		return xp;
	}
private:
	mutable FM::Vec xp;
};

class Linear_invertable_predict_model : public Linear_predict_model
/*
 * Linear invertable prediction model
 */
{
public:
	Linear_invertable_predict_model (size_t x_size, size_t q_size);
	// The inverse model: x(k-1|k-1) = f(x(k|k-1) and Fx,q,G are all matrix inverses
	class inverse_model {
	public:
		inverse_model (size_t x_size, size_t q_size);
		FM::Matrix Fx;		// Model
		FM::Vec q;			// Noise Covariance
		FM::Matrix G;		// Noise Coupling
	} inv;
};



/*
 * Observation of state model
 */

class Likelihood_observe_model : virtual public Bayes_base
/* Likelihood observe model L(z |x)
 *  The most fundamental Bayesian definition of an observation
 * Defines an Interface without data members
 */
{
public:
	Likelihood_observe_model(size_t z_size) : z(z_size)
	{}
	virtual Float L(const FM::Vec& x) const = 0;
	// Likelihood L(z | x)

	virtual void Lz (const FM::Vec& zz)
	// Set the observation zz about which to evaluate the Likelihood function
	{
		z = zz;
	}
protected:
	FM::Vec z;			// z set by Lz
};

class Functional_observe_model : virtual public Bayes_base, public Function_model
/* Functional (non-stochastic) observe model h
 *  z(k) = hx(x(k|k-1))
 * This is a seperate fundamental model and not derived from likelihood because
 * L is a delta function which isn't much use for numerical systems
 * Defines an Interface without data members
 */
{
public:
	Functional_observe_model(size_t /*z_size*/)
	{}
	const FM::Vec& operator()(const FM::Vec& x) const
	{	return fx(x);
	}

};

class Parametised_observe_model : virtual public Bayes_base
/* Observation model parametised with a fixed z size
 * Model is assume to have linear components
 */
{
public:
	Parametised_observe_model(size_t /*z_size*/)
	{}
	virtual void normalise (FM::Vec& /*z_denorm*/, const FM::Vec& /*z_from*/) const
	/* Normalise one observation state (z_denorm) from another if observation model is discontinous.
	    Default for continuous h
	 */
	{}
	Numerical_rcond rclimit;
};

class Uncorrelated_addative_observe_model : public Parametised_observe_model
/* Observation model, uncorrelated addative observation noise
	Z(k) = I * Zv(k) observe noise variance vector Zv
 */
{
public:
	Uncorrelated_addative_observe_model (size_t z_size) :
		Parametised_observe_model(z_size), Zv(z_size)
	{}
	FM::Vec Zv;			// Noise Variance
	virtual const FM::Vec& h(const FM::Vec& x) const = 0;
	// Functional part of addative model
	// Note: Reference return value as a speed optimisation, MUST be copied by caller.
};

class Correlated_addative_observe_model : public Parametised_observe_model
/* Observation model, correlated addative observation noise
	Z(k) = observe noise covariance
 */
{
public:
	Correlated_addative_observe_model (size_t z_size) :
		Parametised_observe_model(z_size), Z(z_size,z_size)
	{}
	FM::SymMatrix Z;	// Noise Covariance

	virtual const FM::Vec& h(const FM::Vec& x) const = 0;
	// Functional part of addative model
	// Note: Reference return value as a speed optimisation, MUST be copied by caller.
};

class Jacobian_observe_model
/* Linrz observation model Hx, h about state x (fixed size)
	z(k) = h(x(k-1|k-1)
	Hx(x(k|k-1) = Jacobian of h with respect to state x
 */
{
protected:
	Jacobian_observe_model (size_t x_size, size_t z_size) :
		Hx(z_size, x_size)
	{}
	FM::Matrix Hx;		// Model
};

class Linrz_correlated_observe_model : public Correlated_addative_observe_model, public Jacobian_observe_model
/* Linrz observation model Hx, h with repespect to state x (fixed size)
    correlated observation noise
	z(k) = h(x(k-1|k-1)
	Hx(x(k|k-1) = Jacobian of f with respect to state x
	Z(k) = observe noise covariance
 */
{
public:
	Linrz_correlated_observe_model (size_t x_size, size_t z_size) :
		Correlated_addative_observe_model(z_size), Jacobian_observe_model(x_size, z_size)
	{}
	using Jacobian_observe_model::Hx;
};

class Linrz_uncorrelated_observe_model : public Uncorrelated_addative_observe_model, public Jacobian_observe_model
/* Linrz observation model Hx, h with repespect to state x (fixed size)
    uncorrelated observation noise
	z(k) = h(x(k-1|k-1)
	Hx(x(k|k-1) = Jacobian of f with respect to state x
	Zv(k) = observe noise covariance
 */
{
public:
	Linrz_uncorrelated_observe_model (size_t x_size, size_t z_size) :
		Uncorrelated_addative_observe_model(z_size), Jacobian_observe_model(x_size, z_size)
	{}
	using Jacobian_observe_model::Hx;
};

class Linear_correlated_observe_model : public Linrz_correlated_observe_model
/* Linear observation model, correlated observation noise
	z(k) = Hx(k) * x(k|k-1)
 */
{
public:
	Linear_correlated_observe_model (size_t x_size, size_t z_size) :
		Linrz_correlated_observe_model(x_size, z_size), hx(z_size)
	{}
	const FM::Vec& h(const FM::Vec& x) const
	{	// Provide a linear implementation of functional h assumes model is already Linrz for Hx
		hx.assign (FM::prod(Hx,x));
		return hx;
	}
private:
	mutable FM::Vec hx;
};

class Linear_uncorrelated_observe_model : public Linrz_uncorrelated_observe_model
/* Linear observation model, uncorrelated observation noise
	z(k) = Hx(k) * x(k|k-1)
 */
{
public:
	Linear_uncorrelated_observe_model (size_t x_size, size_t z_size) :
		Linrz_uncorrelated_observe_model(x_size, z_size), hx(z_size)
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
 * Bayesian Filter
 *
 * A filter that uses Bayes rule to fuse the probability distributions
 * of prior and likelhood 
 */
class Bayes_filter_base : public Bayes_base
{
	// Empty
};

/*
 * Likelihood Filter - Abstract filtering property
 * Represents only the Bayesian Likelihood of a state observation
 */
class Likelihood_filter : virtual public Bayes_filter_base
{
public:
	/* Virtual functions for filter algorithm */

	virtual void observe (Likelihood_observe_model& h, const FM::Vec& z) = 0;
	/* Observation state posterior using likelihood model h at z
	*/
};

/*
 * Functional Filter - Abstract filtering property
 * Represents only filter prediction by a simple functional
 * (non-stochastic) model
 */
class Functional_filter : virtual public Bayes_filter_base
{
public:
	/* Virtual functions for filter algorithm */

	virtual void predict (Functional_predict_model& f) = 0;
	/* Predict state with functional no noise model
		Requires x(k|k), X(k|k) or internal equivilent
		Predicts x(k+1|k), X(k+1|k), using predict model
	*/
private:
	virtual void observe (Functional_observe_model& /*h*/, const FM::Vec& /*z*/)
	/* Observation z(k) with functional no noise model
		Requires x(k|k), X(k|k)
	   PREVENTS use of this model as h needs to be invertable!
	*/
	{}
};

/*
 * State Filter - Abstract filtering property
 * Represents only filter state and an update on that state
 */
class State_filter : virtual public Bayes_filter_base
{
public:
	State_filter (size_t x_size);
	/* Set constant sizes, state must not be empty (must be >=1)
		Exceptions:
			bayes_filter_exception is x_size < 1
	 */

	FM::Vec x;			// expected state

	/* Virtual functions for filter algorithm */

	virtual void update () = 0;
	/* Update filters state
	     Updates x(k|k)
	*/
};


/*
 * Kalman Filter: Linear filter for 1st (mean) and 2nd (covariance) moments of a distribution
 *
 * Probability distributions are represted by state vector (x) and a covariance matix.(X)
 *
 * State (x) sizes is assumed to remain constant.
 * The filter is operated by performing a
 * 	predict*, observe* cycle. (* -> 0 or more)
 * The state and state covariance are public so they can be directly manipulated.
 *  init: Should be called if x or X are altered
 *  update: Guarantees that any internal changes made by predict or observe are
 *	reflected in x,X. This allows considerable flexibility so filter implemtations can
 *  use different numerical representations
 *
 * Because Kalman filters can use a variety of noise models we do not specialise the functional
 * (no noise) form.
 *
 * Derived filters supply definititions for the abstract functions and determine the algorithm used
 * to implement the filter.
 */

class Kalman_filter : public State_filter
{
public:
	FM::SymMatrix X;	// state covariance

	Float rcond_limit;	// Minimum allowable reciprocal condition number for PD Matrix factorisations
						// Applies to state covariance or its derived matrices

	Kalman_filter (size_t x_size);
	/* Initialise filter and set constant sizes
	 */

	/* Virtual functions for filter algorithm */

	virtual void init () = 0;
	/* Initialise from current state and state covariance
	     Requires x(k|k), X(k|k)
	*/
	void init_kalman (const FM::Vec& x, const FM::SymMatrix& X)
	/* Initialise from a state and state covariance
	     Requires x(k|k), X(k|k)
		 Parameters that reference the instance's x and X members is valid
	*/
	{
		Kalman_filter::x = x;
		Kalman_filter::X = X;
		init();
	}
	virtual void update () = 0;
	/* Update filters state and state covariance 
	     Updates x(k|k), X(k|k)
	*/

	Numerical_rcond rclimit;	// Applies to ALL covariance matrices in algorithms
};


/*
 * Generic Kalman type filter for all linearizable model
 *  Linrz == Linearizable (or Extended Kalman Filter)
 *	A linear, or gradient Linearized filter as an Abstract class
 *
 * State (x) size is fixed.
 * Predict uses a Linrz_predict_model that maintains a Jacobian matrix Fx
 * and a prediction covariance Q
 * NOTE: Functional (non-stochastic) predict is NOT possible as predict requires Fx.
 *
 * Observe uses a Linrz_observe_model and a variable size observation (z)
 * There are two variants for correlated and uncorrelated observation noise
 * Derived filters supply the init,predict,observe,update functions and determine
 * the algorithm used to implement the filter.
 */

class Linrz_filter : public Kalman_filter
{ 
public:
	Linrz_filter (size_t x_size);

	/* Virtual functions for filter algorithm */

	virtual Float predict (Linrz_predict_model& f) = 0;
	/* Predict state using a Linrz model
		Requires x(k|k), X(k|k) or internal equivilent
		Returns: Reciprocal condition number of primary matrix used in predict computation (1. if none)
	*/

	virtual Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z) = 0;
	virtual Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z) = 0;
	/* Observation z(k) and with (Un)correlated observation noise model
		Requires x(k|k), X(k|k) or internal equivilent
		Returns: Reciprocal condition number of primary matrix used in observe computation (1. if none)
	*/
};


/*
 * Extended Kalman type filter for all linearizable model with innovation observations
 *
 * As Linrz_filter
 *
 * Observe is implemented using an innovation computed from the non-linear part of the
 * obseve model and linear part of the Linrz_observe_model
 *
 */

class Extended_filter : public Linrz_filter
{ 
public:
	Extended_filter (size_t x_size);

	/* Virtual functions for filter algorithm */

	virtual Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z);
	virtual Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z);
	/* Observation z(k) and with (Un)correlated observation noise model
		Requires x(k|k), X(k|k) or internal equivilent
		Returns: Reciprocal condition number of primary matrix used in observe computation (1. if none)
		Default implementation simple computes innovation for observe_innovation
	*/

	virtual Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s) = 0;
	virtual Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s) = 0;
	/* Observation innovation s(k) and with (Un)correlated observation noise model
		Requires x(k|k), X(k|k) or internal equivilent
		Returns: Reciprocal condition number of primary matrix used in observe computation (1. if none)
	*/
};



/*
 * Sample Filter: Bayes filter using
 *
 * Probability distributions are represted by a finite sampling
 *
 * State (x_size) size and its sampling (s_size) are assumed to remain constant.
 * The state sampling public so they can be directly manipulated.
 *  init: Should be used if they to be altered
 *
 * The filter is operated by performing a
 * 	predict, observe
 * cycle derived from the bayes_filter. observe Likelihoods are merged into a single combined weight.
 *   update: MUST be used to complete a explict resampling of the particles using merged weights
 *
 * Derived filters supply definititions for the abstract functions and determine the algorithm used
 * to implement the filter.
 */

class Sample_filter : public Likelihood_filter, public Functional_filter
{
public:
	FM::ColMatrix S;		// state sampleing (x_size,s_size)

	Sample_filter (size_t x_size, size_t s_size);
	/* Initialise filter and set constant sizes for
		x_size of the state vector
		s_size sample size
		Exceptions:
			bayes_filter_exception is s_size < 1
	*/
	~Sample_filter() = 0;	// Provide unambigues distructor incase S is not distructable

	/* Virtual functions for filter algorithm */

	virtual void init () = 0;
	/* Initialise from current sampleing
	*/

	virtual void init_sample (const FM::ColMatrix& initS)
	/* Initialise from a sampling
	 */
	{
		S.assign (initS);
		init();
	}

	virtual Float update_resample () = 0;
	/* Resampling update
	 *	Returns lcond, Smallest normalised likelihood weight, represents conditioning of resampling solution
	 *          lcond == 1. if no resampling performed
	 *			This should by multipled by the number of samples to get the Likelihood function conditioning
	 */

	virtual void update ()
	// Default update, simple resample
	{
		(void)update_resample ();
	}

	virtual void predict (Functional_predict_model& f);
	/* Predict state posterior with functional no noise model
	*/

	virtual void predict (Sampled_predict_model& f) = 0;
	/* Predict state posterior with sampled noise model
	*/

	virtual void observe (Likelihood_observe_model& h, const FM::Vec& z) = 0;
	/* Observation state posterior using likelihood model h at z
	*/

	virtual void observe_likelihood (const FM::Vec& lw) = 0;
	/* Observation fusion directly from likelihood weights
	 * lw may be smaller then the state sampling. Weights for additional particles are assumed to be 1
	*/

	size_t unique_samples () const;
	/*
	 * Count number of unique (unequal value) samples in S
	 * Default implementation required std::sort on Samples
	 * Post: S is ordered
	 */
};


}//namespace
#endif
