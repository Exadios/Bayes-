#ifndef _BAYES_FILTER_INDIRECT
#define _BAYES_FILTER_INDIRECT

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */
 
/*
 * Indirect Filter Adaptor
 *  Hides the details of the indirect operation of the filter
 *  The error filter uses the same linear models as the direct filter,
 *  observation error computation (subtraction) is linear!
 */

/* Filter namespace */
namespace Bayesian_filter
{


template <typename Error_base>
class Indirect_state_filter : public Expected_state {
/*
 * Indirect state filter
 *  Estimates state using an associated observation error filter
 */
public:
	Indirect_state_filter (Error_base& error_filter)
		: Expected_state(error_filter.x.size()), direct(error_filter)
	{	// Construct and zero initial error filter
		direct.x.clear();
	}

	template <typename P_model>
	void predict (P_model& f)
	{
		x = f.f(x);
		direct.predict(f);				// May be optimised for linear f as direct.x = 0
	};

	void observe (Linear_uncorrelated_observe_model& h, const FM::Vec& z)
	{
				// Observe error (explict temporary)
		FM::Vec z_error (h.h(x) - z);
		observe_error (h, z_error);
	}

	void observe (Linear_correlated_observe_model& h, const FM::Vec& z)
	{
				// Observe error (explict temporary)
		FM::Vec z_error (h.h(x) - z);
		observe_error (h, z_error);
	}

	template <typename O_model>
	void observe_error (O_model& h, const FM::Vec& z_error)
	/*
	 * Observe with precomputed indirect observation error
	 * The observation model here can be non-linear but must then
	 * be modified to return h(x_error) - h(x=0)
	 */
	{
		direct.observe (h, z_error);
		direct.update ();
				// Update State estimate with error
		x -= direct.x;
				// Reset the error
		direct.x.clear();
		direct.init();
	}

	void update ()
	/* Update filters state
	     Updates x(k|k)
	*/
	{}

private:
	Error_base& direct;
};


template <typename Error_base>
class Indirect_kalman_filter : public Kalman_state {
/*
 * Indirect kalman filter
 *  Estimates state using an associated observation error filter
 */
public:
	Indirect_kalman_filter (Error_base& error_filter)
		: Kalman_state(error_filter.x.size()), direct(error_filter)
	{	
	}

	void init ()
	/* Initialise from state and state covariance
	*/
	{
		direct.x.clear();				// Zero initial error
		direct.X = X;
		direct.init();
	}

	template <typename P_model>
	void predict (P_model& f)
	{
		x = f.f(x);
		direct.predict(f);				// May be optimised for linear f as x = 0
	};

	void observe (Linear_uncorrelated_observe_model& h, const FM::Vec& z)
	{
				// Observe error (explict temporary)
		FM::Vec z_error (h.h(x) - z);
		observe_error (h, z_error);
	}

	void observe (Linear_correlated_observe_model& h, const FM::Vec& z)
	{
				// Observe error (explict temporary)
		FM::Vec z_error (h.h(x) - z);
		observe_error (h, z_error);
	}

	template <typename O_model>
	void observe_error (O_model& h, const FM::Vec& z_error)
	/*
	 * Observe with precomputed indirect observation error
	 * The observation model here can be non-linear but must then
	 * be modified to return h(x_error) - h(x=0)
	 */
	{
		direct.observe (h, z_error);
		direct.update();
				// Update State estimate with error
		x -= direct.x;
				// Reset the error
		direct.x.clear();
		direct.init();
	}

	void update ()
	/* Update filters state
	     Updates x(k|k)
	*/
	{
		direct.update();
		X = direct.X;
	}

private:
	Error_base& direct;
};


}//namespace
#endif
