#ifndef _BAYES_FILTER_INDRIECT
#define _BAYES_FILTER_INDIRECT
/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
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
	namespace FM = Bayesian_filter_matrix;


template <typename Error_base>
class Indirect_state_filter : public State_filter {
/*
 * Indirect state filter
 *  Estimates state using an associated observation error filter
 */
public:
	Indirect_state_filter (Error_base& error_filter)
		: State_filter(error_filter.x.size()), direct(error_filter)
	{	// Construct and zero initial error filter
		FM::set(direct.x, 0.);
	}

	template <typename P_model>
	void predict (P_model& f)
	{
		x = f.f(x);
		direct.predict(f);				// May be optimised for linear f as x = 0
	};

	template <typename O_model>
	void observe (O_model& h, const FM::Vec& z)
	{
				// Observe error (explict temporary)
		FM::Vec z_error(z.size());
		z_error = h.h(x);
		z_error -= z;
		direct.observe (h, z_error);
				// Update State estimate with error
		x -= direct.x;
				// Reset the error
		FM::set (direct.x, 0.);
	}

	template <typename O_model>
	void observe_error (O_model& h, const FM::Vec& z_error)
	{
		direct.observe (h, z_error);
				// Update State estimate with error
		x -= direct.x;
				// Reset the error
		FM::set (direct.x, 0.);
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
class Indirect_kalman_filter : public Kalman_filter {
/*
 * Indirect kalman filter
 *  Estimates state using an associated observation error filter
 */
public:
	Indirect_kalman_filter (Error_base& error_filter)
		: Kalman_filter(error_filter.x.size()), direct(error_filter)
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

	template <typename O_model>
	void observe (O_model& h, const FM::Vec& z)
	{
				// Observe error (explict temporary)
		FM::Vec z_error(z.size());
		z_error = h.h(x);
		z_error -= z;
		direct.observe (h, z_error);
		direct.update();
				// Update State estimate with error
		x -= direct.x;
				// Reset the error
		direct.x.clear();
		direct.init ();
	}

	template <typename O_model>
	void observe_error (O_model& h, const FM::Vec& z_error)
	{
		direct.observe (h, z_error);
				// Update State estimate with error
		x -= direct.x;
				// Reset the error
		FM::set (direct.x, 0.);
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
