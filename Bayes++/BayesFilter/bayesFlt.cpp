/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
 * $Header$
 * $NoKeywords: $
 */
 
/*
 * Bayesian_filter Implemention of:
 *  BayesFlt constructor/destructors
 *  exceptions
 *  numeric constants
 */
#include "matSup.h" 
#include <boost/limits.hpp>

#include "bayesFlt.h"

/* Filter namespace */
namespace Bayesian_filter
{


/* Minimum allowable reciprocal condition number for PD Matrix factorisations
 * Initialised default gives 5 decimal digits of headroom */
const Bayes_base::Float Numerical_rcond::limit_PD_init = std::numeric_limits<Bayes_base::Float>::epsilon() * 1e5;


Bayes_base::~Bayes_base()
/*
 * Default definition required for a pure virtual distructor
 */
{}

void Bayes_base::filter_error (const char* errorText) const
/*
 * Generic error
 */
{
	throw Bayes_filter_exception (errorText);
}

Linrz_predict_model::Linrz_predict_model (FM::Subscript x_size, FM::Subscript q_size) :
/*
 * Set the size of things we know about
 */
		Addative_predict_model(x_size, q_size),
		Fx(x_size,x_size)
{
}

Linear_predict_model::Linear_predict_model (FM::Subscript x_size, FM::Subscript q_size) :
/*
 * Set the size of things we know about
 */
		Linrz_predict_model(x_size, q_size),
		xp(x_size)
{
}

Linear_invertable_predict_model::Linear_invertable_predict_model (FM::Subscript x_size, FM::Subscript q_size) :
/*
 * Set the size of things we know about
 */
		Linear_predict_model(x_size, q_size),
		inv(x_size, q_size)
{
}

Linear_invertable_predict_model::inverse_model::inverse_model (FM::Subscript x_size, FM::Subscript q_size) :
		Fx(x_size,x_size),
		q(q_size), G(x_size,q_size)
{
}


State_filter::State_filter (FM::Subscript x_size) :
	x(x_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
	if (x_size < 1)
		filter_error ("Zero state filter constructed");
}


Kalman_filter::Kalman_filter (FM::Subscript x_size) :
/*
 * Initialise filter and set the size of things we know about
 */
		State_filter(x_size), X(x_size,x_size)
{
}


Linrz_filter::Linrz_filter (FM::Subscript x_size) :
/*
 * Initialise filter and set the size of things we know about
 */
		Kalman_filter(x_size)
{
}


Extended_filter::Extended_filter (FM::Subscript x_size) :
/*
 * Initialise filter and set the size of things we know about
 */
		Linrz_filter(x_size)
{
}


Sample_filter::Sample_filter (FM::Subscript x_size, FM::Subscript s_size) :
		Likelihood_filter(),
		S(x_size,s_size)

/*
 * Initialise filter and set the size of things we know about
 * Postcond: s_size >= 1
 */
{
	if (s_size < 1)
		filter_error ("Zero sample filter constructed");
}


}//namespace
