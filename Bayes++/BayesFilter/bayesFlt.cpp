/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */
 
/*
 * Bayesian_filter Implemention of:
 *  BayesFlt constructor/destructors
 *  exceptions
 *  numeric constants
 */
#include "bayesFlt.hpp"
#include <boost/limits.hpp>


/* Filter namespace */
namespace Bayesian_filter
{


/* Minimum allowable reciprocal condition number for PD Matrix factorisations
 * Initialised default gives 5 decimal digits of headroom */
const Bayes_base::Float Numerical_rcond::limit_PD_init = std::numeric_limits<Bayes_base::Float>::epsilon() * Bayes_base::Float(1e5);


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

Linrz_predict_model::Linrz_predict_model (size_t x_size, size_t q_size) :
/*
 * Set the size of things we know about
 */
		Addative_predict_model(x_size, q_size),
		Fx(x_size,x_size)
{
}

Linear_predict_model::Linear_predict_model (size_t x_size, size_t q_size) :
/*
 * Set the size of things we know about
 */
		Linrz_predict_model(x_size, q_size),
		xp(x_size)
{
}

Linear_invertable_predict_model::Linear_invertable_predict_model (size_t x_size, size_t q_size) :
/*
 * Set the size of things we know about
 */
		Linear_predict_model(x_size, q_size),
		inv(x_size, q_size)
{
}

Linear_invertable_predict_model::inverse_model::inverse_model (size_t x_size, size_t q_size) :
		Fx(x_size,x_size)
{
}


State_filter::State_filter (size_t x_size) :
	x(x_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
	if (x_size < 1)
		filter_error ("Zero state filter constructed");
}


Kalman_filter::Kalman_filter (size_t x_size) :
/*
 * Initialise filter and set the size of things we know about
 */
		State_filter(x_size), X(x_size,x_size)
{
}


Linrz_filter::Linrz_filter (size_t x_size) :
/*
 * Initialise filter and set the size of things we know about
 */
		Kalman_filter(x_size)
{
}


Extended_filter::Extended_filter (size_t x_size) :
/*
 * Initialise filter and set the size of things we know about
 */
		Linrz_filter(x_size)
{
}

Information_form_filter::Information_form_filter (size_t x_size) :
/*
 * Initialise filter and set the size of things we know about
 */
		y(x_size), Y(x_size,x_size)
{
}


Sample_filter::Sample_filter (size_t x_size, size_t s_size) :
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
