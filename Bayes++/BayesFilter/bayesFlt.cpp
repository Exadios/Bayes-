/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens
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

void Bayes_base::error (const Numeric_exception& e )
/*
 * Filter error
 */
{
	throw e;
}

void Bayes_base::error (const Logic_exception& e )
/*
 * Filter error
 */
{
	throw e;
}

Gaussian_predict_model::Gaussian_predict_model (size_t x_size, size_t q_size) :
		q(q_size), G(x_size, q_size)
/*
 * Set the size of things we know about
 */
{}

Addative_predict_model::Addative_predict_model (size_t x_size, size_t q_size) :
		q(q_size), G(x_size, q_size)
/*
 * Set the size of things we know about
 */
{}

Linrz_predict_model::Linrz_predict_model (size_t x_size, size_t q_size) :
/*
 * Set the size of things we know about
 */
		Addative_predict_model(x_size, q_size),
		Fx(x_size,x_size)
{}

Linear_predict_model::Linear_predict_model (size_t x_size, size_t q_size) :
/*
 * Set the size of things we know about
 */
		Linrz_predict_model(x_size, q_size),
		xp(x_size)
{}

Linear_invertable_predict_model::Linear_invertable_predict_model (size_t x_size, size_t q_size) :
/*
 * Set the size of things we know about
 */
		Linear_predict_model(x_size, q_size),
		inv(x_size)
{}

Linear_invertable_predict_model::inverse_model::inverse_model (size_t x_size) :
		Fx(x_size,x_size)
{}


State_filter::State_filter (size_t x_size) :
	x(x_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
	if (x_size < 1)
		error (Logic_exception("Zero state filter constructed"));
}


Kalman_state_filter::Kalman_state_filter (size_t x_size) :
/*
 * Initialise state size
 */
		State_filter(x_size), X(x_size,x_size)
{}

void Kalman_state_filter::init_kalman (const FM::Vec& x, const FM::SymMatrix& X)
/*
 * Initialise from a state and state covariance
 *  Parameters that reference the instance's x and X members are valid
 */
{
	Kalman_state_filter::x = x;
	Kalman_state_filter::X = X;
	init ();
}


Information_state_filter::Information_state_filter (size_t x_size) :
/*
 * Initialise state size
 */
		y(x_size), Y(x_size,x_size)
{}

void Information_state_filter::init_information (const FM::Vec& y, const FM::SymMatrix& Y)
/*
 * Initialise from a information state and information
 *  Parameters that reference the instance's y and Y members are valid
 */
{
	Information_state_filter::y = y;
	Information_state_filter::Y = Y;
	init_yY ();
}


Sample_state_filter::Sample_state_filter (size_t x_size, size_t s_size) :
		S(x_size,s_size)

/*
 * Initialise state size
 * Postcond: s_size >= 1
 */
{
	if (s_size < 1)
		error (Logic_exception("Zero sample filter constructed"));
}

void Sample_state_filter::init_sample (const FM::ColMatrix& initS)
/*
 * Initialise from a sampling
 */
{
	S = initS;
	init_S();
}


Sample_filter::Sample_filter (size_t x_size, size_t s_size) :
		Sample_state_filter(x_size,s_size)

/*
 * Initialise filter and set the size of things we know about
 * Postcond: s_size >= 1
 */
{
}


}//namespace
