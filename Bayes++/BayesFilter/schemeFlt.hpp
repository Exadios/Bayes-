#ifndef _BAYES_FILTER_SCHEME
#define _BAYES_FILTER_SCHEME

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Generic Filter
 *  Filter schemes vary in their constructor parameterisation
 *  Filter_scheme derives a generic filter with consistent constructor interface
 *
 *  Provides specialisations for all Bayesian_filter schemes
 */


/* Filter namespace */
namespace Bayesian_filter
{

template <class Scheme>
class Filter_scheme : public Scheme
/*
 * A Generic Filter Scheme
 *  Class template to provide a consistent constructor interface
 *  Specialisations are provided for the Linrz_filter schemes
 */
{
public:
	Filter_scheme(size_t x_size);
	// Common fixed state size

	Filter_scheme(size_t x_size, size_t z_initialsize);
	Filter_scheme(size_t x_size, size_t q_maxsize, size_t z_initialsize);
	// Allow predefined sizes. Useful for some Linrz_filter specialisations
};

template <class Scheme>
Filter_scheme<Scheme>::Filter_scheme(size_t x_size) :
	Scheme (x_size)
{}

// Defaults allowing specialisations for Linrz_filters and Extened_filters
//   They only require x_size, and can make use of z_initalsize
template <class Scheme>
Filter_scheme<Scheme>::Filter_scheme(size_t x_size, size_t z_initalsize) :
	Scheme (x_size, z_initialsize)
{}

template <class Scheme>
Filter_scheme<Scheme>::Filter_scheme(size_t x_size, size_t q_maxsize, size_t z_initialsize) :
	Scheme (x_size, z_initialsize)
{}

// UD_filter specialisation, only one constructor
template <>
class Filter_scheme<UD_scheme> : public UD_scheme
{
public:
	Filter_scheme(size_t x_size, size_t q_maxsize, size_t z_initialsize) :
		UD_scheme (x_size, q_maxsize, z_initialsize)
	{}
};

}//namespace
#endif
