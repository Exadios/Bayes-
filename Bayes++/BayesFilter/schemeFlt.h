#ifndef _BAYES_FILTER_SCHEME
#define _BAYES_FILTER_SCHEME

/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
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
	Filter_scheme(FM::Subscript x_size);
	// Common fixed state size

	Filter_scheme(FM::Subscript x_size, FM::Subscript z_initialsize);
	Filter_scheme(FM::Subscript x_size, FM::Subscript q_maxsize, FM::Subscript z_initialsize);
	// Allow predefined sizes. Useful for some Linrz_filter specialisations
};

template <class Scheme>
Filter_scheme<Scheme>::Filter_scheme(FM::Subscript x_size) :
	Scheme (x_size)
{}

// Defaults allowing specialisations for Linrz_filters
//   They only require x_size, and can make use of z_initalsize
template <class Scheme>
Filter_scheme<Scheme>::Filter_scheme(FM::Subscript x_size, FM::Subscript z_initalsize) :
	Scheme (x_size, z_initialsize)
{}

template <class Scheme>
Filter_scheme<Scheme>::Filter_scheme(FM::Subscript x_size, FM::Subscript q_maxsize, FM::Subscript z_initialsize) :
	Scheme (x_size, z_initialsize)
{}

// UD_filter specialisation, only one constructor
template <>
class Filter_scheme<UD_filter> : public UD_filter
{
public:
	Filter_scheme(FM::Subscript x_size, FM::Subscript q_maxsize, FM::Subscript z_initialsize) :
		UD_filter (x_size, q_maxsize, z_initialsize)
	{}
};

}//namespace
#endif
