/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * "compatibilty.h" - Work around compiler and Library and deficiencies
 * Introduces namespace BayesFiltercompatibility
 */

#ifndef _BAYES_FILTER_COMPATIBILITY
#define _BAYES_FILTER_COMPATIBILITY

/* Filter Namespace */
namespace Bayesian_filter
{

#if (defined(_MSC_VER) && _MSC_VER <= 1200)
// place sqrt in compat namespace

// If CMATH is included
#if defined(_CMATH_)
 namespace compatibility
 {
	using ::sqrt; 
 }
#endif

#else
 namespace compatibility
 {
	using std::sqrt; 
 }
#endif

}//namespace
#endif
