#ifndef _BAYES_FILTER_INFORMATION_ROOT
#define _BAYES_FILTER_INFORMATION_ROOT

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Information Root Filter.
 *	A extended 'Square-root' Information filter as an Abstract class
 *
 * Algorithm: Square-root information propogation using QR factorisation
 * Ref:	P. Dyer and S. McReynolds, "Extension of Square-Root Filtering to Include Process Noise",
 *		Journal of Optimization Theory and Applications, Vol.3 No.6 1969
 * Filter maintains r,R where
 * 		inv(R)*inv(R)' = X
 *		r = R*x
 *  R is upper triangular but not strictly a Cholesky factor as diagonal may be negative
 * Observe algorithm has been extended to include linearised models
 * Discontinous observe models require that prediction is normailised with respect to the observation.
 *
 * The filter is operated by performing a
 * 	predict, observe
 * cycle defined by the base class
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

class Information_root_filter : public Extended_filter
{
public:
	FM::Vec r;			// Information Root state
	FM::UTriMatrix R;	// Information Root

	Information_root_filter (size_t x_size, size_t z_initialsize = 0);

	void init ();
	void init_information (const FM::Vec& y, const FM::SymMatrix& Y);
	void update ();
	Float predict (Linrz_predict_model& f);
	/* Generialised from, requires invserion of Fx */
	Float predict (Linear_invertable_predict_model& f);
	/* Specialised linear prediction avoiding inversion of Fx */

	Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s);
	Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s);
};


/*
 * Information Root Filter with exposed information state
 * Augments Information_root_filter with y,Y in the interface
 */

class Information_root_info_filter : public Information_root_filter
{
public:
	FM::Vec y;				// Information state
	FM::SymMatrix Y;		// Information

	Information_root_info_filter (size_t x_size, size_t z_initialsize = 0);

	void update ();
};


}//namespace
#endif
