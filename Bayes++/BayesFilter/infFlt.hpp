#ifndef _BAYES_FILTER_INFORMATION
#define _BAYES_FILTER_INFORMATION

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Information Filter.
 *  A possibly non-linear Information filter as an Abstract class
 *
 * References
 * [1] "Stochastic Models, Estimation, and Control} Peter S Maybeck, Academic Press, ISBN 0-12-480701-1
 *      Section 5.7
 * [2] "Kalman Filtering, Theory and Practice", Mohinder S. Grewal, Angus P. Andrews ISBN 0-13-211335-X
 * To work with with Linear and Linrz models
 *  a) a seperate state and information (via covariance) prediction is used.
 *  b) a EIF modified innovation update algorithm is used.
 *
 * Two alternative algorithms are used for predict functions:
 *  For linrz models an extended predict form is used so information state 'y' is predicted via the
 *  non-linear function. This requires that X, and Y are invertable so 'x' can be computed.
 * For linear invertable models prediction can be done directly without computing x
 * Discontinous observe models require that prediction is normailised with
 * respect to the observation.
 *
 * The filter is operated by performing a
 *  predict, observe
 * cycle defined by the base class
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

class Information_filter : public Information_form_filter, public Extended_filter
{
public:
	Information_filter (size_t x_size, size_t z_initialsize = 0);
	Information_filter& operator= (const Information_filter&);
	// Optimise copy assignment to only copy filter state

	struct Linear_predict_byproducts
	{
		 Linear_predict_byproducts (size_t x_size, size_t q_size);
		 FM::Matrix tempA;
		 FM::SymMatrix A;
		 FM::Matrix tempG;
		 FM::SymMatrix B, invB;
		 FM::Matrix tempY;
	};

	void init ();
	void update ();
	void init_yY ();
	void update_yY ();
	// Covariance and information form state interface

	Float predict (Linear_invertable_predict_model& f, Linear_predict_byproducts& b);
	// Linear Prediction in information form as in Ref[2]
	Float predict (Linear_invertable_predict_model& f)
	{
		Linear_predict_byproducts b(f.Fx.size1(),f.q.size());
		return predict (f, b);	
	}

	Float predict (Linrz_predict_model& f);
	// Extended Prediction via state

	Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s);
	Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s);

protected:
	bool update_required;	// Postcondition of update is not met

protected:			   		// Permenantly allocated temps
	FM::Vec i;
	FM::SymMatrix I;
	FM::SymMatrix ZI;
					// allow fast operation if z_size remains constant
	size_t last_z_size;
	void observe_size (size_t z_size);
};


}//namespace
#endif
