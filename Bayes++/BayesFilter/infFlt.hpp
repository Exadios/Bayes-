#ifndef _BAYES_FILTER_INFORMATION
#define _BAYES_FILTER_INFORMATION

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Information Filter Scheme.
 *  A possibly non-linear Information filter
 *
 * References
 * [1] "Stochastic Models, Estimation, and Control} Peter S Maybeck, Academic Press, ISBN 0-12-480701-1
 *      Section 5.7
 * [2] "Kalman Filtering, Theory and Practice", Mohinder S. Grewal, Angus P. Andrews ISBN 0-13-211335-X
 * To work with with Linear and Linrz models
 *  a) a seperate state and information (via covariance) predict is used.
 *  b) a EIF modified innovation update algorithm is used.
 *
 * Two alternative algorithms are used for predict functions:
 *  For linrz models an extended predict form is used so information state 'y' is predicted via the
 *  non-linear function. This requires that X, and Y are invertable so 'x' can be computed.
 * For linear invertable models predict can be done directly without computing x
 * Discontinous observe models require that predict is normailised with
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

class Information_scheme : public Extended_kalman_filter, virtual public Information_state_filter
{
public:
	Information_scheme (size_t x_size);
	Information_scheme& operator= (const Information_scheme&);
	// Optimise copy assignment to only copy filter state

	struct Predict_linear_byproduct
	{
		 Predict_linear_byproduct (size_t x_size, size_t q_size);
		 FM::SymMatrix A;
		 FM::Matrix tempG;
		 FM::SymMatrix B;
		 FM::Vec y;
	};

	void init ();
	void update ();
	void init_yY ();
	void update_yY ();
	// Covariance and information form state interface

	Float predict (Linrz_predict_model& f);
	// Extended_kalman_filter predict via state
	Float predict (Linear_invertable_predict_model& f);
	// Linear predict, without byproduct
	Float epredict (Linear_invertable_predict_model& f, Predict_linear_byproduct& b);
	// Linear predict with explict byproduct. In information form as in Ref[2]

	Float observe_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s);
	Float observe_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s);
	// Extended_kalman_filter observe
	Float eobserve_innovation (Linrz_uncorrelated_observe_model& h, const FM::Vec& s,
				State_byproduct& i, Covariance_byproduct& I);
	Float eobserve_innovation (Linrz_correlated_observe_model& h, const FM::Vec& s,
				State_byproduct& i, Covariance_byproduct& I);
	// Observe with explict byproduct

protected:
	bool update_required;	// Postcondition of update is not met

protected:			   		// Permenantly allocated temps
	FM::RowMatrix tempX;
};


}//namespace
#endif
