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
 *	A possibly non-linaer Information filter as an Abstract class
 *
 * Ref:	"Stochastic Models, Estimation, and Control} Peter S Maybeck, Academic Press, ISBN 0-12-480701-1
 *	Section 5.7
 * To work with with Linear and Linrz models
 *  a) a seperate state and information (via covariance) prediction is used.
 *  b) a EIF modified innovation update algorithm is used.
 * Filter maintains y,Y where:
 *		Y = inv(X)
 *		y = Y*x
 *  
 * These require that X, and Y are invertable so state can be computed.
 * Discontinous observe models require that prediction is normailised with
 * respect to the observation.
 *
 * The filter is operated by performing a
 * 	predict, observe
 * cycle defined by the base class
 */
#include "bayesFlt.hpp"

/* Filter namespace */
namespace Bayesian_filter
{

class Information_filter : public Extended_filter
{
public:
	FM::Vec y;				// Information state
	FM::SymMatrix Y;		// Information

	Information_filter (size_t x_size, size_t z_initialsize = 0);
	Information_filter& operator= (const Information_filter&);
	// Optimise copy assignment to only copy filter state

	void init ();
	void init_information (const FM::Vec& y, const FM::SymMatrix& Y);
	void update ();
	Float predict (Linrz_predict_model& f);
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

/*
 * Information Filter with Joseph form prediction
 *	A Linear Predict model Information filter as an Abstract class
 *  Prediction 
 *
 * Ref:	"Stochastic Models, Estimation, and Control} Peter S Maybeck, Academic Press, ISBN 0-12-480701-1
 *	Section 5.7
 * The Joesph form for prediction is used for know inverse prediction model
 * The does NOT require X,Y to be invertable but the inverse prediction model and noise must exist.
 * Discontinous observe models require that prediction is normailised with
 * respect to the observation.
 *
 * Derived filters must supply the plant model function.
 * State and control input sizes should remain constant.
 * The filter is operated by performing a
 * 	predict, observe
 * cycle derived from the linear_filter
 */

class Information_joseph_filter : public Information_filter
{
public:
	Information_joseph_filter (size_t x_size, size_t z_initialsize = 0);

	void predict (Linear_invertable_predict_model& f);

protected:			   		// Permenantly allocated temps
	class Predict_temp
	{	// Space for intermediate computations required by predict
		friend class Information_joseph_filter;
		FM::SymMatrix inv_Q, A, inv_AQ;
		FM::Matrix Chi, IChi;
		FM::Matrix Ywork;
	public:
		Predict_temp (size_t x_size);
	} t;
};


}//namespace
#endif
