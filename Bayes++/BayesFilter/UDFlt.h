#ifndef _BAYES_FILTER_UD
#define _BAYES_FILTER_UD
/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
 * $Header$
 * $NoKeywords: $
 */

/*
 * UD Factorisation of Covariance Filter.
 *	A Covariance (Kalman) filter as an Abstract class
 * 
 * Bierman's UD factorisatised update algorithm using Agee-Turner UdU' factorisation rank 1 update
 * Thornton's MWG-S factorisation prediction  algorithm
 * References
 * "Factorisation Methods for Discrete Sequential Estimation" Gerald J. Bierman ISBN 0-12-097350-2
 * "Real time Kalman Filter Application", Mohinder S. Grewal, Angus P. Andrews ISBN 0-13-211335-X
 *
 * Derived filters must supply predict and observe model functions.
 * State and control input sizes should remain constant.
 * A initial observation size may also be specified for efficiency.
 * 
 * The filter is operated by performing a
 * 	predict, observe
 * cycle derived from the Linrz_filter
 *
 */
#include "bayesFlt.h"

/* Filter namespace */
namespace Bayesian_filter
{

class UD_sequential_observe_model : public Linrz_uncorrelated_observe_model
{
public:
	UD_sequential_observe_model (FM::Subscript x_size, FM::Subscript z_size) :
		Linrz_uncorrelated_observe_model(x_size, z_size), Hx_o(x_size)
	{}
	virtual FM::Vec& ho (const FM::Vec& x, const FM::Subscript o) = 0;
	/* Supplied model (h) for observation using state x, z allows normalisation and model variation
	   Fast model of a single element (o) in observation model
	   Precondition: Hx_o is conformantly dimensioned
	   Postcondition:
		z(k|k-1) = h(x(k|k-1)
		Hx_o(x(k-1|k-1) = Jacobian of h with respect to state x (row o)
	*/
	FM::Vec Hx_o;
};


class UD_filter : public Linrz_filter
{
private:
	FM::Subscript q_max;	// Maxiumum size allocated for noise model, constructed before UD
public:
	FM::Matrix UD;	// UDU factorisation of X with D on diagonal
						// Lower triangle used as workspace
	FM::Vec s;		// Innovation
	FM::Vec Sd;		// Innovation Covariance 

	UD_filter (FM::Subscript x_size, FM::Subscript q_maxsize, FM::Subscript z_initialsize = 0);
	UD_filter& operator= (const UD_filter&);
	// Optimise copy assignment to only copy filter state

	using Linrz_filter::init;
	void init ();
	void update ();
	Float predict (Linrz_predict_model& f);

	Float observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z);
	Float observe (Linrz_correlated_observe_model& h, const FM::Vec& z);
	/* No solution for Correlated noise and Linrz model */

	Float observe (Linear_correlated_observe_model& h, const FM::Vec& z);
	Float observe (UD_sequential_observe_model& h, const FM::Vec& z);
	/* Special observe using with sequential for fast uncorrelated Linrz operation */

protected:
	Float predictGq (const FM::Matrix& Fx, const FM::Matrix& G, const FM::Vec& q);
	FM::Vec d, dv, v;	// predictGQ temporaries
	Float observeUD (FM::Vec& gain, Float& alpha, const FM::Vec& h, const Float r);
	FM::Vec a, b;		// observeUD temporaries
						// Observation temporaies
	void observe_size (FM::Subscript z_size);
	FM::Subscript last_z_size;
	FM::Vec h1;				// Single Observation model
	FM::Vec w;				// Single Gain
	FM::Vec zp;				// prediction
	FM::Matrix UHx;			// Modified Model for linear decorrelation
	FM::Vec zdecol;			// Decorrelated
	FM::Matrix Gz;
};


}//namespace
#endif
