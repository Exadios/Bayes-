#ifndef _BAYES_FILTER_AVERAGE1
#define _BAYES_FILTER_AVERAGE1

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */
 
/*
 * Predefined filter: average1_filter
 *  A single state averager
 */

/* Filter namespace */
namespace Bayesian_filter
{
	namespace FM = Bayesian_filter_matrix;

template <class filter>
class average1_filter : private filter
{
	typedef typename filter::Float Float;
	class Cpredict : public Linear_predict_model
	// Constant prediction model
	{
	public:
		Cpredict(Float qq) : Linear_predict_model(1, 1)
		{
			Fx(0,0) = 1.;
			q[0] = qq;
			G(0,0) = 1.;
		}
	};

	class Cobserve : public Linear_correlated_observe_model
	// Constant observe model
	{
	public:
		Cobserve(Float ZZ) : Linear_correlated_observe_model(1,1)
		{
			Hx(0,0) = 1.;
			Z(0,0) = ZZ;
		}
	};

public:
	average1_filter (Float iQ, Float iZ, Float z);
	average1_filter (Float iQ, Float iZ);
	Float observe (Float zz);
	operator Float () const
	/* Returns filtered estimate
	 */
	{	if (!bInit)	filter_error("average1 not init");
		return x[0];
	}

private:
	Cpredict f;
	Cobserve h;
	
	bool bInit;
	FM::Vec z;
};



template <typename filter>
average1_filter<filter>::average1_filter (Float iQ, Float iZ, Float zz)
	: filter(1), f(iQ), h(iZ), z(1)
/* Initialise noises and set sizes
 * include first observation zz */
{
	bInit = false;
	observe (zz);
}

template <typename filter>
average1_filter<filter>::average1_filter (Float iQ, Float iZ)
	: filter(1), f(iQ), h(iZ), z(1)
// Initialise noises and set sizes
{
	bInit = false;
}

template <typename filter>
average1_filter<filter>::Float average1_filter<filter>::observe(Float zz)
/* Observe z, first call set initial state to z
 * Returns filtered estimate
 */
{
	z[0] = zz;

	if (!bInit) {
		init_kalman (z, h.Z);
		bInit = true;
	}

	filter::predict(f);
	filter::observe(h, z);
	update ();

	return x[0];
}

}//namespace
#endif
