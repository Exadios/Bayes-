/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 */

/*
 * matlabBfilter.h - Matlab MEX support for Bayes++
 * Defines Bayes++ filter functions for the mex enviroment of Matlab
 * Povides:
 *   Dymanic contructions with RTTI
 *	 Handle managment of dymanically contructed filter filters
 *   Implicitly definied Random number generators
 *	 Function object conversion
 * NOTE: Header file contains definitions!
 *  This makes sense as each MEX function is seperately compiled
 */

// Shut up VC++, it has its own problems!
#if defined(_MSC_VER) && _MSC_VER <= 1200
#pragma warning(disable:4786)	// indentifier more 255 chars
#endif

#include <typeinfo.h>
#include <string>
#include <map>

#include <mex.h>				// Matlab MEX suport
								// Need bounds check so we get exceptions
#define BOOST_UBLAS_BOUNDS_CHECK
#include "BayesFilter/allFilters.hpp" // Include all of Bayesian Filtering library
#include <boost/random.hpp>		// Fast and good random numbers

namespace FM = Bayesian_filter_matrix;
#include "matlabConvert.h"




class Run_Filter
/*
 * Runtime abstract filter : RTTI encapsulation for a filter
 *  By encapsulating all filter classes used at runtime we avoid the need for RTTI
 *  support in filtering library
 */
{
private:
	Run_Filter* signature_data;		// Use 'this' as Unique signature to check memory really is a Run_filter
public:
	Run_Filter()
	{	// Fill out signature
		signature_data = this;
	}
	virtual ~Run_Filter() = 0
	{	// Destory signature
		signature_data = 0;
	}; 

		// Matlab handle interface
	static Run_Filter* from_handle( const mxArray* ma );
	mxArray* make_handle();
	static void check_def(const Bayesian_filter::Bayes_filter_base* derived);
};

Run_Filter* Run_Filter::from_handle( const mxArray* fap )
/*
 * Get and check handle in mxArray of UINT32 realy is a pointer to a Run_filter
 */
{
	if (mxGetClassID(fap) != mxUINT32_CLASS || mxIsComplex(fap) || mxGetM(fap)!=1 || mxGetN(fap)!=1)
		mexErrMsgTxt("Filter <handle> is not of the type created by bfilter");

						// ASSUME we can store filter pointer in the mxUINT32 of fap
	Run_Filter* filter = *reinterpret_cast<Run_Filter**>(mxGetPr(fap));
						// Gross check to see we don't have an invalid pointer
	if (!filter)
		mexErrMsgTxt("Filter <handle> value does not represent a filter object");
						// Check memory has correct signature
	if (filter->signature_data != filter)
		mexErrMsgTxt("Filter <handle> value does not represent a filter object");
	return filter;
}

mxArray* Run_Filter::make_handle()
/*
 * Create a numeric array as handle for Run_filter
 * ASSUME  can store filter pointer in the mxUINT32 element of mxArray
 */
{
	mxArray* afp  = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
	*reinterpret_cast<Run_Filter**>(mxGetPr(afp)) = this;
	return afp;
}

void Run_Filter::check_def(const Bayesian_filter::Bayes_filter_base* derived)
/*
 * Check that derived is defined. If it isn't generate error
 */
{
	std::string error;
	if (!derived) {
		error = std::string("Operation undefine on this filter type");
		mexErrMsgTxt( error.c_str() );
	}
}


/*
 * Runtime filter : Encapsulation for a filter schemse
 */
struct Run_SIR : public Run_Filter, public Bayesian_filter::SIR_filter
{
	Run_SIR(unsigned x, unsigned s, Bayesian_filter::SIR_random& r) :
		SIR_filter(x,s,r) {}
};

struct Run_SIR_kalman : public Run_Filter, public Bayesian_filter::SIR_kalman_filter
{
	Run_SIR_kalman(unsigned x, unsigned s, Bayesian_filter::SIR_random& r) :
		SIR_kalman_filter(x,s,r) {}
};

struct Run_Unscented : public Run_Filter, public Bayesian_filter::Unscented_filter
{
	Run_Unscented(unsigned x) : Unscented_filter(x) {}
};

struct Run_Cov : public Run_Filter, public Bayesian_filter::Covariance_filter
{
	Run_Cov(unsigned x) : Covariance_filter(x) {}
};

struct Run_Inf : public Run_Filter, public Bayesian_filter::Information_filter
{
	Run_Inf(unsigned x) : Information_filter(x) {}
};

struct Run_InfJo : public Run_Filter, public Bayesian_filter::Information_joseph_filter
{
	Run_InfJo(unsigned x) : Information_joseph_filter(x) {}
};

struct Run_SRIF : public Run_Filter, public Bayesian_filter::Information_root_filter
{
	Run_SRIF(unsigned x) : Information_root_filter(x) {}
};

struct Run_UD : public Run_Filter, public Bayesian_filter::UD_filter
{
	Run_UD(unsigned x) : UD_filter(x, x) {}
};



class Boost_random
/*
 * Random number distributions
 */
{
public:
	Boost_random() : gen_normal(rng), gen_uniform(rng)
	{
	}
	double normal(const double mean, const double sigma)
	{
		boost::normal_distribution<boost::mt19937> gen(rng, mean, sigma);
		return gen();
	}
	void normal(FM::Vec& v)
	{
		std::generate (v.begin(), v.end(), gen_normal);
	}
	void uniform_01(FM::Vec& v)
	{
		std::generate (v.begin(), v.end(), gen_uniform);
	}
private:
	boost::mt19937 rng;
	boost::normal_distribution<boost::mt19937> gen_normal;
	boost::uniform_01<boost::mt19937> gen_uniform;
};

// Maintain state of random numbers between function innovations
static Boost_random Random;


class Boost_SIR_random_helper : public Bayesian_filter::SIR_random
/*
 * Random number generator for SIR_filter
 * Uses global Random
 */
{
public:
	Boost_SIR_random_helper()
	{}
	void normal(FM::Vec& v)
	{
		::Random.normal(v);
	}
	void uniform_01(FM::Vec& v)
	{
		::Random.uniform_01(v);
	}
};

class Boost_Predict_random_helper : public Bayesian_filter::General_LiAd_predict_model::Random
/*
 * Random number generator for General_LiAd_predict_model
 * Uses global Random
 */
{
public:
	void normal(FM::Vec& v)
	{
		::Random.normal(v);
	}
};

class Matlab_function_model : public Bayesian_filter::Function_model
/*
 * Function model wrapper for a MuPad function
 */
{
public:
	Matlab_function_model(const mxArray* fname);
	virtual const FM::Vec& fx(const FM::Vec& x) const;
	// Note: Reference return value as a speed optimisation, MUST be copied by caller.
private:
	char* function_name;
	mutable FM::Vec rfx;
	mutable mxArray *plhs[1], *prhs[1];
};


Matlab_function_model::Matlab_function_model(const mxArray* fname) :
	rfx(FM::Empty)
{
	function_name = Matlab_convert::cstring(fname);
    /* must be a string. */
	if (!function_name)
		mexErrMsgTxt("Function name must be a string");
	prhs[0] = 0;
}

const FM::Vec& Matlab_function_model::fx(const FM::Vec& x) const
{
					// Create x argument, assume size is constant
	if (!prhs[0])
		prhs[0] = mxCreateDoubleMatrix(x.size(),1, mxREAL);
	Matlab_convert::Array(prhs[0], x);

					// Call Matlab function by name
	int bad = mexCallMATLAB(1, plhs, 1, prhs, function_name);
	if (bad)
		mexErrMsgTxt( "Bad call to Matlab function");
	else {
		FM::Vec& fx = Matlab_convert::Vector(plhs[0]);
		mxDestroyArray(plhs[0]);

		rfx.resize(fx.size());
		rfx.assign(fx);
	}
	return rfx;
};
