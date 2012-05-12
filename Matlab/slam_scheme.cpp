/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_scheme: Constuct a Bayes++ filter scheme and return its handle
 */

#include <mex.h>
#include "matlabBfilter.h"
#include "SLAM/SLAM.hpp"
#include "SLAM/fastSLAM.hpp"
#include "SLAM/kalmanSLAM.hpp"

using namespace::Matlab_convert;

/*
 * Runtime filter : Encapsulation for a filter schemse
 */
struct Run_Fast_SLAM_Kstatistics : public Run_Filter, public SLAM_filter::Fast_SLAM_Kstatistics
{
	Run_Fast_SLAM_Kstatistics( Bayesian_filter::SIR_kalman_filter& L_filter ) :
		Fast_SLAM_Kstatistics(L_filter)
	{}
};

struct Run_Kalman_SLAM : public Run_Filter, public SLAM_filter::Kalman_SLAM
{
	Run_Kalman_SLAM( Bayesian_filter::Linrz_filter& full_filter, unsigned full_nL) :
		Kalman_SLAM(full_filter, full_nL)
	{}
};


class Slam_maker
/*
 * Matlab support for constructing Bayes++ SLAM filters with arguments from a MEX function
 */
{
public:
	typedef Run_Filter* Filter;
	Filter make (const char* filter_scheme, const mxArray* make_args[], int make_nargs);
private:
	typedef Filter (Slam_maker::* pMake_mem)();	// MemFun retuning a Filter
	typedef std::map<std::string, pMake_mem> Make_map;
	Make_map filterMakers;
	const mxArray** args;	// Current arguments to make
	int nargs;
	Boost_SIR_random_helper randomHelper;
	// SLAM Filter makers
	Filter make_Fast();
	Filter make_Kalman();
};

Slam_maker::Slam_maker()
/*
 * Construct name-> make mapping
 */
{
	filterMakers["Fast"] = &Slam_maker::make_Fast;
	filterMakers["Kalman"] = &Slam_maker::make_Kalman;
}

Slam_maker::Filter
 Slam_maker::make (const char* filter_scheme, const mxArray* make_args[], int make_nargs)
/*
 * Use Matlab args to Make instance of a filter scheme
 *  args starts at scheme specific arguments
 */
{
								// Find the entry in filterMakers map
	Make_map::iterator i = filterMakers.find (filter_scheme);

	Filter f = 0;
	if (i != filterMakers.end())
	{							// Rememeber args for scheme maker
		args = make_args; nargs = make_nargs;
								// Call the make member function found to create filter scheme
		pMake_mem pmake = i->second;
		f = (this->*pmake)();
	}
	return f;
};


Slam_maker::Filter
 Slam_maker::make_Fast()
{
	if (nargs < 0 || nargs > 3)
		mexErrMsgTxt("Filter expects: <state_size> [,<sample_size> [,<rougheningK]]");

	int state_size = cint(args[0]);
	if (!(state_size > 0))
		mexErrMsgTxt ( "state_size must > 0" );
				// Initialise Sample S with state and covariance
	FM::Vec init_x = Vector(args[0]);
	FM::SymMatrix init_X = SymMatrix(args[1]);
				// Defaultable parameters
	int sample_size = 1000;
	double rougheningK = -1;
	if (nargs > 1)
	{
		sample_size = cint(args[1]);
		if (!(sample_size > 0))
			mexErrMsgTxt ( "sample_size must > 0" );
	}
	if (nargs > 2)
	{
		rougheningK = cdouble(args[1]);
		if (!(rougheningK >= 0))
			mexErrMsgTxt ( "rougheningK must >= 0" );
	}

				// Make support filter scheme and Fast SLAM filter
	Run_SIR_kalman* lf = new Run_SIR_kalman(state_size, sample_size, randomHelper);
				// Overide filter default
	if (rougheningK != -1.)
		lf->rougheningK = rougheningK;
	Run_Fast_SLAM_Kstatistics* f = new Run_Fast_SLAM_Kstatistics(*lf);

	return f;
}

Slam_maker::Filter
 Slam_maker::make_Kalman()
{

	if (nargs != 2)
		mexErrMsgTxt("Filter expects: <init_x>, <init_X>");

	// Initialise with state and covariance
	FM::Vec init_x = Vector(args[0]);
	FM::SymMatrix init_X = SymMatrix(args[1]);

	Run_Unscented* f = new Run_Unscented(init_x.size());
	f->init_kalman (init_x,init_X);
	return f;
}

// Maintain a set of default filters to make
static Slam_maker DefaultFilters;


/*
 * MEX Create a Bayesfilter
 */

void mexFunction(
		 int          nlhs,
		 mxArray      *plhs[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
	/* Check for proper usage */
	bool usageOk = true;
	if (nrhs < 1)	// rhs always includes filter scheme
		usageOk = false;
	if (nlhs != 1)	// lhs only filter handle
		usageOk = false;

	char* filter_scheme;
	if (usageOk) {
		filter_scheme = cstring(prhs[0]);
		if (!filter_scheme)
			usageOk = false;
	}
	if (!usageOk)
		mexErrMsgTxt("usage: <handle> = slam_scheme(<filter_scheme_string>, ...)");

	/*
	 * Create a filter from its scheme name and remaining arguments
	 * - Must catch exceptions for Matlab kernel
	 */
	Run_Filter* filter = 0;
	try
	{
		filter = DefaultFilters.make (filter_scheme, &prhs[1], nrhs-1);
		if (!filter)
			mexErrMsgTxt( "<filter_scheme> name unknown" );
	}
	catch (std::exception& se)
	{
		mexErrMsgTxt( se.what() );
	}
	catch (...)
	{
		mexErrMsgTxt( "Filter operation caused an unknown exception" );
	}

	// Use return value as filter handler
	if (filter)
	{
		plhs[0] = filter->make_handle();
	}
	else
		mexErrMsgTxt("Internal error: never reached");
}
