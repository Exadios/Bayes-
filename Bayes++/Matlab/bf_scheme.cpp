/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens
 * See Bayes++.htm for copyright license details
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_scheme: Construct a Bayes++ filter scheme and return its handle
 */

#include "matlabBfilter.h"

using namespace::Matlab_convert;


class Bf_maker
/*
 * Matlab support for constructing Bayes++ filters with arguments from a MEX function
 */
{
public:
	typedef Run_Filter* Filter;
	Filter make (const char* filter_scheme, const mxArray* make_args[], int make_nargs);
private:
	typedef Filter (Bf_maker::* pMake_mem)();
	typedef std::map<std::string, pMake_mem> Make_map;
	Make_map filterMakers;
	const mxArray** args;	// Current arguments to make
	int nargs;
	Boost_SIR_random_helper randomHelper;
	// Filter makers
	Filter make_SIR();
	Filter make_SIR_Kalman();
	Filter make_Unscented();
	Filter make_Covariance();
	Filter make_Information();
	Filter make_Information_joseph();
	Filter make_SRIF();
	Filter make_UD();
};


Bf_maker::Bf_maker()
/*
 * Construct name-> make mapping
 */
{
	filterMakers["SIR"] = &Bf_maker::make_SIR;
	filterMakers["SIR_Kalman"] = &Bf_maker::make_SIR_Kalman;
	filterMakers["Unscented"] = &Bf_maker::make_Unscented;
	filterMakers["Covariance"] = &Bf_maker::make_Covariance;
	filterMakers["Information"] = &Bf_maker::make_Information;
	filterMakers["InformationJoseph"] = &Bf_maker::make_Information_joseph;
	filterMakers["SRIF"] = &Bf_maker::make_SRIF;
	filterMakers["UD"] = &Bf_maker::make_UD;
}


Bf_maker::Filter
 Bf_maker::make (const char* filter_scheme, const mxArray* make_args[], int make_nargs)
/*
 * Use Matlab args to Make instance of a filter scheme
 *  args starts at scheme specific arguments
 */
{
								// Find the entry in filterMakers map
	Make_map::iterator i = filterMakers.find (filter_scheme);

	Filter f = NULL;
	if (i != filterMakers.end())
	{							// Rememeber args for scheme maker
		args = make_args; nargs = make_nargs;
								// Call the make member function found to create filter scheme
		pMake_mem pmake = i->second;
		f = (this->*pmake)();
	}
	return f;
};


Bf_maker::Filter
 Bf_maker::make_SIR()
{
	if (nargs < 2 || nargs > 4)
		mexErrMsgTxt("Filter expects: <init_x>, <init_X> [,<sample_size> [,<rougheningK]]");

				// Initialise Sample S with state and covariance
	FM::Vec init_x = Vector(args[0]);
	FM::SymMatrix init_X = SymMatrix(args[1]);
				// Defaultable parameters
	int sample_size = 1000;
	double rougheningK = -1;
	if (nargs > 2)
	{
		sample_size = cint(args[2]);
		if (!(sample_size > 0))
			mexErrMsgTxt ( "sample_size must > 0" );
	}
	if (nargs > 3)
	{
		rougheningK = cdouble(args[3]);
		if (!(rougheningK >= 0))
			mexErrMsgTxt ( "rougheningK must >= 0" );
	}

	Run_SIR* f = new Run_SIR(init_x.size(), sample_size, randomHelper);

	// Use a temporary filter to create samples
	Bayesian_filter::SIR_kalman_filter tempS(f->S.size1(), f->S.size2(), randomHelper);
	tempS.init_kalman (init_x,init_X);
	f->init_sample(tempS.S);
				// Overide filter default
	if (rougheningK != -1.)
		f->rougheningK = rougheningK;

	return f;
}

Bf_maker::Filter
 Bf_maker::make_SIR_Kalman()
{
	if (nargs < 2 || nargs > 4)
		mexErrMsgTxt("Filter expects: <init_x>, <init_X> [,<sample_size> [,<rougheningK]]");

				// Initialise Sample S with state and covariance
	FM::Vec init_x = Vector(args[0]);
	FM::SymMatrix init_X = SymMatrix(args[1]);
				// Defaultable parameters
	int sample_size = 1000;
	double rougheningK = -1;
	if (nargs > 2)
	{
		sample_size = cint(args[2]);
		if (!(sample_size > 0))
			mexErrMsgTxt ( "sample_size must > 0" );
	}
	if (nargs > 3)
	{
		rougheningK = cdouble(args[3]);
		if (!(rougheningK >= 0))
			mexErrMsgTxt ( "rougheningK must >= 0" );
	}

	Run_SIR_kalman* f = new Run_SIR_kalman(init_x.size(), sample_size, randomHelper);

	f->init_kalman (init_x,init_X);
				// Overide filter default
	if (rougheningK != -1.)
		f->rougheningK = rougheningK;

	return f;
}

Bf_maker::Filter
 Bf_maker::make_Unscented()
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

Bf_maker::Filter
 Bf_maker::make_Covariance()
{
	if (nargs != 2)
		mexErrMsgTxt("Filter expects: <init_x>, <init_X>");

	// Initialise with state and covariance
	FM::Vec init_x = Vector(args[0]);
	FM::SymMatrix init_X = SymMatrix(args[1]);

	Run_Cov* f = new Run_Cov(init_x.size());
	f->init_kalman (init_x,init_X);
	return f;
}

Bf_maker::Filter
 Bf_maker::make_Information()
{
	if (nargs != 2)
		mexErrMsgTxt("Filter expects: <init_x>, <init_X>");

	// Initialise with state and covariance
	FM::Vec init_x = Vector(args[0]);
	FM::SymMatrix init_X = SymMatrix(args[1]);

	Run_Inf* f = new Run_Inf(init_x.size());
	f->init_kalman (init_x,init_X);
	return f;
}

Bf_maker::Filter
 Bf_maker::make_Information_joseph()
{
	if (nargs != 2)
		mexErrMsgTxt("Filter expects: <init_x>, <init_X>");

	// Initialise with state and covariance
	FM::Vec init_x = Vector(args[0]);
	FM::SymMatrix init_X = SymMatrix(args[1]);

	Run_InfJo* f = new Run_InfJo(init_x.size());
	f->init_kalman (init_x,init_X);
	return f;
}

Bf_maker::Filter
 Bf_maker::make_SRIF()
{
	if (nargs != 2)
		mexErrMsgTxt("Filter expects: <init_x>, <init_X>");

	// Initialise with state and covariance
	FM::Vec init_x = Vector(args[0]);
	FM::SymMatrix init_X = SymMatrix(args[1]);

	Run_SRIF* f = new Run_SRIF(init_x.size());
	f->init_kalman (init_x,init_X);
	return f;
}

Bf_maker::Filter
 Bf_maker::make_UD()
{
	if (nargs != 2)
		mexErrMsgTxt("UD filter expects:  <init_x>, <init_X>");

	// Initialise with state and covariance
	FM::Vec init_x = Vector(args[0]);
	FM::SymMatrix init_X = SymMatrix(args[1]);

	Run_UD* f = new Run_UD(init_x.size());
	f->init_kalman (init_x,init_X);
	return f;
}


// Maintain a set of default filters to make
static Bf_maker DefaultFilters;


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
		mexErrMsgTxt("usage: <handle> = bf_scheme(<filter_scheme_string>, ...)");

	/*
	 * Create a filter from its scheme name and remaining arguments
	 * - Must catch exceptions for Matlab kernel
	 */
	Run_Filter* filter = NULL;
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
