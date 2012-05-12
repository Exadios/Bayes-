/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_unique_samples:  Get no of unique (different value) samples in Filter
 */

#include <mex.h>
#include "matlabBfilter.h"

using namespace::Matlab_convert;


void mexFunction(
		 int          nlhs,
		 mxArray      *plhs[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
	/* Check for proper usage */
	bool usageOk = true;
	if (nrhs != 1)
		usageOk = false;
	if (nlhs > 1)
		usageOk = false;

	if (!usageOk)
		mexErrMsgTxt("usage: [<number_unique>=] bf_unique_samples(<handle>)");

	/*
	 * Check handle and get filter pointer
	 */
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	/*
	 * Must catch exceptions for Matlab kernel
	 */
	try
	{
		// Must be a Sample filter
		Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter);
		Run_Filter::check_def(f);

		unsigned nsample = f->unique_samples();
		plhs[0] = Array(nsample);
	}
	catch (std::exception& se)
	{
		mexErrMsgTxt( se.what() );
	}
	catch (...)
	{
		mexErrMsgTxt( "Filter operation caused an unknown exception" );
	}
}
