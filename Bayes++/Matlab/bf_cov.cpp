/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens
 * See Bayes++.htm for copyright license details
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_cov:  Get filter state mean vector and covariance matrix
 */

#include <mex.h>
#include <exception>
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
	if (nlhs != 2)
		usageOk = false;

	if (!usageOk)
		mexErrMsgTxt("usage: [<state_mean_vec>,<state_covariance_matrix>] = cov(<handle>)");

	/*
	 * Check handle and get filter pointer
	 */
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	/*
	 * Must catch exceptions for Matlab kernel
	 */

	try
	{
		// Must be a Kalman filter
		Bayesian_filter::Kalman_filter* f = dynamic_cast<Bayesian_filter::Kalman_filter*>(filter);
		Run_Filter::check_def(f);

		// Update and get state
		f->update();
		plhs[0] = Array(f->x);
		plhs[1] = Array(f->X);
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
