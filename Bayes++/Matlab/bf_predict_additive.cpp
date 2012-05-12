/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_predict_additive:  Predict filter through function f(x) with Gq additive noise
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
	if (nrhs != 4)
		usageOk = false;
	if (nlhs > 1)
		usageOk = false;

	if (!usageOk)
		mexErrMsgTxt("usage: [ <rcond>= ]predict_additive(<handle>, <function>, <G>, <q>)");

	// Check handle and get filter pointer
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	// Need filter state size
	Bayesian_filter::State_filter* f = dynamic_cast<Bayesian_filter::State_filter*>(filter);
	Run_Filter::check_def(f);

	// Must catch exceptions for Matlab kernel
	try
	{
		// Arguments
		Matlab_function_model matlab_f(prhs[1]);
		FM::Vec    q = Vector(prhs[2]);
		FM::Matrix G = Matrix(prhs[3]);
		// Check matrix conformance
		if (f->x.size() != G.size1())
			mexErrMsgTxt("Mismatch in x and G size");
		if (q.size() != G.size2())
			mexErrMsgTxt("Mismatch in q and G size");

		bool bFilterOp = false;	// Flag operation complete
		double rcond;

					// Predict filter using Matlab function
		// TODO Requires General_LzAd_predict_model
		//		f->predict (FF(matlab_f));
		mexErrMsgTxt("Not implemented");
		if (bFilterOp) {
			if (nlhs == 1)		// Optional return value
				plhs[0] = Array(rcond);
		}
		else {
			Run_Filter::check_def(0);
		}
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
