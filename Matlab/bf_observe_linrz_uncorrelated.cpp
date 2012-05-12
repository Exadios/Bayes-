/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_observe_linrz_uncorrelated:  Observe z with noise Zv through non-linear h(x)
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
	if (nrhs != 5)
		usageOk = false;
	if (nlhs > 1)
		usageOk = false;

	if (!usageOk)
		mexErrMsgTxt("usage: [ <rcond>= ] observe_linrz_uncorrelated(<handle>, <h_function>, <Hx>, <z>, <Z_vector>)");

	// Check handle and get filter pointer
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	// Need filter state size
	Bayesian_filter::State_filter* f = dynamic_cast<Bayesian_filter::State_filter*>(filter);
	Run_Filter::check_def(f);

	// Must catch exceptions for Matlab kernel
	try
	{
		// Arguments
		Matlab_function_model matlab_h(prhs[1]);
		FM::Matrix Hx = Matrix(prhs[2]);
		FM::Vec    z  = Vector(prhs[3]);
		FM::Vec    Zv = Vector(prhs[4]);
		// Check matrix conformance
		if (f->x.size() != Hx.size2())
			mexErrMsgTxt("Mismatch in x and H matrix");
		if (z.size() != Hx.size1())
			mexErrMsgTxt("Mismatch in z and H matrix");
		if (Zv.size() != z.size())
			mexErrMsgTxt("Mismatch in z and Zv size");

		// Observe model
		Bayesian_filter::Simple_linrz_uncorrelated_observe_model model(matlab_h, Hx, Zv);

		bool bFilterOp = false;			// Flag operation complete
		double rcond;
		if (Bayesian_filter::Linrz_filter* f = dynamic_cast<Bayesian_filter::Linrz_filter*>(filter)) {
			bFilterOp = true;
			rcond = f->observe (model, z);
		}
				
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
