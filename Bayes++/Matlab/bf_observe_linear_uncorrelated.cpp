/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_observe_linear_uncorrelated:  Observe z with noise Zv through linear H
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
		mexErrMsgTxt("usage: [ <rcond>= ] observe_linear_uncorrelated(<handle>, <Hx>, <z>, <Z_vector>)");

	// Check handle and get filter pointer
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	// Need filter state size
	Bayesian_filter::State_filter* f = dynamic_cast<Bayesian_filter::State_filter*>(filter);
	Run_Filter::check_def(f);

	// Must catch exceptions for Matlab kernel
	try
	{
		// Arguments
		FM::Matrix Hx = Matrix(prhs[1]);
		FM::Vec    z  = Vector(prhs[2]);
		FM::Vec    Zv = Vector(prhs[3]);
		// Check matrix conformance
		if (f->x.size() != Hx.size2())
			mexErrMsgTxt("Mismatch in x and H matrix");
		if (z.size() != Hx.size1())
			mexErrMsgTxt("Mismatch in z and H matrix");
		if (Zv.size() != z.size())
			mexErrMsgTxt("Mismatch in z and Zv size");

		// Observe model
		Bayesian_filter::General_LiUnAd_observe_model model(Hx.size2(), z.size());
		model.Hx = Hx;
		model.Zv = Zv;

		bool bFilterOp = false;			// Flag operation complete
		double rcond;
		// Must be a Sample, Linrz filter
		if (Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter)) {
			bFilterOp = true;
			f->observe (model, z);
			rcond = f->update_resample ();
		}
		else
		if (Bayesian_filter::Linrz_filter* f = dynamic_cast<Bayesian_filter::Linrz_filter*>(filter)) {
			bFilterOp = true;
			rcond = f->observe (model, z);
		}
				
		if (bFilterOp) {
			if (nlhs == 1)		// Optional return value
				plhs[0] = Array(rcond);
		}
		else {
			Run_Filter::check_def(NULL);
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
