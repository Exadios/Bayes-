/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens
 * See accompanying Bayes++.htm for terms and conditions of use.
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_observe_new:  Observe new SLAM feature
 */

#include <mex.h>
#include "matlabBfilter.h"
#include "SLAM/SLAM.hpp"


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
	if (nlhs != 0)
		usageOk = false;

	if (!usageOk)
		mexErrMsgTxt("usage: slam_observe_new(<handle>, <feature_no>, <h>, <z>, <Z_vector>)");

	// Check handle and get filter pointer
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	// Must catch exceptions for Matlab kernel
	try
	{
		// Arguments
		unsigned feature = cunsigned(prhs[1]);
		Matlab_function_model matlab_h(prhs[2]);
		FM::Vec    z  = Vector(prhs[3]);
		FM::Vec    Zv = Vector(prhs[4]);
		// Check matrix conformance
		if (1 != Hx.size2())
			mexErrMsgTxt("Mismatch in x and H matrix");
		if (z.size() != Hx.size1())
			mexErrMsgTxt("Mismatch in z and H matrix");
		if (Zv.size() != z.size())
			mexErrMsgTxt("Mismatch in z and Zv size");

		// Observe model
		Bayesian_filter::Linear_observe_model model(Hx.size2(), z.size());
		model.Hx = Hx;
		model.Zv = Zv;

		// Must be a SLAM filter
		SLAM_filter::SLAM* f = dynamic_cast<SLAM_filter::SLAM*>(filter);
		Run_Filter::check_def(f);
		f->observe_new (feature, model, z);
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
