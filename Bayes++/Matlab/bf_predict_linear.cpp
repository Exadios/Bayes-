/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_predict_linear:  Predict filter through Linear F with Gq addative noise
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
		mexErrMsgTxt("usage: [ <rcond>= ]predict_linear(<handle>, <Fx>, <G>, <q>)");

	// Check handle and get filter pointer
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	// Need filter state size
	Bayesian_filter::State_filter* f = dynamic_cast<Bayesian_filter::State_filter*>(filter);
	Run_Filter::check_def(f);

	// Must catch exceptions for Matlab kernel
	try
	{
		// Arguments
		FM::Matrix Fx = Matrix(prhs[1]);
		FM::Matrix G = Matrix(prhs[2]);
		FM::Vec    q = Vector(prhs[3]);
		// Check matrix conformance
		if (f->x.size() != Fx.size1())
			mexErrMsgTxt("Mismatch in x and Fx matrix");
		if (Fx.size1() != Fx.size2())
			mexErrMsgTxt("Mismatch Fx not square");
		if (f->x.size() != G.size1())
			mexErrMsgTxt("Mismatch in x and G size");
		if (q.size() != G.size2())
			mexErrMsgTxt("Mismatch in q and G size");

		// Predict model
		Boost_Predict_random_helper random;
		Bayesian_filter::General_LiAd_predict_model model(Fx.size1(), q.size(), random);
		model.Fx = Fx;
		model.q = q;
		model.G = G;

		bool bFilterOp = false;			// Flag operation complete
		double rcond;
		// Must be a Sample filter or Linrz
		if (Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter)) {
			bFilterOp = true;
			f->predict (model);
			rcond = 0.;
		}
		else
		if (Bayesian_filter::Linrz_filter* f = dynamic_cast<Bayesian_filter::Linrz_filter*>(filter)) {
			bFilterOp = true;
			rcond = f->predict (model);
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
