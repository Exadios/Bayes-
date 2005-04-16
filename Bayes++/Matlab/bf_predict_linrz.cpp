/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_predict_linrz:  Predict filter through non-linear f(x) with Gq addative noise
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
		mexErrMsgTxt("usage: [ <rcond>= ] predict_linrz(<handle>, <f_function>, <Fx>, <G>, <q>)");

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
		FM::Matrix Fx = Matrix(prhs[2]);
		FM::Matrix G = Matrix(prhs[3]);
		FM::Vec    q = Vector(prhs[4]);
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
		Bayesian_filter::Simple_linrz_predict_model model(matlab_f, Fx, G, q);
char warn[1000];
const FM::Vec& fx = matlab_f.fx(f->x);
sprintf(warn, "x%d, fx%d, Fx%dX%d, G%dX%d, q%d", f->x.size(), fx.size(), Fx.size1(),Fx.size2(), G.size1(),G.size2(), q.size());
mexWarnMsgTxt(warn);

		bool bFilterOp = false;			// Flag operation complete
		double rcond;
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
