/*
 * MATLAB MEX interface for Bayes++
 *  bf_predict_functional:  Predict filter with zero noise through function f(x)
 */

#include <mex.h>
#include "bfilter.h"

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
	if (nrhs != 2)
		usageOk = false;
	if (nlhs != 0)
		usageOk = false;

	if (!usageOk)
		mexErrMsgTxt("usage: predict_functional(<handle>, <function>)");

	// Check handle and get filter pointer
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	// Must be a Functional filter
	Bayesian_filter::Functional_filter* f = dynamic_cast<Bayesian_filter::Functional_filter*>(filter);
	Run_Filter::check_def(f);

	class FF : public Bayesian_filter::Functional_predict_model
	{	// Functional form to turn a Function_model into a Functional_predict_model
		Bayesian_filter::Function_model& ff;
	public:
		FF (Bayesian_filter::Function_model& f_init) :
			ff(f_init)
		{}
		virtual const FM::Vec& fx(const FM::Vec& x) const
		{	return ff.fx(x);
		}
	};

	// Must catch exceptions for Matlab kernel
	try
	{
		// Arguments
		Matlab_function_model matlab_f(prhs[1]);

		// Predict filter using Matlab function
		f->predict (FF(matlab_f));
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
