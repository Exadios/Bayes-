/*
 * MATLAB MEX interface for Bayes++
 *  bf_mean:  Get filter state mean vector
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
	if (nrhs != 1)
		usageOk = false;
	if (nlhs > 1)
		usageOk = false;

	if (!usageOk)
		mexErrMsgTxt("usage: [<state_mean_vec>=] mean(<handle>)");

	/*
	 * Check handle and get filter pointer
	 */
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	/*
	 * Must catch exceptions for Matlab kernel
	 */
	try
	{
		// Must be a State filter
		Bayesian_filter::State_filter* f = dynamic_cast<Bayesian_filter::State_filter*>(filter);
		Run_Filter::check_def(f);

		// Update and get state
		f->update();
		plhs[0] = Array(f->x);
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
