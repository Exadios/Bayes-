/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 */

/*
 * MATLAB MEX interface for Bayes++
 *  sample:
 *  Return sample array reresenting filters state probability distribution.
 *  If the filter is represented by samples these direcly used
 *  otherwise a set of samples are generated to represent it
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
	if (nlhs != 1)
		usageOk = false;

	if (!usageOk)
		mexErrMsgTxt("usage: <sample_matrix> = sample(<handle>)");

	/*
	 * Check handle and get filter pointer
	 */
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	/*
	 * Must catch exceptions for Matlab kernel
	 */
	try
	{
		// Sample filter is easy
		{
			Bayesian_filter::Sample_filter* sf = dynamic_cast<Bayesian_filter::Sample_filter*>(filter);
			if (sf) {
				plhs[0] = ArrayTranspose(sf->S);
				return;
			}
		}

		// Kalman filter requres more work. Create Samples from a mean and covariance of Kalman filter
		{
			unsigned GenerateSamples = 1000;
			Bayesian_filter::Kalman_filter* kf = dynamic_cast<Bayesian_filter::Kalman_filter*>(filter);
			if (kf) {
				kf->update();
				Boost_SIR_random_helper randomHelper;
				Bayesian_filter::SIR_kalman_filter tempS(kf->x.size(), GenerateSamples, randomHelper);
				tempS.init_kalman (kf->x, kf->X);
				plhs[0] = ArrayTranspose(tempS.S);
				return;
			}
		}

		// Fall through error if type unknown
		Run_Filter::check_def(NULL);
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
