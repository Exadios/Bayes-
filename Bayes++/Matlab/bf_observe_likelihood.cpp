/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 */

/*
 * MATLAB MEX interface for Bayes++
 *  bf_observe_likelihood:  Observe filter with Likelihood L(z,x)
 */

#include <mex.h>
#include "matlabBfilter.h"

using namespace::Matlab_convert;


class Matlab_Likelihood_observe_model : public Bayesian_filter::Likelihood_observe_model
/*
 * Likelihood observe model using a MuPad expression
 */
{
public:
	Matlab_Likelihood_observe_model(unsigned z_size, const mxArray* fname);
	~Matlab_Likelihood_observe_model();
	Float L(const FM::Vec& x) const;
	void Lz(const FM::Vec& z);
private:
	char* function_name;
	mxArray* M_z;
};

Matlab_Likelihood_observe_model::Matlab_Likelihood_observe_model(unsigned z_size, const mxArray* fname) :
	Likelihood_observe_model(z_size)
{
	function_name = Matlab_convert::cstring(fname);
    /* must be a string. */
	if (!function_name)
		mexErrMsgTxt("Function name must be a string");

	M_z = NULL;
}

Matlab_Likelihood_observe_model::~Matlab_Likelihood_observe_model()
{
	// No need to free M_z: mxArray automatically garbage collect
}

Matlab_Likelihood_observe_model::Float
Matlab_Likelihood_observe_model::L(const FM::Vec& x) const
{
					// Create argument
	mxArray *plhs[1], *prhs[2];
	prhs[0] = M_z;
	prhs[1] = Matlab_convert::Array(x);
					// Call Matlab function by name
	int bad = mexCallMATLAB(1, plhs, 2, prhs, function_name);
	if (bad)
		mexErrMsgTxt( "Unable to call Matlab function");

					// Extract returned likelihood
	Float L = cdouble(plhs[0]);
	return L;
};

void Matlab_Likelihood_observe_model::Lz(const FM::Vec& z)
{
	// No need to free M_z: mxArray automatically garbage collect
	M_z = Matlab_convert::Array(z);
};


void mexFunction(
		 int          nlhs,
		 mxArray      *plhs[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
	/* Check for proper usage */
	bool usageOk = true;
	if (nrhs != 3)
		usageOk = false;
	if (nlhs > 1)
		usageOk = false;

	if (!usageOk)
		mexErrMsgTxt("usage: [ <lcond>= ] observe_likelihood(<handle>, <likelihood_function>, <observation_vector>)");

	// Check handle and get filter pointer
	Run_Filter* filter = Run_Filter::from_handle(prhs[0]);

	// Must be a Sample filter
	Bayesian_filter::Sample_filter* f = dynamic_cast<Bayesian_filter::Sample_filter*>(filter);
	Run_Filter::check_def(f);

	// Must catch exceptions for Matlab kernel
	try
	{
					// Arguments
		FM::Vec z = Vector(prhs[2]);
		Matlab_Likelihood_observe_model model(z.size(), prhs[1]);
					// Observe filter using Matlab function
		double lcond;
		f->observe(model, z);
		lcond = f->update_resample();
				
		if (nlhs == 1)		// Optional return value
			plhs[0] = Array(lcond);
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
