#ifndef _BAYES_FILTER_SIR
#define _BAYES_FILTER_SIR

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $NoKeywords: $
 */

/*							
 * Sampling Importance Resampleing Filter.
 *  Bootstrap filter as an Abstract class
 *
 * References
 *  [1] "Novel approach to nonlinear-non-Guassian Bayesian state estimation"
 *   NJ Gordon, DJ Salmond, AFM Smith IEE Proceeding-F Vol.140 No.2 April 1993
 *  [2] Building Robust Simulation-based Filter for Evolving Data Sets"
 *   J Carpenter, P Clifford, P Fearnhead Technical Report Unversity of Oxford
 *
 *  A variety of reampling algorithms can be used for the SIR filter.
 *  There are implementations for two algorihtms:
 *	 standard_resample: Standard resample algorithm from [1]
 *	 systematic_resample: A Simple stratified resampler from [2]
 *  A virtual 'weighted_resample' provides an standard interface to these and defaults
 *	to the standard_resample.
 *
 * NOTES:
 *  SIR algorithm is sensative to random generator
 *  In particular random uniform must be [0..1) NOT [0..1]
 *  Quantisation in the random number generator must not approach the sample size.
 *  This will result in quantisation of the resampling.
 *  For example if random identically equal to 0 becomes highly probable due to quantisation
 *  this will result in the first sample being selectively draw whatever its likelihood.
 *
 *  Numerics
 *   Resampling requires comparisons of normalised weights. These may
 *   become insignificant if Likelihoods have a large range. Resampling becomes ill conditioned
 *   for these samples.
 */
#include "bayesFlt.h"
#include "models.h"

/* Filter namespace */
namespace Bayesian_filter
{

struct SIR_random
/* Random number generators interface
 *  Helper to allow polymorthic use of random number generators 
 */
{
	virtual void normal(FM::Vec& v) = 0;
	virtual void uniform_01(FM::Vec& v) = 0;
};


class Importance_resampler : public Bayes_base
/*
 * Importance resampler
 *  Represents a function that computes the posterior resampling baes on importance weights
 *  Polymorphic function object use to parameterise the resampling operation
 */
{
public:
	typedef std::vector<size_t> Resamples_t;	// resampling counts

	virtual Float resample (Resamples_t& presamples, size_t& uresamples, FM::Vec& w, SIR_random& r) const = 0;
	/*
	 * The resampling function
	 *  Weights w are proportional to the posterior Likelihood of a state
	 * Sideeffect
	 *  w becomes a normalised cumulative sum
	 *  Random draws can be made from 'r'
	 *
	 * Exceptions:
	 *  bayes_filter_exception for numerical problems with weights including
	 *   any w < 0, all w == 0
	 * Return
	 *	lcond, smallest normalised weight, represents conditioning of resampling solution
	 *	presamples: posterior resample, the number of times each sample should appear posterior baes on it weight
	 *  uresample: the number of unique resamples in the posterior == number of non zero elements in presamples
	 * Preconditions
	 *  wresample,w must have same size
	 */
};

class Standard_resampler : public Importance_resampler
// Standard resample algorithm from [1]
{
	Float resample (Resamples_t& presamples, size_t& uresamples, FM::Vec& w, SIR_random& r) const;
};

class Systematic_resampler : public Importance_resampler
// Systematic resample algorithm from [2]
{
	Float resample (Resamples_t& presamples, size_t& uresamples, FM::Vec& w, SIR_random& r) const;
};


class SIR_filter : public Sample_filter
/*
 * Sampling Importance Resampleing Filter.
 *  Implement a general form of SIR filter.
 *  Importance resampling is delayed until an update is required. The sampler used
 *  is a parameter of update to allow a wide variety of usage.
 *  A stochastic sample is defined as a sample with a unqiue stochastic history other then roughening
 */
{
	friend class SIR_kalman_filter;
public:
	size_t stochastic_samples;	// Number of stochastic samples in S

	SIR_filter (size_t x_size, size_t s_size, SIR_random& random_helper);
	SIR_filter& operator= (const SIR_filter&);
	// Optimise copy assignment to only copy filter state

	/* Specialisations for filter algorithm */
	void init ();

	Float update_resample ()
	// Default resampling update
	{
		return update_resample (Standard_resampler());
	}

	virtual Float update_resample (const Importance_resampler& resampler);
	/* Update: resample particles using weights and then roughen
	 *  Return: lcond
	 */

	void predict (Sampled_predict_model& f);
	// Predict particles with sampled noise model

	void observe (Likelihood_observe_model& h, const FM::Vec& z);
	// Weight particles using likelihood model h and z

	void observe_likelihood (const FM::Vec& lw);
	// Observation fusion directly from likelihood weights

	Float rougheningK;			// Current roughening value (0 implies no roughening)
	virtual void roughen()
	{	// Generalised roughening:  Default to roughen_minmax
		if (rougheningK != 0.)
			roughen_minmax (S, rougheningK);
	}

	static void copy_resamples (FM::ColMatrix& P, const Importance_resampler::Resamples_t& presamples);
	// Update P by selectively copying based on presamples 

	SIR_random& random;			// Reference random number generator helper

protected:
	void roughen_minmax (FM::ColMatrix& P, Float K) const;	// Roughening using minmax of P distribution
	Importance_resampler::Resamples_t resamples;		// resampling counts
	FM::Vec wir;				// resamping weights
	bool wir_update;			// weights have been updated requring a resampling on update
private:
	static const Float rougheningKinit;
	size_t x_size;
};



class SIR_kalman_filter : public SIR_filter, public Kalman_filter
/*
 * SIR implementation of a Kalman filter
 *  Updates Kalman statistics of SIR_filter
 *  These statistics are use to provide a specialised correlated roughening procedure
 */
{
public:
	using Kalman_filter::x;
	SIR_kalman_filter (size_t x_size, size_t s_size, SIR_random& random_helper);

	/* Specialisations for filter algorithm */
	void Kalman_filter::init();
	// Initialisation from kalman statistics

	// Init from S is already defined by SIR_filter

	virtual void update ()
	// Implement Default Kalman_filter update
	{
		Float ignore_lcond = SIR_filter::update_resample();
	}

	Float update_resample (const Importance_resampler& resampler);
	// Modified SIR_filter update implementation: update mean and covariance of sampled distribution

	void update_statistics ();
	// Update kalman statistics without resampling

	void roughen()
	{	// Specialised correlated roughening
		if (rougheningK != 0.)
			roughen_correlated (S, rougheningK);
	}


protected:
	void roughen_correlated (FM::ColMatrix& P, Float K);	// Roughening using covariance of P distribution
	// Roughening model and random normal
	class Roughen_random : public General_LiAd_predict_model::Random
	{
	public:
		Roughen_random(SIR_random& random_helper) :
			random(random_helper)
		{}
		virtual void normal(FM::Vec& v)
		{
			random.normal(v);
		}
	
	private:
		SIR_random& random;
	} rough_random;
	General_LiAd_predict_model rough;
private:
	static Float scaled_vector_square(const FM::Vec& v, const FM::SymMatrix& S);
	void mean();
};


}//namespace
#endif
