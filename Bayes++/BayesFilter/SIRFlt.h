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
 * Bootstap filter (Sequential Importance Resampleing).
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

class SIR_filter : public Sample_filter
/*
 * Sampling Importance Resampleing Filter.
 *  A stochastic sample := a sample with a unqiue stochastic history other then roughening
 */
{
	friend class SIR_kalman_filter;
public:
	size_t stochastic_samples;	// Number of stochastic samples in S
	typedef std::vector<size_t> resamples_t;	// resampling counts

	struct Random
	/* Random number generators interface
	 *  Helper to allow polymorthic use of random number generators 
	 */
	{
		virtual void normal(FM::Vec& v) = 0;
		virtual void uniform_01(FM::Vec& v) = 0;
	};

	SIR_filter (size_t x_size, size_t s_size, Random& random_helper);
	SIR_filter& operator= (const SIR_filter&);
	// Optimise copy assignment to only copy filter state

	/* Specialisations for filter algorithm */
	void init (const FM::ColMatrix& S);

	void update (Float& lcond);
	/* Resample particles using weights and roughen
	 *  Return: lcond smallest normalised weight, represents conditioning of resampling solution
	 *          lcond == 1. if no resampling preformed
	 *			This should by multipled by the number of samples to get the Likelihood function conditioning
	 */

	void predict (Sampled_predict_model& f);
	// Predict particles with sampled noise model

	void observe (Likelihood_observe_model& h, const FM::Vec& z);
	// Weight particles using likelihood model h and z

	Float rougheningK;			// Current roughening value (0 implies no roughening)
	virtual void roughen()
	{	// Generalised roughening:  Default to roughen_minmax
		if (rougheningK != 0.)
			roughen_minmax (S, rougheningK);
	}

	virtual Float weighted_resample (resamples_t& Presamples, size_t& uresamples, FM::Vec& w) const
	/*
	* Resample from P(rior) using on importance weights w
	*  Weights w are proportional to the posterior Likelihood of the state represented by a particle
	*  Selects particles from priot to represent the posterior distribution based on their weights w
	* Sideeffect
	*  w becomes a normalised cumulative sum
	*  Draws are made from this instances 'random'
	*
	* Exceptions:
	*  bayes_filter_exception for numerical problems with weights including
	*   any w < 0, all w == 0
	* Return
	*	Smallest normalised weight, represents conditioning of resampling solution
	*	Presamples: the number of times each sample in P is resampled in posterior
	*   uresample: the number of unique resamples in the posterior == number of non zero elements in Presamples
	* Preconditions
	*  Presample,w must  have same size
	*/
	{				// Default to use standard resample algorithm from [1]
		return standard_resample(Presamples, uresamples, w);
	}

	Float standard_resample (resamples_t& Presamples, size_t& uresamples, FM::Vec& w) const;
	// Standard resample algorithm from [1]

	Float systematic_resample (resamples_t& Presamples, size_t& uresamples, FM::Vec& w) const;
	// Systematic resample algorithm from [2]

	static void copy_resamples (FM::ColMatrix& P, const resamples_t& Presamples);
	// Update P by selectively copying Presamples 

protected:			   		// Permenantly allocated temps
	void roughen_minmax (FM::ColMatrix& P, Float K) const;	// Roughening using minmax of P distribution
	resamples_t resamples;		// resampling counts
	FM::Vec wir;				// resamping weights
	bool wir_update;			// weights have been updated requring a resampling on update
private:
	static const Float rougheningKinit;
	Random& random;		// Reference Random number generator
	size_t x_size;
};


class SIR_kalman_filter : public SIR_filter, public Kalman_filter
/*
 * SIR implementation of a Kalman filter
 * Note: 
 *  Provided same interface but is not a Sample_filter, samples only used for layered implementation
 */
{
public:
	using Kalman_filter::x;
	SIR_kalman_filter (size_t x_size, size_t s_size, Random& random_helper);

	/* Specialisations for filter algorithm */
	void init ();
	void update ()
	// Modified update with correlated roughening
	{
		Float lcond_ignore;
		update (lcond_ignore);
	}

	void update (Float& lcond);
	// Modified update with correlated roughening


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
		Roughen_random(SIR_filter::Random& random_helper) :
			random(random_helper)
		{}
		virtual void normal(FM::Vec& v)
		{
			random.normal(v);
		}
	
	private:
		SIR_filter::Random& random;
	} rough_random;
	General_LiAd_predict_model rough;
private:
	static Float scaled_vector_square(const FM::Vec& v, const FM::SymMatrix& S);
	void mean();
	void covariance();
};


}//namespace
#endif
