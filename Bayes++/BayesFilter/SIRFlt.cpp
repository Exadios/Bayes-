/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Sampling Importance Resampleing Filter.
 *
 * Bootstap filter (Sequential Importance Resampleing).
 */
#include "bayesFlt.hpp"
#include "matSup.hpp"
#include "SIRFlt.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <boost/limits.hpp>

namespace {

template <class scalar>
inline scalar sqr(scalar x)
// Square 
{
	return x*x;
}
};//namespace


/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;



Standard_resampler::Float
 Standard_resampler::resample (Resamples_t& presamples, size_t& uresamples, FM::DenseVec& w, SIR_random& r) const
/*
 * Standard resampler from [1]
 * Algorithm:
 *	Complexity O(n*log(n)) complexity required to sort the random draws made.
 *	This allows comparing of the two ordered lists w and ur.
 * Output:
 *  Presamples number of times this particle should be resampled
 *  uresamples number of unqiue particles (number of non zeros in Presamples)
 *  w becomes a normalised cumulative sum
 * Sideeffects:
 *  A draw is made from this instances 'random' for each particle
 */
{
	assert (presamples.size() == w.size());
						// Normalised cumulative sum of likelihood weights, find smallest weight
	Float wcum = 0;
	Float wmin = std::numeric_limits<Float>::max();
	DenseVec::iterator wi, wi_end = w.end();
	for (wi = w.begin(); wi != wi_end; ++wi) {
		if (*wi < wmin) {
			wmin = *wi;
		}
		wcum = *wi = wcum + *wi;
	}
	if (wmin < 0)		// Bad weights
		error (Numeric_exception("negative weight"));
	if (wcum <= 0)		// Bad cumulative weights (previous check should actually prevent -ve
		error (Numeric_exception("zero cumulative weight sum"));
						// Any numerical failure should cascade into cummulative sum
	if (wcum != wcum)	// Inequality due to NaN
		error (Numeric_exception("NaN cumulative weight sum"));

						// Sorted uniform random distribution [0..1) for each resample
	DenseVec ur(w.size());
	r.uniform_01(ur);
	std::sort (ur.begin(), ur.end());
	assert(ur[0] >= 0 && ur[ur.size()-1] < 1);		// Very bad if random is incorrect
						// Scale ur to cummulative sum
	ur *= wcum;

						// Resamples based on cumulative weights from sorted resample random values
	Resamples_t::iterator pri = presamples.begin();
	wi = w.begin();
	DenseVec::iterator ui = ur.begin(), ui_end = ur.end();
	size_t unique = 0;

	while (wi != wi_end)
	{
		size_t Pres = 0;			// assume P not resampled until find out otherwise
		if (ui != ui_end && *ui < *wi)
		{
			++unique;
			do						// count resamples
			{
				++Pres;
				++ui;
			} while (ui != ui_end && *ui < *wi);
		}
		++wi;
		*pri++ = Pres;
	}
	assert(pri==presamples.end());
	if (ui != ui_end)				// Resample failed due no non numeric weights
		error (Numeric_exception("weights are not numeric and cannot be resampled"));
	
	uresamples = unique;
	return wmin / wcum;
}


Systematic_resampler::Float
 Systematic_resampler::resample (Resamples_t& presamples, size_t& uresamples, FM::DenseVec& w, SIR_random& r) const
/*
 * Systematic resample algorithm from [2]
 * Algorithm:
 *	Particles are chosesn from those whose cumulative weight intersect with an equidistant grid
 *	A uniform random draw is chosen to position the grid with the cumulative weights
 *	Complexity O(n)
 * Output:
 *  Presamples number of times this particle should be resampled
 *  uresamples number of unqiue particles (number of non zeros in Presamples)
 *  w becomes a normalised cumulative sum
 * Sideeffects:
 *  A single draw is made from this instances 'random'
 */
{
	size_t nParticles = presamples.size();
	assert (nParticles == w.size());
						// Normalised cumulative sum of likelihood weights, find smallest weight
	Float wcum = 0;
	Float wmin = std::numeric_limits<Float>::max();
	DenseVec::iterator wi, wi_end = w.end();
	for (wi = w.begin(); wi != wi_end; ++wi) {
		if (*wi < wmin) {
			wmin = *wi;
		}
		wcum = *wi = wcum + *wi;
	}
	if (wmin < 0)		// Bad weights
		error (Numeric_exception("negative weight"));
	if (wcum <= 0)		// Bad cumulative weights (previous check should actually prevent -ve
		error (Numeric_exception("total likelihood zero"));
						// Any numerical failure should cascade into cummulative sum
	if (wcum != wcum)
		error (Numeric_exception("total likelihood numerical error"));

						// Stratified step
	Float wstep = wcum / Float(nParticles);
						
	DenseVec ur(1);			// Single uniform for initialisation
	r.uniform_01(ur);
	assert(ur[0] >= 0 && ur[0] < 1);		// Very bad if random is incorrect

						// Resamples based on cumulative weights
	Importance_resampler::Resamples_t::iterator pri = presamples.begin();
	wi = w.begin();
	size_t unique = 0;
	Float s = ur[0] * wstep;	// Random initialisation

	while (wi != wi_end)
	{
		size_t Pres = 0;			// assume P not resampled until find out otherwise
		if (s < *wi)
		{
			++unique;
			do						// count resamples
			{
				++Pres;
				s += wstep;
			} while (s < *wi);
		}
		++wi;
		*pri++ = Pres;
	}
	assert(pri==presamples.end());

	uresamples = unique;
	return wmin / wcum;
}


/*
 * SIR filter implementation
 */
const SIR_scheme::Float SIR_scheme::rougheningKinit = 1;
		// use 1 std.dev. per sample as default roughening

SIR_scheme::SIR_scheme (size_t x_size, size_t s_size, SIR_random& random_helper) :
		Sample_filter (x_size, s_size),
		Sample_state_filter (x_size, s_size),
		random (random_helper),
		resamples (s_size), wir (s_size)
/*
 * Initialise filter and set the size of things we know about
 */
{
	SIR_scheme::x_size = x_size;
	rougheningK = rougheningKinit;
}

SIR_scheme& SIR_scheme::operator= (const SIR_scheme& a)
/* Optimise copy assignment to only copy explict filter state (particles)
 * random helper is not part of filter state
 * Precond: matrix size conformance
 */
{
	*static_cast<Sample_filter*>(this) = a;
	stochastic_samples = S.size2();
	std::fill (wir.begin(), wir.end(), Float(1));		// Initial uniform weights
	wir_update = false;
	return *this;
}


void
 SIR_scheme::init ()
/*
 * Initialise sampling
 *		Pre: S
 *		Post: stochastic_samples := samples in S
 */
{
	stochastic_samples = S.size2();
	std::fill (wir.begin(), wir.end(), Float(1));		// Initial uniform weights
	wir_update = false;
}


SIR_scheme::Float
 SIR_scheme::update_resample (const Importance_resampler& resampler)
/*
 * Resample particles using weights and roughen
 * Pre : S represent the predicted distribution
 * Post: S represent the fused distribution, n_resampled from weighted_resample
 * Exceptions:
 *  Bayes_filter_exception from resampler
 *		unchanged: S, stochastic_samples
 * Return
 *		lcond, Smallest normalised weight, represents conditioning of resampling solution
 *		lcond == 1 if no resampling preformed
 *		This should by multipled by the number of samples to get the Likelihood function conditioning
 */
{
	Float lcond = 1;
	if (wir_update)		// Resampleing only required if weights have been updated
	{
		// Resample based on likelihood weights
		size_t R_unique;
		lcond = resampler.resample (resamples, R_unique, wir, random);

							// No resampling exceptions: update S
		copy_resamples (S, resamples);
		stochastic_samples = R_unique;

		roughen ();			// Roughen samples

		std::fill (wir.begin(), wir.end(), Float(1));		// Resampling results in uniform weights
		wir_update = false;
	}
	return lcond;
}


void
 SIR_scheme::predict (Sampled_predict_model& f)
/*
 * Predict state posterior with sampled noise model
 *		Pre : S represent the prior distribution
 *		Post: S represent the predicted distribution, stochastic_samples := samples in S
 */
{
						// Predict particles S using supplied predict model
	const size_t nSamples = S.size2();
	for (size_t i = 0; i != nSamples; ++i) {
		FM::ColMatrix::Column Si(S,i);
		Si .assign (f.fw(Si));
	}
	stochastic_samples = S.size2();
}


void
 SIR_scheme::observe (Likelihood_observe_model& h, const FM::Vec& z)
/*
 * Observation fusion using Likelihood at z
 * Pre : wir previous particle likelihood weights
 * Post: wir fused (multiplicative) particle likehood weights
 */
{
	h.Lz (z);			// Observe likelihood at z

						// Weight Particles. Fused with previous weight
	const size_t nSamples = S.size2();
	for (size_t i = 0; i != nSamples; ++i) {
		wir[i] *= h.L( FM::column(S,i) );
	}
	wir_update = true;
}

void
 SIR_scheme::observe_likelihood (const FM::Vec& lw)
/*
 * Observation fusion directly from likelihood weights
 * lw may be smaller then the number of particles. Weights for additional particles are assumed to be 1
 * Pre : wir previous particle likelihood weights
 * Post: wir fused (multiplicative) particle likehood weights
 */
{
					// Weight Particles. Fused with previous weight
	Vec::const_iterator lw_end = lw.end();
	for (Vec::const_iterator lw_i = lw.begin(); lw_i != lw_end; ++lw_i) {
		wir[lw_i.index()] *= *lw_i;
	}
	wir_update = true;
}


void SIR_scheme::copy_resamples (FM::ColMatrix& P, const Importance_resampler::Resamples_t& presamples)
/*
 * Update P by selectively copying presamples
 * Uses a in-place copying alogorithm
 * Algorithm: In-place copying
 *  First copy the live samples (those resampled) to end of P
 *  Replicate live sample in-place
 */
{
							// reverse_copy_if live
	size_t si = P.size2(), livei = si;
	Importance_resampler::Resamples_t::const_reverse_iterator pri, pri_end = presamples.rend();
	for (pri = presamples.rbegin(); pri != pri_end; ++pri) {
		--si;
		if (*pri > 0) {
			--livei;
			FM::column(P,livei) .assign (FM::column(P,si));
		}
	}
	assert(si == 0);
							// Replicate live samples
	si = 0;
	Importance_resampler::Resamples_t::const_iterator pi, pi_end = presamples.end();
	for (pi = presamples.begin(); pi != pi_end; ++pi) {
		size_t res = *pi;
		if (res > 0) {
			do  {
				FM::column(P,si) .assign (FM::column(P,livei));
				++si; --res;
			} while (res > 0);
			++livei;
		}
	}
	assert(si == P.size2()); assert(livei == P.size2());
}


void SIR_scheme::roughen_minmax (FM::ColMatrix& P, Float K) const
/*
 * Roughening
 *  Uses algorithm from Ref[1] using max-min in each state of P
 *  K is scaleing factor for roughening noise
 *  unique_samples is unchanged as roughening is used to postprocess observe resamples
 * Numerical colapse of P
 *  P with very small or zero range result in minimal roughening
 * Exceptions:
 *  none
 *		unchanged: P
 */
{
	using namespace std;
						// Scale Sigma by constant and state dimensions
	Float SigmaScale = K * pow (Float(P.size2()), -1/Float(x_size));

						// Find min and max states in all P, precond P not empty
	Vec xmin(x_size); noalias(xmin) = column(P,0);
	Vec xmax(x_size); noalias(xmax) = xmin;
	ColMatrix::iterator2 pi = P.begin2();
	while (pi != P.end2())		// Loop includes 0 to simplify code
	{
		Vec::iterator ni = xmin.begin();
		Vec::iterator xi = xmax.begin();

#ifdef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		for (ColMatrix::iterator1 xpi = begin(pi, ublas::iterator2_tag()); xpi != pi.end(); ++xpi)
#else
		for (ColMatrix::iterator1 xpi = pi.begin(); xpi != pi.end(); ++xpi)
#endif
		{
			if (*xpi < *ni) *ni = *xpi;
			if (*xpi > *xi) *xi = *xpi;
			++ni; ++xi;

		}
		++pi;
	}
   						// Roughening st.dev max-min
	Vec rootq(x_size);
	rootq = xmax;
	noalias(rootq) -= xmin;
	rootq *= SigmaScale;
   						// Apply roughening prediction based on scaled variance
	DenseVec n(x_size);
	for (pi = P.begin2(); pi != P.end2(); ++pi) {
		random.normal (n);			// independant zero mean normal
									// multiply elements by std dev
		for (FM::DenseVec::iterator ni = n.begin(); ni != n.end(); ++ni) {
			*ni *= rootq[ni.index()];
		}
		FM::ColMatrix::Column Pi(P,pi.index2());
		noalias(n) += Pi;			// add to P
		Pi = n;
	}
}


/*
 * SIR implementation of a Kalman filter
 */
SIR_kalman_scheme::SIR_kalman_scheme (size_t x_size, size_t s_size, SIR_random& random_helper) :
	SIR_scheme (x_size, s_size, random_helper),
	Sample_state_filter (x_size, s_size), Kalman_state_filter(x_size), 
	rough_random (random_helper),
	rough (x_size,x_size, rough_random)
{
	FM::identity (rough.Fx);
}

void SIR_kalman_scheme::init_from_kalman ()
/*
 * Initialise sampling from kalman statistics
 *	Pre: x,X
 *	Post: x,X,S
 */
{
						// Samples at mean
	const size_t nSamples = S.size2();
	for (size_t i = 0; i != nSamples; ++i) {
		FM::ColMatrix::Column Si(S,i);
		Si.assign (x);
	}
						// Decorrelate init state noise
	Matrix UD(x_size,x_size);
	Float rcond = UdUfactor (UD, X);
	rclimit.check_PSD(rcond, "Roughening X not PSD");

						// Use predict roughening to apply initial noise
	UdUseperate (rough.G, rough.q, UD);
	rough.init_GqG();
	predict (rough);

	SIR_scheme::init();
}


void SIR_kalman_scheme::mean ()
/*
 * Update state mean
 *		Pre : S
 *		Post: x,S
 */
{
						// Mean of distribution: mean of particles
	x.clear();
	const size_t nSamples = S.size2();
	for (size_t i = 0; i != nSamples; ++i) {
		FM::ColMatrix::Column Si(S,i);
		x.plus_assign (Si);
	}
	x /= Float(S.size2());
}


Bayes_base::Float
 SIR_kalman_scheme::update_resample (const Importance_resampler& resampler)
/*
 * Modified SIR_scheme update implementation
 *  update mean and covariance of sampled distribution with update_statistics
 */
{
	Float lcond = SIR_scheme::update_resample(resampler);	// Resample particles

	update_statistics();			// Estimate sample mean and covariance

	// No assert_isPSD (X) as it may fail due to normal numerical circumstances
	return lcond;
}

void SIR_kalman_scheme::update_statistics ()
/*
 * Update kalman statistics without resampling
 * Update mean and covariance of sampled distribution => mean and covariance of particles
 *		Pre : S (s_size >=1 enforced by state_filter base)
 *		Post: x,X,S	(X may be non PSD)
 *
 * Sample Covariance := (Sum_i [transpose(S[i])*S[i]) - s_size*transpose(mean)*mean) / (s_size-1)
 *  The definition is the unbiased estimate of covariance given samples with unknown (estimated) mean
 *  Implies X is indeterminate when size == 1
 * Semi-definate covariance X due to colapse of S:
 *  If the number of unique samples in S is small relative to the state size then
 *  X will become semi-definate and numerically may appear be negative
 *  As this situation is quite common and due to the normal circumstances and
 *  not an algorithm failure no assertion is made.
 *  Use with care and check the results of any algorithms relying on X
 */
{
	mean ();
	X.clear();				// Covariance

	const size_t nSamples = S.size2();
	for (size_t i = 0; i != nSamples; ++i) {
		FM::ColMatrix::Column Si(S,i);
		X.plus_assign (FM::outer_prod(Si, Si));
	}
	X.minus_assign (FM::outer_prod(x, x) * Float(nSamples));
	X /= Float(nSamples-1);
}


void SIR_kalman_scheme::roughen_correlated (FM::ColMatrix& P, Float K)
/*
 * Roughening
 *  Uses a roughening noise based on covariance of P
 *  This is a more sophisticated algorithm then Ref[1] as it takes
 *  into account the correlation of P
 *  K is scaleing factor for roughening noise
 * Numerical colapse of P
 *  Numerically when covariance of P semi definate (or close), X's UdU factorisation
 *  may be negative.
 * Exceptions:
 *   Bayes_filter_exception due colapse of P
 *		unchanged: P
 */
{
	using namespace std;
						// Scale variance by constant and state dimensions
	Float VarScale = sqr(K) * pow (Float(P.size2()), Float(-2.)/Float(x_size));

	update_statistics();	// Estimate sample mean and covariance

						// Decorrelate states
	Matrix UD(x_size,x_size);
	Float rcond = UdUfactor (UD, X);
	rclimit.check_PSD(rcond, "Roughening X not PSD");

	UdUseperate (rough.G, rough.q, UD);
						// Apply roughening prediction based on scaled variance
	rough.q *= VarScale;
	rough.init_GqG();

						// Predict particles P using rough model
	const size_t nSamples = P.size2();
	for (size_t i = 0; i != nSamples; ++i) {
		FM::ColMatrix::Column Pi(P,i);
		Pi.assign (rough.f(Pi));
	}
}


}//namespace
