/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
 * $Header$
 * $NoKeywords: $
 *
 */

/*
 * Sampling Importance Resampleing Filter.
 *  Bootstrap filter as an Abstract class
 *
 * Bootstap filter (Sequential Importance Resampleing).
 */
#include "matSup.h"
#include <algorithm>
#include <iterator>
#include <vector>
#include <limits>

#include "SIRFlt.h"

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


// use 1 std.dev. per sample as default roughening
const SIR_filter::Float SIR_filter::rougheningKinit = 1.;

SIR_filter::SIR_filter (Subscript x_size, Subscript s_size, Random& random_helper) :
		Sample_filter(x_size, s_size),
		resamples(s_size), wir(s_size),
		random(random_helper)
/*
 * Initialise filter and set the size of things we know about
 */
{
	SIR_filter::x_size = x_size;
	rougheningK = rougheningKinit;
}

SIR_filter& SIR_filter::operator= (const SIR_filter& a)
/* Optimise copy assignment to only copy explict filter state (particles)
 * random helper is not part of filter state
 * Precond: matrix size conformance
 */
{
	*static_cast<Sample_filter*>(this) = a;
	stochastic_samples = S.size2();
	std::fill (wir.begin(), wir.end(), 1.);		// Initial uniform weights
	wir_update = false;
	return *this;
}


void
 SIR_filter::init (const ColMatrix& SI)
/*
 * Initialise sampling
 *		Post: S, stochastic_samples := samples in S
 */
{
	S = SI;			// Simple initialisation
	stochastic_samples = S.size2();
	std::fill (wir.begin(), wir.end(), 1.);		// Initial uniform weights
	wir_update = false;
}


void
 SIR_filter::update (Float& lcond)
/*
 * Resample particles using weights and roughen
 * Pre : S represent the predicted distribution
 * Post: S represent the fused distribution, n_resampled from weighted_resample
 * Exceptions:
 *  Bayes_filter_exception from resampler
 *		unchanged: S, stochastic_samples
 * Return
 *		lcond: Smallest normalised weight, represents conditioning of resampling solution
 */
{
	lcond = 1.;
	if (wir_update)		// Resampleing only required if weights have been updated
	{
		// Resample based on likelihood weights
		unsigned R_unique; 
		Float lcond = weighted_resample (resamples, R_unique, S, wir);

							// No resampling exceptions: update S
		copy_resamples (S, resamples);
		stochastic_samples = R_unique;

		roughen ();			// Roughen samples

		std::fill (wir.begin(), wir.end(), 1.);		// Initial uniform weights
		wir_update = false;
	}
}


void
 SIR_filter::predict (Sampled_predict_model& f)
/*
 * Predict state posterior with sampled noise model
 *		Pre : S represent the prior distribution
 *		Post: S represent the predicted distribution, stochastic_samples := samples in S
 */
{
						// Predict particles S using supplied predict model
	for (ColMatrix::iterator si = S.begin(); si != S.end(); ++si) {
		(*si).assign (f.fx(*si));
	}
	stochastic_samples = S.size2();
}


void
 SIR_filter::observe (Likelihood_observe_model& h, const Vec& z)
/*
 * Observation fusion using Likelihood at z
 * Pre : wir previous particle likelihood weights
 * Post: wir fused (multiplicative) particle likehood weights 
 */
{
	h.Lz (z);			// Observe likelihood at z

	{					// Weight Particles. Fused with previous weight
		ColMatrix::iterator si, si_end = S.end();
		Vec::iterator wi = wir.begin();
		for (si = S.begin(); si != si_end; ) {
			*wi *= h.L(*si);
			++wi; ++si;
		}
	}
	wir_update = true;
}

  
SIR_filter::Float
SIR_filter::standard_resample (resamples_t& Presamples, unsigned& uresamples, const ColMatrix& P, Vec& w) const
/*
 * Standard resampler from [1]
 * This algorithm has O(n*log(n)) complexity required to sort the random draws made.
 * This allows comparing of two ordered lists.
 * Sideeffects
 *  w becomes a normalised cumulative sum
 *  A draw is made from this instances 'random' for each particle
 */
{
	assert (P.size2() == Presamples.size());
	assert (P.size2() == w.size());
						// Normalised cumulative sum of likelihood weights, find smallest weight
	Float wcum = 0.;
	Float wmin = std::numeric_limits<Float>::max();
	Vec::iterator wi;
	for (wi = w.begin(); wi != w.end(); ++wi) {
		if (*wi < wmin) {
			wmin = *wi;
		}
		wcum = *wi = wcum + *wi;
	}
	if (wmin < 0.)		// Bad weights
		filter_error("negative weight");
	if (wcum <= 0.)		// Bad cumulative weights (previous check should actually prevent -ve
		filter_error("total likelihood zero");
						// Any numerical failure should cascade into cummulative sum
	if (wcum != wcum)
		filter_error("total likelihood numerical error");

						// Sorted uniform random distribution [0..1) for each resample
	Vec ur(w.size());
	random.uniform_01(ur);
	std::sort (ur.begin(), ur.end());
	assert(ur[0] >= 0. && ur[ur.size()-1] < 1.);		// Very bad if random is incorrect
						// Scale ur to cummulative sum
	ur *= wcum;

						// Resamples based on cumulative weights from sorted resample random values
	ColMatrix::const_iterator pi = P.begin(), pi_end = P.end();
	resamples_t::iterator pri = Presamples.begin();
	wi = w.begin();
	Vec::iterator ui = ur.begin(), ui_end = ur.end();
	unsigned unique = 0;

	while (pi != pi_end)
	{
		unsigned Pres = 0;			// assume P not resampled until find out otherwise
		assert(wi != w.end());
		if (ui != ui_end && *ui < *wi)
		{
			++unique;
			do						// count resamples
			{
				++Pres;
				++ui;
			} while (ui != ui_end && *ui < *wi);
		}
		++wi; ++pi;
		*pri++ = Pres;
	}
	assert(ui==ui_end); assert(wi==w.end());

	uresamples = unique;
	return wmin / wcum;
}


SIR_filter::Float
SIR_filter::systematic_resample (resamples_t& Presamples, unsigned& uresamples, const ColMatrix& P, Vec& w) const
/*
 * Systematic resample algorithm from [2]
 * This algorithm has O(n) complexity
 * Particles are chosesn from those whose cumulative weight intersect with an equidistant grid
 * A uniform random draw is chosen to position the grid with the cumulative weights
 * Sideeffects
 *  w becomes a normalised cumulative sum
 *  A simple draw is made from this instances 'random'
 */
{
	unsigned nParticles = P.size2();
	assert (nParticles == Presamples.size());
	assert (nParticles == w.size());
						// Normalised cumulative sum of likelihood weights, find smallest weight
	Float wcum = 0.;
	Float wmin = std::numeric_limits<Float>::max();
	Vec::iterator wi;
	for (wi = w.begin(); wi != w.end(); ++wi) {
		if (*wi < wmin) {
			wmin = *wi;
		}
		wcum = *wi = wcum + *wi;
	}
	if (wmin < 0.)		// Bad weights
		filter_error("negative weight");
	if (wcum <= 0.)		// Bad cumulative weights (previous check should actually prevent -ve
		filter_error("total likelihood zero");
						// Any numerical failure should cascade into cummulative sum
	if (wcum != wcum)
		filter_error("total likelihood numerical error");

						// Stratified step
	Float wstep = wcum / Float(nParticles);
						
	Vec ur(1);			// Single uniform for initialisation
	random.uniform_01(ur);
	assert(ur[0] >= 0. && ur[0] < 1.);		// Very bad if random is incorrect

						// Resamples based on cumulative weights
	ColMatrix::const_iterator pi = P.begin(), pi_end = P.end();
	resamples_t::iterator pri = Presamples.begin();
	wi = w.begin();
	unsigned unique = 0;
	Float s = ur[0] * wstep;	// Random initialisation

	while (pi != pi_end)
	{
		unsigned Pres = 0;			// assume P not resampled until find out otherwise
		if (s < *wi)
		{
			++unique;
			do						// count resamples
			{
				++Pres;
				s += wstep;
			} while (s < *wi);
		}
		++wi; ++pi;
		*pri++ = Pres;
	}
	assert(wi==w.end());

	uresamples = unique;
	return wmin / wcum;
}


void SIR_filter::copy_resamples (ColMatrix& P, const resamples_t& Presamples)
/*
 * Update P by selectively copying Presamples 
 * Uses a in-place copying alogorithm
 * Algorithm: In-place copying
 *  First copy the live samples (those resampled) to end of P
 *  Replicate live sample in-place
 */
{
							// reverse_copy_if live
	ColMatrix::reverse_iterator liveri = P.rbegin(), sri = P.rbegin();
	resamples_t::const_reverse_iterator pri, pri_end = Presamples.rend();
	for (pri = Presamples.rbegin(); pri != pri_end; ++pri) {
		if (*pri > 0) {
			(*liveri).assign (*sri);
			++liveri;
		}
		++sri;
	}
	assert(sri == P.rend());
							// Replicate live samples
	ColMatrix::iterator si = P.begin();
							// ISSUE livei = liveri; should be used here but is not allowed for MTL reverse_iterators
	ColMatrix::iterator livei = si + (liveri.index()+1);
	resamples_t::const_iterator pi, pi_end = Presamples.end();
	for (pi = Presamples.begin(); pi != pi_end; ++pi) {
		unsigned res = *pi;
		if (res > 0) {
			do  {
				(*si).assign (*livei);
				++si; --res;
			} while (res > 0);
			++livei;
		}
	}
	assert(si == P.end()); assert(livei == P.end());
}


void SIR_filter::roughen_minmax (ColMatrix& P, Float K) const
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
	Float SigmaScale = K * pow (Float(P.size2()), -1./Float(x_size));

						// Find min and max states in all P, precond P not empty
	ColMatrix::iterator pi = P.begin();
	Vec xmin(x_size); xmin = *pi;
	Vec xmax(x_size); xmax = *pi;
	while (pi != P.end())
	{
		Vec::iterator ni = xmin.begin();
		Vec::iterator xi = xmax.begin();

		for (ColMatrix::Column::iterator xpi = (*pi).begin(); xpi != (*pi).end(); ++xpi)
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
	rootq.minus_assign (xmin);
	rootq *= SigmaScale;
   						// Apply roughening prediction based on scaled variance
	Vec n(x_size);
	for (ColMatrix::iterator pi = P.begin(); pi != P.end(); ++pi) {
		random.normal (n);			// independant zero mean normal
									// multiply elements by std dev
		for (FM::Vec::iterator ni = n.begin(); ni != n.end(); ++ni) {
			*ni *= rootq[ni.index()];
		}

		n.plus_assign (*pi);					// add to P
		*pi = n;
	}
}


/*
 * SIR implementation of a Kalman filter
 */
SIR_kalman_filter::SIR_kalman_filter (Subscript x_size, Subscript s_size, Random& random_helper) :
	SIR_filter(x_size, s_size, random_helper), Kalman_filter(x_size),
	rough_random(random_helper),
	rough(x_size,x_size, rough_random)
{
	FM::identity (rough.Fx);
}

void SIR_kalman_filter::init ()
/*
 * Initialise sampling
 *		Post: S
 */
{
						// Samples at mean
	for (ColMatrix::iterator si = S.begin(); si != S.end(); ++si) {
		*si = x;
	}
						// Decorrelate init state noise
	Matrix UD(x_size,x_size);
	Float rcond = UdUfactor (UD, X);
	rclimit.check_PSD(rcond, "Roughening X not PSD");

						// Use predict roughening to apply initial noise
	UdUseperate (rough.G, rough.q, UD);
	rough.init_GqG();
	predict (rough);

	SIR_filter::init(S);
}


void SIR_kalman_filter::mean ()
/*
 * Update state mean
 *		Pre : S
 *		Post: x,S
 */
{
						// Mean of distribution: mean of particles
	x.clear();
	for (ColMatrix::iterator Si = S.begin(); Si != S.end(); ++Si) {
		x.plus_assign (*Si);
	}
	x /= Float(S.size2());
}

void SIR_kalman_filter::covariance ()
/*
 * Update state mean and covariance
 *		Pre : S
 *		Post: x,X,S
 * See numerical note on update()
 */
{
	mean ();
	X.clear();				// Covariance
	const Subscript s_size = S.size2();
	for (ColMatrix::iterator Si = S.begin(); Si != S.end(); ++Si) {
		X.plus_assign (FM::outer_prod(*Si, *Si));
	}
	X /= Float(s_size-1);
	X.minus_assign (FM::outer_prod(x, x));
}


void SIR_kalman_filter::update (Float& lcond)
/*
 * Update mean and covariance of sampled distribution => mean and covariance of particles
 *		Pre : S (s_size >=1 enforced by state_filter base)
 *		Post: x,X,S	(X may be non PSD)
 *
 * Sample Covariance := (1/(s_size-1)* Sum_i [transpose(S[i])*S[i] - transpose(mean)*mean]
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
	SIR_filter::update(lcond);	// Resample particles

	covariance();			// Estimate sample mean and covariance

	// No assert_isPSD (X) as it may fail due to normal numerical circumstances
}


void SIR_kalman_filter::roughen_correlated (ColMatrix& P, Float K)
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
	Float VarScale = sqr(K) * pow (Float(P.size2()), -2./Float(x_size));

	covariance();		// Estimate sample mean and covariance

						// Decorrelate states
	Matrix UD(x_size,x_size);
	Float rcond = UdUfactor (UD, X);
	rclimit.check_PSD(rcond, "Roughening X not PSD");

	UdUseperate (rough.G, rough.q, UD);
						// Apply roughening prediction based on scaled variance
	rough.q *= VarScale;
	rough.init_GqG();

						// Predict particles P using rough model
	for (ColMatrix::iterator pi = P.begin(); pi != P.end(); ++pi) {
		*pi = rough.fx(*pi);
	}
}


}//namespace
