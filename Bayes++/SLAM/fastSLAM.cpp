/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 */

/*
 * SLAM : Simultaneous Locatization and Mapping
 *  FastSLAM augmented particle algorithm
 *  Direct implementation without log(n) tree pruneing for associated observations
 */

		// Bayes++ Bayesian filtering schemes
#include "BayesFilter/SIRFlt.hpp"
		// Types required for SLAM classes
#include <map>
		// Bayes++ SLAM
#include "SLAM.hpp"
#include "fastSLAM.hpp"
#include <cmath>


namespace {
	template <class scalar>
	inline scalar sqr(scalar x)
	{
		return x*x;
	}
}

namespace SLAM_filter
{

Fast_SLAM::Fast_SLAM( BF::SIR_scheme& L_filter ) :
	SLAM(),
	L(L_filter),
	wir(L.S.size2())
// Construct filter using referenced SIR_filter for resampling
{
	std::fill (wir.begin(), wir.end(), Float(1.));		// Initial uniform weights
	wir_update = false;
}

void Fast_SLAM::observe_new( unsigned feature, const Feature_observe_inverse& fom, const FM::Vec& z )
/*
 * SLAM New Feature observation (overwrite)
 * Assumes there is no prior information about the feature (strictly a uniform un-informative prior)
 *
 * This implies
 *  a) There has no information about location so no resampling is requires
 *  b) Feature Posterior estimated directly from observation using inormation form
 *
 * fom: must have a the special from required for SLAM::obeserve_new
 */
{
	assert(z.size() == 1);		// Only single state observation supported

	const size_t nL = L.S.size1();	// No of location states
	const size_t nparticles = L.S.size2();
	FeatureCondMap fmap(nparticles);

	FM::Vec sz(nL+1);						// Location state

	for (size_t pi = 0; pi < nparticles; ++pi)
	{
		sz.sub_range(0,nL) = L.S.column(pi);
		sz[nL] = z[0];
		FM::Vec t = fom.h(sz);
		fmap[pi].x = t[0];
		fmap[pi].X = fom.Zv[0];
	}
	M.insert (std::make_pair(feature, fmap));
}

void Fast_SLAM::observe_new( unsigned feature, const FM::Float& t, const FM::Float& T )
/*
 * SLAM New observation directly of state statistics (overwrite)
 */
{
	const size_t nparticles = L.S.size2();
	Feature_1 m1;  				// single map feature
	FeatureCondMap fmap(nparticles);

	m1.x = t;			   // Initial particle conditional map is sample
	m1.X = T;		    // independant
	std::fill(fmap.begin(),fmap.end(), m1);
	M.insert (std::make_pair(feature, fmap));
}

void Fast_SLAM::observe( unsigned feature, const Feature_observe& fom, const FM::Vec& z )
/*
 * SLAM Feature observation
 *  Uses Extended Fast_SLAM observation equations
 */
{
	assert(z.size() == 1);		// Only single state observation supported

	const AllFeature::iterator inmap = M.find(feature);
	if (inmap == M.end())
	{
		error (BF::Logic_exception("Observe non existing feature"));
		return;
	}
								// Existing feature
	FeatureCondMap& afm = (*inmap).second;	// Reference the associated feature map
	const size_t nL = L.S.size1();	// No of location states
	const size_t nparticles = L.S.size2();

	Float Ht = fom.Hx(0,nL);
	if (Ht == 0.)
		error (BF::Numeric_exception("observe Hx feature component zero"));

								// Observation weight for each particle
	FM::Vec x2(nL+1);					// Augmented state (particle + feature mean)
	FM::Vec znorm(z.size());
	for (size_t pi = 0; pi < nparticles; ++pi)
	{
		Feature_1& m1 = afm[pi];		// Associated feature's map particle
							
		x2.sub_range(0,nL) = L.S.column(pi);		// Build Augmented state x2
		x2[nL] = m1.x;
		const FM::Vec& zp = fom.h(x2);	// Observation model
		znorm = z;						// Normalised observation
		fom.normalise(znorm, zp);

		const Float zpp = zp[0] - Ht*m1.x;	// Extended State observation for non-linear h
		const Float p = (znorm[0] - zpp) / Ht;
		const Float P = fom.Zv[0] /sqr(Ht) ;
		const Float q = m1.x;
		const Float Q = m1.X;
											// Multiplicative fussion of observation weights
 		wir[pi] *= exp(Float(-0.5) *sqr(p-q) / (P+Q)) / sqrt(P+Q);		// Integral( g(p,P)*g(q,Q) )
	}
	wir_update = true;			// Weights have been updated requiring a resampling

								// Estimate associated features conditional map for resampled particles
	for (size_t pi = 0; pi < nparticles; ++pi)
	{
		Feature_1& m1 = afm[pi];		// Associated feature's map particle
		x2.sub_range(0,nL) = L.S.column(pi);		// Build Augmented state x2
		x2[nL] = m1.x;
		FM::Vec zp = fom.h(x2);			// Observation model
		znorm = z;						// Normalised observation
		fom.normalise(znorm, zp);

								// EKF for conditional feature observation - specialised for 1D and zero state uncertianty
		Float Sf = m1.X*sqr(Ht) + fom.Zv[0];	// Innovation variance
		if (Sf <= 0.)
			error (BF::Numeric_exception("Conditional feature estimate not PD"));
		Float Wf = m1.X*Ht / Sf;
		
		m1.x += Wf*(znorm[0] - zp[0]);
		m1.X -= sqr(Wf) * Sf;
	}
}

Fast_SLAM::Float
 Fast_SLAM::update_resample( const Bayesian_filter::Importance_resampler& resampler )
/* Resampling Update
 *  Resample particles using weights
 *  Propogate resampling to All features
 *  Only resamples if weights have been updated
 */
{
	if (wir_update)
	{
		const size_t nparticles = L.S.size2();
		Resamples_t presamples(nparticles);
		size_t R_unique;			// Determine resamples of S
		Float lcond = resampler.resample (presamples, R_unique, wir, L.random);

									// Initial uniform weights
		std::fill (wir.begin(), wir.end(), Float(1.));
		wir_update = false;
									// Update S bases on resampling, and init filter
		L.copy_resamples (L.S, presamples);
		L.init_S ();
									// Propogate resampling to All features
		FeatureCondMap fmr(nparticles);		// Resampled feature map
		for (AllFeature::iterator fi = M.begin(); fi != M.end(); ++fi)	// All Features
		{
			FeatureCondMap& fm = (*fi).second;		// Reference the feature map
										// Iterate over All feature particles
			FeatureCondMap::iterator fmi, fmi_begin = fm.begin(), fmi_end = fm.end();
			FeatureCondMap::iterator fmri = fmr.begin();
			for (fmi = fmi_begin; fmi < fmi_end; ++fmi)
			{							// Multiple copies of this resampled feature
				for (size_t res = presamples[fmi-fmi_begin]; res > 0; --res) {
					*fmri = *fmi;
					++fmri;
				}
			}
			fm = fmr;				// Copy in resamples feature map
		}

		L.roughen ();				// Roughen location
		L.stochastic_samples = R_unique;
		return lcond;
	}
	else
		return 1.;		// No resampling
}

void Fast_SLAM::forget( unsigned feature, bool must_exist )
// Forget all feature information, feature no can be reused for a new feature
{
	AllFeature::size_type n = M.erase(feature);
	if (n==0 && must_exist)
		error (BF::Logic_exception("Forget non existing feature"));
}

size_t Fast_SLAM::feature_unique_samples( unsigned feature )
/*
 * Count the number of unique samples in S associated with a feature
 */
{
	const AllFeature::iterator inmap = M.find(feature);
	if (inmap == M.end())
	{
		error (BF::Logic_exception("feature_unique_samples non existing feature"));
		return 0;
	}
								// Existing feature
	FeatureCondMap& afm = (*inmap).second;	// Reference the associated feature map

	typedef FeatureCondMap::iterator Sref;
	// Provide a ordering on feature sample means
	struct order {
		static bool less(Sref a, Sref b)
		{
			return (*a).x < (*b).x;
		}
	};

						// Sorted reference container
	typedef std::vector<Sref> SRContainer;
	SRContainer sortR(afm.size());

						// Reference each element in S
	{	Sref elem = afm.begin();
		SRContainer::iterator ssi = sortR.begin();
		for (; ssi < sortR.end(); ++ssi)
		{
			*ssi = elem; ++elem;
		}
	}

	std::sort (sortR.begin(), sortR.end(), order::less);

						// Count element changes, precond: sortS not empty
	size_t u = 1;
	SRContainer::const_iterator ssi= sortR.begin();
	SRContainer::const_iterator ssp = ssi;
	++ssi;
	while (ssi < sortR.end())
	{
		if (order::less(*ssp, *ssi))
			++u;
		ssp = ssi;
		++ssi;
	}
	return u;
}


/*
 * Fast_SLAM_Kstatistics
 */
Fast_SLAM_Kstatistics::Fast_SLAM_Kstatistics( BF::SIR_kalman_scheme& L_filter ) :
	Fast_SLAM(L_filter), L(L_filter)
// Construct filter using referenced SIR_filter for resampling
{
}


void Fast_SLAM_Kstatistics::statistics_feature(
		BF::Kalman_state_filter& kstat, size_t fs,
		const AllFeature::const_iterator& fi, const AllFeature::const_iterator& fend ) const
/*
 * Compute sample mean and covariance statistics of feature.
 *  fs is subscript in kstat to return statisics
 *  fi iterator of feature
 *  fend end of map iterator (statistics are computed for fi with a map subset)
 * Precond: kstat contains location mean
 */
{
	const size_t nL = L.S.size1();	// No of location states
	const size_t nparticles = L.S.size2();

	const FeatureCondMap& fm = (*fi).second;		// Reference the feature map

									// Iterate over All feature particles
	FeatureCondMap::const_iterator fpi, fpi_begin = fm.begin(), fpi_end = fm.end();

	Float mean = 0.;				// Feature mean
	for (fpi = fpi_begin; fpi < fpi_end; ++fpi)	{
		mean += (*fpi).x;
	}
	mean /= Float(nparticles);

	Float var = 0.;					// Feature variance: 1/(n-1) is unbiased estimate given estimated mean
	for (fpi = fpi_begin; fpi < fpi_end; ++fpi) {
		var += (*fpi).X + sqr((*fpi).x);
	}
	var = (var - Float(nparticles)*sqr(mean)) / Float(nparticles-1);

	kstat.x[fs] = mean;				// Copy into Kalman statistics
	kstat.X(fs,fs) = var;

									// Location,feature covariance
	for (size_t si=0; si < nL; ++si)
	{
		Float covar = 0.;
		size_t spi=0;
		for (fpi = fpi_begin; fpi < fpi_end; ++fpi) {
			covar += (*fpi).x*L.S(si,spi);
			++spi;
		}
		covar = (covar - Float(nparticles)*mean*kstat.x[si]) / Float(nparticles);
		kstat.X(si,fs) = covar;
	}

									// Feature,feature covariance. Iterate over previous features with means already computed
	size_t fsj = nL;				// Feature subscript
	for (AllFeature::const_iterator fj = M.begin(); fj != fend; ++fj, ++fsj)
	{
		Float covar = 0.;
		FeatureCondMap::const_iterator fpj = (*fj).second.begin();
		for (fpi = fpi_begin; fpi < fpi_end; ++fpi) {
			covar += (*fpi).x*(*fpj).x;
			++fpj;
		}
		covar = (covar - Float(nparticles)*mean*kstat.x[fsj]) / Float(nparticles);
		kstat.X(fs,fsj) = covar;
	}
}


void Fast_SLAM_Kstatistics::statistics_compressed( BF::Kalman_state_filter& kstat )
/*
 * Compute sample mean and covariance statistics of filter
 *  
 *  kstat elements are filled first with Location statistics and then the Map feature statistics
 *  Feature statisics are are computed in feature number order and only for those for which there is space in kstat
 * Note: Covariance values are indeterminate for nparticles ==1
 * Precond:
 *   nparticles >=1 (enforced by Sample_filter contruction)
 *   kstat must have space for Location statistics
 * Postcond:
 *  kstat compressed sample statisics of filter
 */
{	
	const size_t nL = L.S.size1();	// No of location states

	kstat.x.clear();			// Zero everything (required only for non existing feature states
	kstat.X.clear();			// Zero everything (required only for non existing feature states

								// Get Location statistics
	if (nL > kstat.x.size())
		error (BF::Logic_exception("kstat to small to hold filter locatition statistics"));
	L.update_statistics();
	kstat.x.sub_range(0,nL) = L.x;
	kstat.X.sub_matrix(0,nL, 0,nL) = L.X;

								// Iterated over feature statistics (that there is space for in kstat)
	size_t fs = nL;						// Feature subscript
	for (AllFeature::const_iterator fi = M.begin(); fi != M.end() && fs < kstat.x.size(); ++fi, ++fs)
	{
		statistics_feature(kstat, fs, fi, fi);	// build statistics of fi with other features up to fi
	}
}//statistics_compressed


void Fast_SLAM_Kstatistics::statistics_sparse( BF::Kalman_state_filter& kstat )
/*
 * Compute sample mean and covariance statistics of filter
 *  
 *  kstat elements are filled first with Location statistics and then the Map feature statistics
 *  Feature statisics are are computed in feature number as index (after location) and only for those for which there is space in kstat
 * Note: Covariance values are indeterminate for nparticles ==1
 * Precond:
 *   nparticles >=1 (enforced by Sample_filter contruction)
 *   kstat must have space for Location statistics
 * Postcond:
 *  kstat sparse sample statisics of filter
 */
{	
	const size_t nL = L.S.size1();	// No of location states

	kstat.x.clear();			// Zero everything (required only for non existing feature states
	kstat.X.clear();			// Zero everything (required only for non existing feature states

								// Get Location statistics
	if (nL > kstat.x.size())
		error (BF::Logic_exception("kstat to small to hold filter locatition statistics"));
	L.update_statistics();
	kstat.x.sub_range(0,nL) = L.x;
	kstat.X.sub_matrix(0,nL, 0,nL) = L.X;

								// Iterated over feature statistics (that there is space for in kstat)
	for (AllFeature::const_iterator fi = M.begin(); fi != M.end(); ++fi)
	{
		size_t fs = nL + (*fi).first;		// Feature subscript
		if (fs < kstat.x.size())			// Space in kstat
		{
		statistics_feature(kstat, fs, fi, fi);	// build statistics of fi with other features up to fi
		}
	}//all feature
}//statistics_sparse


}//namespace SLAM
