/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 */

/*
 * Test the FastSLAM alogorithm
 *  A linear filter with one state and constant noises
 */

// Types Required for SLAM classes
#include <vector>
#include <map>

// Include all the Bayes++ Bayesian filtering library
#include "BayesFilter/allFilters.hpp"
#include "SLAM.hpp"
#include "fastSLAM.hpp"
#include "kalmanSLAM.hpp"

#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>




using namespace SLAM_filter;


// Random numbers for filters from Boost
class Boost_random : public BF::SIR_random, public BF::General_LiInAd_predict_model::Random
/*
 * Random number distributions
 */
{
public:
	Boost_random() : dist_normal(), dist_uniform()
	{}
	FM::Float normal(const FM::Float mean, const FM::Float sigma)
	{
		boost::normal_distribution<FM::Float> dist(mean, sigma);
		boost::variate_generator<boost::mt19937& ,boost::normal_distribution<FM::Float> > gen(rng, dist);
		return gen();
	}
	void normal(FM::DenseVec& v)
	{
		boost::variate_generator<boost::mt19937& ,boost::normal_distribution<FM::Float> > gen(rng, dist_normal);
		std::generate (v.begin(), v.end(), gen);
	}
	void uniform_01(FM::DenseVec& v)
	{
		boost::variate_generator<boost::mt19937& ,boost::uniform_real<FM::Float> > gen(rng, dist_uniform);
		std::generate (v.begin(), v.end(), gen);
	}
	void reseed()
	{
		rng.seed();
	}
private:
	boost::mt19937 rng;
	boost::normal_distribution<FM::Float> dist_normal;
	boost::uniform_real<FM::Float> dist_uniform;
};





/*
 * Demonstrate a SLAM example
 */
struct SLAMDemo
{
	void OneDExperiment ();
	void InformationLossExperiment ();
	
	Boost_random goodRandom;

	// Relative Observation with  Noise model
	struct Simple_observe : BF::Linear_uncorrelated_observe_model
	{
		Simple_observe (Float i_Zv) : Linear_uncorrelated_observe_model(2,1)
		// Construct a linear model with const Hx
		{
			Hx(0,0) = -1.;	// Location
			Hx(0,1) = 1.;	// Map
			Zv[0] = i_Zv;
		}
	};
	struct Simple_observe_inverse : BF::Uncorrelated_addative_observe_model
	{
		Simple_observe_inverse (Float i_Zv) : Uncorrelated_addative_observe_model(1), t(1)
		{
			Zv[0] = i_Zv;
		}
		const FM::Vec& h(const FM::Vec& lz) const
		{
			t[0] = lz[1]+lz[0];
			return t;
		}
		mutable FM::Vec t;
	};

	class Kalman_statistics : public BF::Kalman_state_filter
	// Kalman_statistics without any filtering
	{	public:
		Kalman_statistics (size_t x_size) : Kalman_state_filter(x_size) {}
		void init() {}
		void update() {}
	};

	void display( const std::string label, const BF::Kalman_state_filter& stats)
	{
		std::cout << label << stats.x << stats.X << std::endl;
	}
};


void SLAMDemo::OneDExperiment()
// Experiment with a one dimensional problem
{
	// State size
	const unsigned nL = 1;	// Location
	const unsigned nM = 2;	// Map

	// Construct simple Prediction models
	BF::General_LiAd_predict_model location_predict(nL,1, goodRandom);
	BF::General_LiAd_predict_model all_predict(nL+nM,1, goodRandom);
	// Stationary Prediction model (Identity)
	FM::identity(location_predict.Fx);	FM::identity(all_predict.Fx);
				// Constant Noise model
	location_predict.q[0] = 1000.;	all_predict.q[0] = 1000.;
	location_predict.G.clear(); all_predict.G.clear();
	location_predict.G(0,0) = 1.; all_predict.G(0,0) = 1.;

	// Relative Observation with  Noise model
	Simple_observe observe0(5.), observe1(3.);
	Simple_observe_inverse observe_new0(5.), observe_new1(3.);

	// Setup the initial state and covariance
	// Location with no uncertainty
	FM::Vec x_init(nL); FM::SymMatrix X_init(nL, nL);
	x_init[0] = 20.;
	X_init(0,0) = 0.;

	// Truth model : location plus one map feature
	FM::Vec true0(nL+1), true1(nL+1);
	true0(0,nL) = x_init; true0[nL] = 50.;
	true1(0,nL) = x_init; true1[nL] = 70.;
	FM::Vec z(1);

	// Kalman_SLAM filter
	BF::Unscented_scheme full_filter(nL+nM);
	full_filter.x.clear(); full_filter.X.clear();
	full_filter.x[0] = x_init[0];
	full_filter.X(0,0) = X_init(0,0);
	full_filter.init();
	Kalman_SLAM kalm(full_filter, nL);

	// Fast_SLAM filter
#ifdef NDEBUG
	unsigned nParticles = 1000000;
#else
	unsigned nParticles = 20;
#endif
	BF::SIR_kalman_scheme fast_location(nL, nParticles, goodRandom);
	fast_location.init_kalman(x_init, X_init);
	Fast_SLAM_Kstatistics fast(fast_location);

	Kalman_statistics stat(nL+nM);

	// Initial feature states
	z = observe0.h(true0);		// Observe a relative position between location and map landmark
	z[0] += 0.5;			
	kalm.observe_new (0, observe_new0, z);
	fast.observe_new (0, observe_new0, z);

	z = observe1.h(true1);
	z[0] += -1.0;		
	kalm.observe_new (1, observe_new1, z);
	fast.observe_new (1, observe_new1, z);

	// Experiment with highly correlated features
	//
	kalm.update(); display("Initial Kalm", full_filter);
	fast.update(); fast.statistics(stat); display("Initial Fast", stat);

	// Predict the filter forward
	full_filter.predict (all_predict);
	kalm.update(); display("Predict Kalm", full_filter);
	fast_location.predict (location_predict);
	fast.update(); fast.statistics(stat); display("Predict Fast", stat);

	// Observation feature 0
	z = observe0.h(true0);
	z[0] += 0.5;			// Observe a relative position between location and map landmark
	kalm.observe( 0, observe0, z );
	kalm.update(); display("ObserveA Kalm", full_filter);
	fast.observe( 0, observe0, z );
	fast.update(); fast.statistics(stat); display("ObserveA Fast", stat);

	// Observation feature 1
	z = observe1.h(true1);
	z[0] += 1.0;			// Observe a relative position between location and map landmark
	kalm.observe( 1, observe1, z );
	kalm.update(); display("ObserveB Kalm", full_filter);
	fast.observe( 1, observe1, z );
	fast.update(); fast.statistics(stat); display("ObserveB Fast", stat);

	// Observation feature 0
	z = observe0.h(true0);
	z[0] += 0.5;			// Observe a relative position between location and map landmark
	kalm.observe( 0, observe0, z );
	kalm.update(); display("ObserveC Kalm", full_filter);
	fast.observe( 0, observe0, z );
	fast.update(); fast.statistics(stat); display("ObserveC Fast", stat);

	// Forget feature 0
	kalm.forget(0);
	kalm.update(); display("Forget Kalm", full_filter);
	fast.forget(0);
	fast.update(); fast.statistics(stat); display("Forget Fast", stat);
}


void SLAMDemo::InformationLossExperiment()
// Experiment with information loss due to resampling
{
	// State size
	const unsigned nL = 1;	// Location
	const unsigned nM = 2;	// Map

	// Construct simple Prediction models
	BF::General_LiAd_predict_model location_predict(nL,1, goodRandom);
	BF::General_LiAd_predict_model all_predict(nL+nM,1, goodRandom);
	// Stationary Prediction model (Identity)
	FM::identity(location_predict.Fx);	FM::identity(all_predict.Fx);
				// Constant Noise model
	location_predict.q[0] = 1000.;	all_predict.q[0] = 1000.;
	location_predict.G.clear(); all_predict.G.clear();
	location_predict.G(0,0) = 1.; all_predict.G(0,0) = 1.;

	// Relative Observation with  Noise model
	Simple_observe observe0(5.), observe1(3.);
	Simple_observe_inverse observe_new0(5.), observe_new1(3.);

	// Setup the initial state and covariance
	// Location with no uncertainty
	FM::Vec x_init(nL); FM::SymMatrix X_init(nL, nL);
	x_init[0] = 20.;
	X_init(0,0) = 0.;

	// Truth model : location plus one map feature
	FM::Vec true0(nL+1), true1(nL+1);
	true0(0,nL) = x_init; true0[nL] = 50.;
	true1(0,nL) = x_init; true1[nL] = 70.;
	FM::Vec z(1);

	// Kalman_SLAM filter
	BF::Unscented_scheme full_filter(nL+nM);
	full_filter.x.clear(); full_filter.X.clear();
	full_filter.x[0] = x_init[0];
	full_filter.X(0,0) = X_init(0,0);
	full_filter.init();
	Kalman_SLAM kalm(full_filter, nL);

	// Fast_SLAM filter
	unsigned nParticles = 1000;
	BF::SIR_kalman_scheme fast_location(nL, nParticles, goodRandom);
	fast_location.init_kalman(x_init, X_init);
	Fast_SLAM_Kstatistics fast(fast_location);

	Kalman_statistics stat(nL+nM);

	// Initial feature states
	z = observe0.h(true0);		// Observe a relative position between location and map landmark
	z[0] += 0.5;			
	kalm.observe_new (0, observe_new0, z);
	fast.observe_new (0, observe_new0, z);

	z = observe1.h(true1);
	z[0] += -1.0;		
	kalm.observe_new (1, observe_new1, z);
	fast.observe_new (1, observe_new1, z);

	unsigned it = 0;
	for (;;) {
		++it;
		std::cout << it << std::endl;
		
		// Groups of observations without resampling
		{
			// Predict the filter forward
			full_filter.predict (all_predict);
			fast_location.predict (location_predict);

			// Observation feature 0 with bias
			z = observe0.h(true0);		// Observe a relative position between location and map landmark
			z[0] += 0.5;
			kalm.observe( 0, observe0, z );
			fast.observe( 0, observe0, z );

			// Predict the filter forward
			full_filter.predict (all_predict);
			fast_location.predict (location_predict);

			// Observation feature 1 with bias
			z = observe1.h(true1);		// Observe a relative position between location and map landmark
			z[0] += -1.0;		
			kalm.observe( 1, observe1, z );
			fast.observe( 1, observe1, z );
		}

		// Update and resample
		kalm.update();
		try {
			fast.update();
		}
		catch (BF::Filter_exception ne)
		{
			std::cout << ne.what() << std::endl;
			std::cout.flush();
			throw;	// re-throw
		}
		
		display("Kalm", full_filter);
		fast.statistics(stat); display("Fast", stat);
		std::cout << fast_location.stochastic_samples <<','<< fast_location.unique_samples()
			<<' '<< fast.feature_unique_samples(0) <<','<< fast.feature_unique_samples(1) <<std::endl;
		std::cout.flush();
	}
}


int main()
{
	// Global setup for test output
	std::cout.flags(std::ios::fixed); std::cout.precision(4);

	// Create test and run experiments
	SLAMDemo test;
	test.OneDExperiment();
	//test.InformationLossExperiment();
}
