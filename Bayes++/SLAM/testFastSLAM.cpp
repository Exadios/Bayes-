/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.htm for terms and conditions of use.
 *
 * $Id$
 */

/*
 * Test the FastSLAM algorithm
 */

		// Bayes++ Bayesian filtering schemes
#include "BayesFilter/SIRFlt.hpp"
#include "BayesFilter/covFlt.hpp"
#include "BayesFilter/unsFlt.hpp"
#include "BayesFilter/models.hpp"
		// Types required for SLAM classes
#include <vector>
#include <map>
		// Bayes++ SLAM
#include "SLAM.hpp"
#include "fastSLAM.hpp"
#include "kalmanSLAM.hpp"

#include "Test/random.hpp"
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lexical_cast.hpp> 

using namespace SLAM_filter;


class SLAM_random : public Bayesian_filter_test::Boost_random, public BF::SIR_random
/*
 * Random numbers for SLAM test
 */
{
public:
	FM::Float normal (const FM::Float mean, const FM::Float sigma)
	{
		return Boost_random::normal (mean, sigma);
	}
	void normal (FM::DenseVec& v)
	{
		Boost_random::normal (v);
	}
	void uniform_01 (FM::DenseVec& v)
	{
		Boost_random::uniform_01 (v);
	}
	void seed ()
	{
		Boost_random::seed();
	}
};


/*
 * Demonstrate a SLAM example
 */
struct SLAMDemo
{
	const unsigned nParticles;
	
	SLAMDemo (unsigned setnParticles) : nParticles(setnParticles)
	{}
	void OneDExperiment ();
	void InformationLossExperiment ();
	
	SLAM_random goodRandom;

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
	struct Simple_observe_inverse : BF::Linear_uncorrelated_observe_model
	{
		Simple_observe_inverse (Float i_Zv) : Linear_uncorrelated_observe_model(2,1)
		{
			Hx(0,0) = 1.;	// location
			Hx(0,1) = 1.;	// observation
			Zv[0] = i_Zv;
		}
	};

	struct Kalman_statistics : public BF::Kalman_state_filter
	// Kalman_statistics without any filtering
	{
		Kalman_statistics (std::size_t x_size) : Kalman_state_filter(x_size) {}
		void init() {}
		void update() {}
	};

	template <class Filter>
	struct Generic_kalman_generator : public Kalman_filter_generator
	// Generate and dispose of generic kalman filter type
	{
		Filter_type* generate( unsigned full_size )
		{
			return new Filter(full_size);
		}
		void dispose( Filter_type* filter )
		{
			delete filter;
		}
	};


	void display( const std::string label, const BF::Kalman_state_filter& stats)
	{
		std::cout << label << stats.x << stats.X << std::endl;
	}
};


void SLAMDemo::OneDExperiment ()
// Experiment with a one dimensional problem
//  Use to look at implication of highly correlated features
{
	// State size
	const unsigned nL = 1;	// Location
	const unsigned nM = 2;	// Map

	// Construct simple Prediction models
	BF::Sampled_LiAd_predict_model location_predict(nL,1, goodRandom);
	// Stationary Prediction model (Identity)
	FM::identity(location_predict.Fx);
				// Constant Noise model
	location_predict.q[0] = 1000.;
	location_predict.G.clear();
	location_predict.G(0,0) = 1.;

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
	true0.sub_range(0,nL) = x_init; true0[nL] = 50.;
	true1.sub_range(0,nL) = x_init; true1[nL] = 70.;
	FM::Vec z(1);

	// Filter statistics for display
	Kalman_statistics stat(nL+nM);

	// Kalman_SLAM filter:
	Generic_kalman_generator<BF::Covariance_scheme> full_gen;
	Kalman_SLAM kalm (full_gen);
	kalm.init_kalman (x_init, X_init);

	// Fast_SLAM filter
	BF::SIR_kalman_scheme fast_location (nL, nParticles, goodRandom);
	fast_location.init_kalman (x_init, X_init);
	Fast_SLAM_Kstatistics fast (fast_location);

	// Initial feature states
	z = observe0.h(true0);		// Observe a relative position between location and map landmark
	z[0] += 0.5;			
	kalm.observe_new (0, observe_new0, z);
	fast.observe_new (0, observe_new0, z);

	z = observe1.h(true1);
	z[0] += -1.0;		
	kalm.observe_new (1, observe_new1, z);
	fast.observe_new (1, observe_new1, z);

	fast.update(); fast.statistics_sparse(stat); display("Feature Fast", stat);
	kalm.update(); kalm.statistics_sparse(stat); display("Feature Kalm", stat);

	// Predict the location state forward
	fast_location.predict (location_predict);
	kalm.predict (location_predict);
	fast.update(); fast.statistics_sparse(stat); display("Predict Fast", stat);
	kalm.update(); kalm.statistics_sparse(stat); display("Predict Kalm", stat);

	// Observation feature 0
	z = observe0.h(true0);
	z[0] += 0.5;			// Observe a relative position between location and map landmark
	fast.observe( 0, observe0, z );
	kalm.observe( 0, observe0, z );
	fast.update(); fast.statistics_sparse(stat); display("ObserveA Fast", stat);
	kalm.update(); kalm.statistics_sparse(stat); display("ObserveA Kalm", stat);

	// Observation feature 1
	z = observe1.h(true1);
	z[0] += 1.0;			// Observe a relative position between location and map landmark
	fast.observe( 1, observe1, z );
	kalm.observe( 1, observe1, z );
	fast.update(); fast.statistics_sparse(stat); display("ObserveB Fast", stat);
	kalm.update(); kalm.statistics_sparse(stat); display("ObserveB Kalm", stat);

	// Observation feature 0
	z = observe0.h(true0);
	z[0] += 0.5;			// Observe a relative position between location and map landmark
	fast.observe( 0, observe0, z );
	kalm.observe( 0, observe0, z );
	fast.update(); fast.statistics_sparse(stat); display("ObserveC Fast", stat);
	kalm.update(); kalm.statistics_sparse(stat); display("ObserveC Kalm", stat);

	// Forget feature 0
	fast.forget(0);
	kalm.forget(0);
	fast.update(); fast.statistics_sparse(stat); display("Forget Fast", stat);
	kalm.update(); kalm.statistics_sparse(stat); display("Forget Kalm", stat);
}


void SLAMDemo::InformationLossExperiment ()
// Experiment with information loss due to resampling
{
	// State size
	const unsigned nL = 1;	// Location
	const unsigned nM = 2;	// Map

	// Construct simple Prediction models
	BF::Sampled_LiAd_predict_model location_predict(nL,1, goodRandom);
	// Stationary Prediction model (Identity)
	FM::identity(location_predict.Fx);
				// Constant Noise model
	location_predict.q[0] = 1000.;
	location_predict.G.clear();
	location_predict.G(0,0) = 1.;

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
	true0.sub_range(0,nL) = x_init; true0[nL] = 50.;
	true1.sub_range(0,nL) = x_init; true1[nL] = 70.;
	FM::Vec z(1);

	// Filter statistics for display
	Kalman_statistics stat(nL+nM);

	// Kalman_SLAM filter
	Generic_kalman_generator<BF::Unscented_scheme> full_gen;
	Kalman_SLAM kalm (full_gen);
	kalm.init_kalman (x_init, X_init);

	// Fast_SLAM filter
	BF::SIR_kalman_scheme fast_location (nL, nParticles, goodRandom);
	fast_location.init_kalman (x_init, X_init);
	Fast_SLAM_Kstatistics fast (fast_location);

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
			kalm.predict (location_predict);
			fast_location.predict (location_predict);

			// Observation feature 0 with bias
			z = observe0.h(true0);		// Observe a relative position between location and map landmark
			z[0] += 0.5;
			kalm.observe( 0, observe0, z );
			fast.observe( 0, observe0, z );

			// Predict the filter forward
			kalm.predict (location_predict);
			fast_location.predict (location_predict);

			// Observation feature 1 with bias
			z = observe1.h(true1);		// Observe a relative position between location and map landmark
			z[0] += -1.0;		
			kalm.observe( 1, observe1, z );
			fast.observe( 1, observe1, z );
		}

		// Update and resample
		kalm.update();
		fast.update();
		
		kalm.statistics_sparse(stat); display("Kalm", stat);
		fast.statistics_sparse(stat); display("Fast", stat);
		std::cout << fast_location.stochastic_samples <<','<< fast_location.unique_samples()
			<<' '<< fast.feature_unique_samples(0) <<','<< fast.feature_unique_samples(1) <<std::endl;
		std::cout.flush();
	}
}


int main (int argc, char* argv[])
{
	// Global setup for test output
	std::cout.flags(std::ios::fixed); std::cout.precision(4);

	unsigned nParticles = 1000;
	if (argv[1])
	{
    	try {
			nParticles = boost::lexical_cast<unsigned>(argv[1]);
		}
		catch (boost::bad_lexical_cast) {
			// ignore error and use default
		}
	}
	std::cout << "nParticles = " << nParticles << std::endl;

	// Create test and run experiments
	try {
		SLAMDemo test(nParticles);
		test.OneDExperiment();
		//test.InformationLossExperiment();
	}
	catch (const BF::Filter_exception& ne)
	{
		std::cout << ne.what() << std::endl;
	}
}
