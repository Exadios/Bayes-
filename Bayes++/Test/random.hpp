/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 */

/*
 * Good random numbers from Boost
 *  Provides a common class  for all random number requirements to test Bayes++
 */

#include <boost/version.hpp>
#include <boost/random.hpp>


namespace Bayesian_filter_test
{

#if (BOOST_VERSION >= 103100)

class Boost_random
/*
 * Random numbers from Boost 1_31_0
 */
{
public:
	Boost_random() : dist_normal(), dist_uniform()
	{}
	Bayesian_filter_matrix::Float normal(const Bayesian_filter_matrix::Float mean, const Bayesian_filter_matrix::Float sigma)
	{
		boost::normal_distribution<Bayesian_filter_matrix::Float> dist(mean, sigma);
		boost::variate_generator<boost::mt19937& ,boost::normal_distribution<Bayesian_filter_matrix::Float> > gen(rng, dist);
		return gen();
	}
	void normal(Bayesian_filter_matrix::DenseVec& v, const Bayesian_filter_matrix::Float mean, const Bayesian_filter_matrix::Float sigma)
	{
		boost::normal_distribution<Bayesian_filter_matrix::Float> dist(mean, sigma);
		boost::variate_generator<boost::mt19937&, boost::normal_distribution<Bayesian_filter_matrix::Float> > gen(rng, dist);
		std::generate (v.begin(), v.end(), gen);
	}
	void normal(Bayesian_filter_matrix::DenseVec& v)
	{
		boost::variate_generator<boost::mt19937& ,boost::normal_distribution<Bayesian_filter_matrix::Float> > gen(rng, dist_normal);
		std::generate (v.begin(), v.end(), gen);
	}
	void uniform_01(Bayesian_filter_matrix::DenseVec& v)
	{
		boost::variate_generator<boost::mt19937& ,boost::uniform_real<Bayesian_filter_matrix::Float> > gen(rng, dist_uniform);
		std::generate (v.begin(), v.end(), gen);
	}
	void seed()
	{
		rng.seed();
	}
private:
	boost::mt19937 rng;
	boost::normal_distribution<Bayesian_filter_matrix::Float> dist_normal;
	boost::uniform_real<Bayesian_filter_matrix::Float> dist_uniform;
};

#else

// Random numbers from Boost 1_30_0 or earlier
class Boost_random : public BF::SIR_random, public BF::General_LiInAd_predict_model::Random
/*
 * Random numbers from Boost 1_30_0 or earlier
 */
{
public:
	Boost_random() : gen_normal(rng), gen_uniform(rng)
	{}
	double normal(const double mean, const double sigma)
	{
		boost::normal_distribution<boost::mt19937,FM::Float> gen(rng, mean, sigma);
		return gen();
	}
	void normal(Vec& v, const FM::Float mean, const FM::Float sigma)
	{
		boost::normal_distribution<boost::mt19937,FM::Float> gen(rng, mean, sigma);
		std::generate (v.begin(), v.end(), gen);
	}
	void normal(FM::Vec& v)
	{
		std::generate (v.begin(), v.end(), gen_normal);
	}
	void uniform_01(FM::Vec& v)
	{
		std::generate (v.begin(), v.end(), gen_uniform);
	}
	void reseed()
	{
		rng.seed();
	}
private:
	boost::mt19937 rng;
	boost::normal_distribution<boost::mt19937,FM::Float> gen_normal;
	boost::uniform_01<boost::mt19937,FM::Float> gen_uniform;
};

#endif

}//namespace
