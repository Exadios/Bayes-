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

/*
 * Random numbers from Boost 1_31_0
 */

namespace
{
	// ISSUE variate_generator cannot be used without Partial template specialistion
	template<class Engine, class Distribution>
	class simple_generator
	{
	public:
		typedef typename Distribution::result_type result_type;
		simple_generator(Engine& e, Distribution& d)
			: _eng(e), _dist(d)
		{}
		result_type operator()()
		{	return _dist(_eng); }
	private:
		Engine& _eng;
		Distribution& _dist;
	};
}//namespace

class Boost_random
{
public:
	typedef Bayesian_filter_matrix::Float Float;
	typedef boost::uniform_01<boost::mt19937, Float> UGen;
	Boost_random() : gen01(boost::mt19937()), dist_normal()
	{}
	Bayesian_filter_matrix::Float normal(const Float mean, const Float sigma)
	{
		boost::normal_distribution<Float> dist(mean, sigma);
		simple_generator<UGen, boost::normal_distribution<Float> > gen(gen01, dist);
		return gen();
	}
	void normal(Bayesian_filter_matrix::DenseVec& v, const Float mean, const Float sigma)
	{
		boost::normal_distribution<Float> dist(mean, sigma);
		simple_generator<UGen, boost::normal_distribution<Float> > gen(gen01, dist);
		std::generate (v.begin(), v.end(), gen);
	}
	void normal(Bayesian_filter_matrix::DenseVec& v)
	{
		simple_generator<UGen, boost::normal_distribution<Float> > gen(gen01, dist_normal);
		std::generate (v.begin(), v.end(), gen);
	}
	void uniform_01(Bayesian_filter_matrix::DenseVec& v)
	{
		std::generate (v.begin(), v.end(), gen01);
	}
	void seed()
	{
		gen01.base().seed();
	}
private:
	UGen gen01;
	boost::normal_distribution<Float> dist_normal;
};


#else

/*
 * Random numbers from Boost 1_30_0 or earlier
 */
class Boost_random : public BF::SIR_random, public BF::General_LiInAd_predict_model::Random
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
