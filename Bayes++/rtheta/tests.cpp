/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 */
#include "BayesFilter/allFilters.hpp"
#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/limits.hpp>
#include <boost/random.hpp>


using namespace Bayesian_filter;
using namespace Bayesian_filter_matrix;


// Instantiate complete fileters to check the templates
#include "BayesFilter/filters/average1.hpp"
template average1_filter<Covariance_filter>;

#include "BayesFilter/filters/indirect.hpp"
template Indirect_state_filter<State_filter>;
template Indirect_kalman_filter<Kalman_filter>;


// Square
template <class scalar>
inline scalar sqr(scalar x)
{
	return x*x;
}


template <typename Float>
class Test_random : public SIR_random
/*
 * Random number distributions
 */
{
public:
	Test_random() : gen_normal(rng), gen_uniform(rng), gen_exponential(rng, 1.)
	{}
	Float normal(const Float mean, const Float sigma)
	{
		boost::normal_distribution<boost::mt19937, Float> gen(rng, mean, sigma);
		return gen();
	}
	void normal(DenseVec& v)
	{
		std::generate (v.begin(), v.end(), gen_normal);
	}
	void uniform_01(DenseVec& v)
	{
		std::generate (v.begin(), v.end(), gen_uniform);
	}
	void exponential_1(DenseVec& v)
	{
		std::generate (v.begin(), v.end(), gen_exponential);
	}
	void reseed()
	{
		rng.seed();
	}
private:
	typedef boost::mt19937 good_random;
	good_random rng;
	boost::normal_distribution<good_random,Float> gen_normal;
	boost::uniform_01<good_random,Float> gen_uniform;
	boost::exponential_distribution<good_random,Float> gen_exponential;
};

void test_inverse()
{
	Matrix U (3,3);
	U.clear();
	U(0,0) = 1*1;
	U(1,1) = 1*1;
	U(2,2) = 1*1;
	U(0,1) = 2;
	U(0,2) = 3;
	U(1,2) = 5;
/*	
	Vec rv(U.size1()*U.size2());
	Test_random rand;
	rand.exponential_1(rv);

	Vec::iterator rv_i = rv.begin();
	for (size_t r = 0; r < U.size1(); ++r)
		for (size_t c = 0; c < U.size2(); ++c)
		{
			U(r,c) = *rv_i; ++rv_i;
		}
	*/
	std::cout << U << std::endl;

	SymMatrix X(U.size1(),U.size2());
	X = prod(U, trans(U));
	std::cout << X << std::endl;
	
	SymMatrix XI(X.size1(), X.size2());
	SymMatrix XII(X.size1(), X.size2());
//	XI = X;
//	UdUfactor(XI, 2);
//	std::cout << XI << std::endl;
//	XI = X;
//	LdLfactor(XI, 2);
//	std::cout << XI << std::endl;

	UdUinversePD (XI, X);
	std::cout << XI << std::endl;
	std::cout << prod(X,XI) << std::endl;

	UdUinversePD (XII, XI);
	std::cout << XII << std::endl;
};

void test_SPD()
{
	SymMatrix S(1,1);
	S(0,0) = 5.;
	Vec Stemp(1);
	SymMatrix P(2,2);

	RowMatrix R(2,1);
	R(0,0) = 1.;
	R(1,0) = 2.;

	// Try prod_SPD
	RowMatrix RStemp(R.size1(), S.size2());
	P.assign( prod_SPD(R,S, RStemp) );
	std::cout << P << std::endl;

	// Try prod_SPD
	ColMatrix RT(1,2);
	P.assign (prod_SPDT (RT, Stemp));
}

void test_SPD_all()
{
	SymMatrix S(1,1);
	Vec x(1);
	S(0,0) = 5.;
	x[0] = 1.;

	Float v1 = ublas::inner_prod (x, ublas::prod(S,x) );
	Float v2 = prod_SPDT(x,S);
	(void)v1,(void)v2;

	v2 = prod_SPD(x,x);

	test_SPD();
}

void test_temp_prod()
{
	Matrix A(2,2);
	Matrix R (prod( ublas::prod<Matrix>(A,A), ublas::prod<Matrix>(A,A)));
}


void test_UdU()
{
	SymMatrix S(3,3);
	S(0,0) = 2.;
	S(2,2) = 3.;
	S(0,1) = S(1,0) = 5.;
//	S(2,1) = S(1,2) = 5.;
	S(0,2) = S(2,0) = 7.;
	std::cout << S << std::endl;

	UTriMatrix UD(3,3);
	Float rcond = UCfactor (UD, S);
	std::cout <<rcond << std::endl;

	std::cout << UD << std::endl;

//	UdUrecompose (UD);
//	std::cout << UD << std::endl;
}

void test_ident()
{
	SymMatrix S(3,3);
	FM::identity(S);
	std::cout << S << std::endl;
}

void test_info_init()
{
	Information_root_info_filter iri(2);
	Information_filter i(2);
	SymMatrix X(2,2);
	X(0,0) = 2.;
	X(1,0) = X(0,1) = 1.1;
	X(1,1) = 3.;
	Vec x(2);
	x[0] = 5.;
	x[1] = 7.;

	i.init_kalman (x,X);
	iri.init_kalman (x,X);
	i.update(); iri.update();
	std::cout << i.x << i.X <<std::endl;;
	std::cout << i.y << i.Y <<std::endl;;
	std::cout << iri.x << iri.X <<std::endl;;
	std::cout << iri.y << iri.Y <<std::endl;;

	i.init_information (x,X);
	iri.init_information (x,X);
	i.update(); iri.update();
	std::cout << i.x << i.X <<std::endl;;
	std::cout << i.y << i.Y <<std::endl;;
	std::cout << iri.x << iri.X <<std::endl;;
	std::cout << iri.y << iri.Y <<std::endl;;
}


void test_unique()
{
	SymMatrix X(2,2);
	X(0,0) = 2.;
	X(1,0) = X(0,1) = 1.1;
	X(1,1) = 3.;
	Vec x(2);
	x[0] = 5.;
	x[1] = 7.;

	Test_random<Float> r;
	SIR_kalman_filter sf(2, 100, r);
	sf.init_kalman (x,X);

	std::cout << sf.unique_samples() <<',' << sf.stochastic_samples << std::endl;

	General_LiUnAd_observe_model obs(2,1);
	obs.Hx(0,0) = 1;
	obs.Hx(0,1) = 1;
	obs.Zv[0] = 1;
	Vec z(1);
	z[0] = 12.5;

	sf.observe (obs,z);
	sf.update();
	std::cout << sf.unique_samples() <<',' << sf.stochastic_samples << std::endl;

	sf.rougheningK = 0;
	sf.observe (obs,z);
	sf.update();
	std::cout << sf.unique_samples() <<',' << sf.stochastic_samples << std::endl;
}


void other_tests()
{
	test_SPD();
}


/* Boost Random
namespace {
template boost::uniform_smallint<boost::mt19937, int>;
template boost::uniform_int<boost::mt19937, int>;
typedef float Float;
template boost::uniform_01<boost::mt19937, Float>;
template boost::uniform_real<boost::mt19937, Float>;
template boost::triangle_distribution<boost::mt19937, Float>;
template boost::bernoulli_distribution<boost::mt19937, Float>;
template boost::cauchy_distribution<boost::mt19937, Float>;
template boost::exponential_distribution<boost::mt19937, Float>;
template boost::geometric_distribution<boost::mt19937, Float>;
template boost::normal_distribution<boost::mt19937, Float>;
template boost::lognormal_distribution<boost::mt19937, Float>;
template boost::uniform_on_sphere<boost::mt19937, Float>;
}
*/
