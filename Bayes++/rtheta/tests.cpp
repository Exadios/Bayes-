/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 */
#include <ctime>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <boost/limits.hpp>
#include <boost/random.hpp>
#include "timing.h"
#include "angle.h"

#include "BayesFilter/allFlt.h"


using namespace Bayesian_filter;
using namespace Bayesian_filter_matrix;
using namespace angleArith;


// Instantiate Average1_filter to check the template
#include "BayesFilter/filters/average1.h"
template average1_filter<Covariance_filter>;

// Square
template <class scalar>
inline scalar sqr(scalar x)
{
	return x*x;
}


class Test_random
/*
 * Random number distributions
 */
{
public:
	Test_random() : gen_normal(rng), gen_uniform(rng), gen_exponential(rng, 1.)
	{}
	double normal(const double mean, const double sigma)
	{
		boost::normal_distribution<boost::mt19937> gen(rng, mean, sigma);
		return gen();
	}
	void normal(Vec& v)
	{
		std::generate (v.begin(), v.end(), gen_normal);
	}
	void uniform_01(Vec& v)
	{
		std::generate (v.begin(), v.end(), gen_uniform);
	}
	void exponential_1(Vec& v)
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
	boost::normal_distribution<good_random> gen_normal;
	boost::uniform_01<good_random> gen_uniform;
	boost::exponential_distribution<good_random> gen_exponential;
};

void test_inverse()
{	// DEBUG Just for Debugging!

	Matrix U (3,3);
	U.clear();
	U(0,0) = sqr(1.);
	U(1,1) = sqr(1.);
	U(2,2) = sqr(1.);
	U(0,1) = 2.;
	U(0,2) = 3.;
	U(1,2) = 5.;
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

	Matrix X(U.size1(),U.size2());
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

	// Try mult_SPD
	P.clear();
	mult_SPD(R,S,P,Stemp);
	std::cout << P << std::endl;

	// Try prod_SPD
	RowMatrix RStemp(R.size1(), S.size2());
	P.assign( prod_SPD(R,S, RStemp) );
	std::cout << P << std::endl;

	// Try mult_SPDT
	RowMatrix RT(1,2);
	RT(0,0) = 1.;
	RT(0,1) = 2.;
	std::cout << RT << std::endl;

	P.clear();
	mult_SPDT(RT,S,P,Stemp);
	std::cout << P << std::endl;
}

void test_SPD_all()
{
	SymMatrix S(1,1);
	Vec x(1);
	S(0,0) = 5.;
	x[0] = 1.;

	double v1 = ublas::inner_prod (x, ublas::prod(S,x) );

	double v2 = prod_SPDT(x,S);

	test_SPD();
}

void test_temp_prod()
{
	Matrix A(2,2);
	Matrix R = prod( ublas::prod<Matrix>(A,A), ublas::prod<Matrix>(A,A));
}


#ifdef REMOVED
// Test MTL storage problems
void test_storage()
{
	UTriMatrix A(2,2);
	A(0,0) = 1.;
	A(1,0) = 9.;
	A(1,1) = 2.;
	std::cout << A << std::endl;

	// Show effect of banded storage
	std::cout << A[0] << std::endl;
	std::cout << A[1] << std::endl;

	SymMatrix S(2,2);
	set (S, 0.);

	// Try mult_SPD
	mult_SPDi(A, S);
	std::cout << S << std::endl;

	// Try mult
	set (S, 0.);
	mult (A, A, S);
	std::cout << S << std::endl;

	// Try mult
	set (S, 0.);
	mult (A, trans(A), S);
	std::cout << S << std::endl;
}
#endif

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
	double rcond = UCfactor (UD, S);
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

#ifdef REMOVED
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

	i.init(x,X);
	iri.init(x,X);
	i.update(); iri.update();
	FM::print_vector(i.x); FM::print_all_matrix(i.X);
	FM::print_vector(i.y); FM::print_all_matrix(i.Y);
	FM::print_vector(iri.x); FM::print_all_matrix(iri.X);
	FM::print_vector(iri.y); FM::print_all_matrix(iri.Y);

	i.init_information(x,X);
	iri.init_information(x,X);
	i.update(); iri.update();
	FM::print_vector(i.x); FM::print_all_matrix(i.X);
	FM::print_vector(i.y); FM::print_all_matrix(i.Y);
	FM::print_vector(iri.x); FM::print_all_matrix(iri.X);
	FM::print_vector(iri.y); FM::print_all_matrix(iri.Y);
}
#endif

/*
void test_addative()
{
	Matrix G(2,1);
	Vec q(1);
	class ff : public Function_model
	{
	virtual const FM::Vec& fx(const FM::Vec& x) const {return x;7}
	} f;
	Simple_addative_predict_model s(f,G,q);
}
*/

void other_tests()
{
	//	test_();
	test_SPD_all();
	test_temp_prod();
}
