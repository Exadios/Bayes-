/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 */
#define NO_TESTS

#ifndef NO_TESTS
#include "BayesFilter/allFilters.hpp"
#include "BayesFilter/matSup.hpp"
#include "Test/random.hpp"
#include <cmath>
#include <iostream>
#include <exception>
#include <boost/numeric/ublas/io.hpp>
#include <boost/limits.hpp>
#include <boost/bind.hpp>


using namespace Bayesian_filter;
using namespace Bayesian_filter_matrix;


// Instantiate complete schemes to check the templates
#include "BayesFilter/filters/average1.hpp"
template Average1_filter<Covariance_scheme>;

#include "BayesFilter/filters/indirect.hpp"
template Indirect_state_filter<Covariance_scheme>;
template Indirect_kalman_filter<Covariance_scheme>;

void nested_prod ()
{
	Matrix R(1,1), A(1,1),B(1,1),C(1,1);
	R = prod(A, prod(B,C));
	R = prod(A, B);
}

// Square
template <class scalar>
inline scalar sqr(scalar x)
{
	return x*x;
}


class Test_random : public Bayesian_filter_test::Boost_random, public SIR_random
/*
 * Random numbers for filters from Boost
 */
{
public:
	Float normal (const Float mean, const Float sigma)
	{
		return Boost_random::normal (mean, sigma);
	}
	void normal (DenseVec& v)
	{
		Boost_random::normal (v);
	}
	void uniform_01 (DenseVec& v)
	{
		Boost_random::uniform_01 (v);
	}
	void seed ()
	{
		Boost_random::seed();
	}
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

void test_small_inverse ()
{
	SymMatrix X(1,1);
	SymMatrix XI(X.size1(), X.size2());

	Float a = 0;
	Float rcond;
	X(0,0) = a;
	std::cout << X << std::endl;
	rcond = UdUinversePD (XI, X);
	std::cout << rcond << XI << std::endl;

	X(0,0) = 1./a;
	std::cout << X << std::endl;
	rcond = UdUinversePD (XI, X);
	std::cout << rcond << XI << std::endl;
}

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
	noalias(P) = prod_SPD(R,S, RStemp);
	std::cout << P << std::endl;

	// Try prod_SPDT
	ColMatrix RT(1,2);
	RowMatrix SRtemp(S.size1(), RT.size2());
	noalias(P) = prod_SPDT (RT, S, SRtemp);
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
	Information_root_info_scheme iri(2);
	Information_scheme i(2);
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

	Test_random r;
	SIR_kalman_scheme sf(2, 100, r);
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

void numeric_tested()
{
	float a = 0;
	float b = 1.f/a;
	float c = 1.f/b;
	float d = -1.f/b;
	std::cout << a <<','<< b <<','<<c <<','<< d <<std::endl;
	std::cout << a*1.f <<','<< b*1.f <<','<<c*1.f <<','<< d*1.f <<std::endl;
	std::cout << (c==0) <<','<< (c!=0) << std::endl;
	std::cout << (c<0) <<','<< (c>0) << std::endl;
	std::cout << (c>=0) <<','<< (c<=0) << std::endl;

	std::cout << (d==0) <<','<< (d!=0) << std::endl;
	std::cout << (d<0) <<','<< (d>0) << std::endl;
	std::cout << (d>=0) <<','<< (d<=0) << std::endl;
}

void test_sym_proxy()
/*
 * Problems accessing the symmetric half via a proxy
 * IF packed proxy assign is used then only the represented half is copied
 * in the assign. The internal equality test will pickup the error
 */
{
	ublas::symmetric_matrix<double,ublas::upper> S(2,2);
	ublas::vector<double> v(2);
	S(0,0) = 1; S(1,0) = 2; S(1,1) = 3;
	v(0) = 5; v(1) = 6;
	noalias(ublas::row(S, 0)) = v;		// OK top row is all in upper
	noalias(ublas::row(S, 1)) = v;		// ERROR only one element in upper for packed proxy assign
}

void utinverse_test()
{
	UTriMatrix m(3,3);
	m(0,0) = 2; m(0,1) = 4; m(0,2) = 5;
	m(1,1) = 3; m(1,2) = 7;
	m(2,2) = 5;
	std::cout << m << std::endl;
	std::cout << UTinverse(m) << std::endl;
	std::cout << m << std::endl;
	std::cout << UTinverse(m) << std::endl;
	std::cout << m << std::endl;
}

void other_tests()
{
	try {
		// Tests go here
		test_SPD ();
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
}


#else
void other_tests()
{
}
#endif
