/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Implement a NON-LINEAR range angle observer for many filter schemes
 * The model parameters allow model sizes and the non-linearity.
 * This provides an excellent vehicle to test each scheme. Use for regression testing
 */

#include "BayesFilter/allFilters.hpp"
#include "BayesFilter/schemeFlt.hpp"
#include "Test/random.hpp"
#include <angle.hpp>			// Angle arithmatic header from Ms
#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Bayesian_filter;
using namespace Bayesian_filter_matrix;
using namespace angleArith;


const std::size_t NX = 2+1;			// State dimension (x,y) and some more to show noise coupling or singular X
const std::size_t NQ = 2;			// State dimension of noise
const std::size_t NZ = 2+2;			// Observation dimension (r,a) and some empty dummies
const std::size_t NS = 1000;			// Number of samples for a sampled representation

const bool RA_MODEL = true;		// Use Range angle NON-linear model (requires normalising angle)
const bool NOISE_MODEL = true;		// Add noise to truth model
const bool TRUTH_STATIONARY = false;// Truth model setup
const Float INIT_XY[2] = {1.,-0.2};	// XY initial position
const Float TARGET[2] = {-11.,0.};	// XY position of target

const Float RANGE_NOISE = NOISE_MODEL ? Float(0.1) : Float(1e-6);
const Float ANGLE_NOISE = NOISE_MODEL ? Float(5. * angle<Float>::Deg2Rad) : Float(1e-6);
const Float Z_CORRELATION = Float(0e-1);	// (Un)Correlated observation model

const Float X_NOISE = Float(0.05);		// predict model
const Float Y_NOISE = Float(0.09);
const Float XY_NOISE_COUPLING = Float(0.05);
const Float Q_NOISE = Float(1.0);		// Noise in addition Q terms
const Float G_COUPLING = Float(1.0);	// Coupling in addition G terms

const Float INIT_X_NOISE = Float(0.07);
const Float INIT_Y_NOISE = Float(0.10);
const Float INIT_XY_NOISE_CORRELATION = Float(0.4);
const Float INIT_2_NOISE = Float(0.09);	// Use zero for singular X	
const Float INIT_2_NOISE_CORRELATION = Float(0.5);


// Square 
template <class scalar>
inline scalar sqr(scalar x)
{
	return x*x;
}

class Rtheta_random : public Bayesian_filter_test::Boost_random, public SIR_random
/*
 * Random numbers for filters from Boost
 */
{
public:
	Float normal (const Float mean, const Float sigma)
	{
		Float f = Boost_random::normal (mean, sigma);
		return f;
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
} Random, Random2;


/*
 * Linear predict model
 *  x static and y tending to x
 *  Correlated additive noise
 */
class pred_model : public Sampled_LiInAd_predict_model
{
public:
	pred_model();
};


pred_model::pred_model () :
	Sampled_LiInAd_predict_model(NX,NQ, Random2)
// Construct constant model
{
	// Build Fx, Identity all except active partition
	FM::identity (Fx);
	Fx(0,0) = 1.;
	Fx(0,1) = 0.;
	Fx(1,0) = Float(0.1);
	Fx(1,1) = Float(0.9);

	// Build q,G, Test addition parts using coupled noise
	q = ublas::scalar_vector<Float>(NQ, Q_NOISE);
	q[0] = sqr(X_NOISE);
	q[1] = sqr(Y_NOISE);
	G = ublas::scalar_matrix<Float>(NX,NQ, G_COUPLING);
	G(0,0) = 1;
	G(0,1) = XY_NOISE_COUPLING;
	G(1,0) = XY_NOISE_COUPLING;
	G(1,1) = 1;

	// Build inverse Fx, Identity all except active partition
	DenseColMatrix denseinvFx (Fx.size1(), Fx.size2());
	Information_root_scheme::inverse_Fx (denseinvFx, Fx);
	inv.Fx = denseinvFx;
}


/*
 * Observation model: gradient linearised versions
 *  Linearisation state must be set before use
 *  Both uncorrelated and correlated representation
 *  Correlated represtation is just the uncorrelated adapted so the
 *  correlation term are added.
 */
#ifndef MODEL_LINEAR
	typedef General_LzUnAd_observe_model Observe_un_special;
	typedef General_LzCoAd_observe_model Observe_co_special;
#else
	typedef General_LiUnAd_observe_model Observe_un_special;
	typedef General_LiCoAd_observe_model Observe_co_special;
#endif

class uobs_model : public Observe_un_special
{
public:
	uobs_model();
	const Vec& h(const Vec& x) const;
	void normalise (Vec& z_denorm, const Vec& z_from) const;
	void state (const Vec& x);
private:
	mutable Vec z_pred;
};

class cobs_model : public Observe_co_special
{
public:
	cobs_model(uobs_model& u);
	const Vec& h(const Vec& x) const;
	void normalise (Vec& z_denorm, const Vec& z_from) const;
	void state (const Vec& x);
private:
	uobs_model& uobs;
};

uobs_model::uobs_model() :
	Observe_un_special(NX,NZ),
	z_pred(NZ)
{
	// Observation covariance uncorrelated
	Zv = ublas::scalar_vector<Float>(Zv.size(), Float(1));
	Zv[0] = sqr(RANGE_NOISE);
	Zv[1] = (RA_MODEL ? sqr(ANGLE_NOISE): sqr(RANGE_NOISE));
}

cobs_model::cobs_model(uobs_model& u) :
	Observe_co_special(NX,NZ),
	uobs(u)
{
	FM::identity (Z);
	// Create the correlation in Z
	Z(0,0) = Float(uobs.Zv[0]); Z(1,1) = Float(uobs.Zv[1]);	// ISSUE mixed type proxy assignment
	Z(1,0) = Float(Z(0,1) = sqrt(Z(0,0))*sqrt(Z(1,1))*Z_CORRELATION);
}

void uobs_model::state (const Vec& x)
{
	Float dx = TARGET[0] - x[0];
	Float dy = TARGET[1] - x[1];

	Hx.clear();
	if (RA_MODEL) {
		Float distSq = dx*dx + dy*dy;
		Float dist = sqrt (distSq);
		Hx(0,0) = -dx / dist;
		Hx(0,1) = -dy / dist;
		Hx(1,0) = +dy / distSq;
		Hx(1,1) = -dx / distSq;
	}
	else {
		Hx(0,0) = -1.;
		Hx(0,1) = 0.;
		Hx(1,0) = 0.;
		Hx(1,1) = -1.;
	}
}

void cobs_model::state (const FM::Vec& x)
{
	uobs.state (x);
	Hx = uobs.Hx;
}

const Vec& uobs_model::h (const Vec& x) const
{
	Float dx = TARGET[0] - x[0];
	Float dy = TARGET[1] - x[1];

	z_pred.clear();
	if (RA_MODEL) {
		Float distSq = dx*dx + dy*dy;
		Float dist = sqrt (distSq);

		z_pred[0] = dist;
		z_pred[1] = std::atan2 (dy, dx);
	}
	else {
		z_pred[0] = dx;
		z_pred[1] = dy;
	}
	return z_pred;
}

const Vec& cobs_model::h (const Vec& x) const
{
	return uobs.h(x);
}

void uobs_model::normalise (Vec& z_denorm, const Vec& z_from) const
{
	if (RA_MODEL) {
		z_denorm[1] = angle<Float>(z_denorm[1]).from (z_from[1]);
	}
}

void cobs_model::normalise (Vec& z_denorm, const Vec& z_from) const
{
	uobs.normalise(z_denorm, z_from);
}

/*
 * A dynamic system with noise, or a fixed additive noise
 */
class walk : private pred_model
{
public:
	walk (const Vec start, const bool fixed = false);
	void predict ();
	Vec x, x_pred, rootq;
private:
	bool m_fixed;
	Vec m_base;
};

walk::walk (const Vec start, const bool fixed) : x(start), x_pred(NX), rootq(NQ), m_fixed(fixed), m_base(start)
{
	for (Vec::const_iterator qi = q.begin(); qi != q.end(); ++qi) {
		rootq[qi.index()] = std::sqrt(*qi);
	}
}

void walk::predict ()
{
						// Correlated additive random noise
	DenseVec n(rootq.size()), nc(x.size());
	::Random.normal (n);		// independant zero mean normal
								// multiply elements by std dev
	for (DenseVec::iterator ni = n.begin(); ni != n.end(); ++ni) {
		*ni *= rootq[ni.index()];
	}
	nc = prod(G,n);				// correlate

	if (m_fixed) {				// Randomize based on assumed noise
		x = m_base + nc;
	}
	else {
		noalias(x_pred) = f(x);
		x = x_pred;
		noalias(x) += nc;
	}
}


// Special for Information scheme using Linrz predict
class Information_linrz_scheme : public Information_scheme
{
public:
	Information_linrz_scheme (std::size_t x_size) :
		Kalman_state_filter (x_size),
		Information_state_filter (x_size),
		Information_scheme (x_size)
	{}
	Float predict (Linrz_predict_model& f)
	// Enforce use of Linrz predict
	{
		return Information_scheme::predict (f);
	}
};

// Filter_scheme Information_linrz_scheme specialisation
template <>
Filter_scheme<Information_linrz_scheme>::Filter_scheme(std::size_t x_size, std::size_t q_maxsize) :
	Kalman_state_filter (x_size),
	Information_state_filter (x_size),
	Information_linrz_scheme (x_size)
{}

/*
 * Filter under test. Initialised for state and covariance
 */
template <class TestScheme>
class Filter
{
public:
	Filter (const Vec& x_init, const SymMatrix& X_init);
	Filter_scheme<TestScheme> ts;
	template <class P>
	void predict (P& pmodel)
	{
		ts.predict (pmodel);
	}
	template <class O>
	void observe (O& omodel, const Vec& z)
	{
		ts.observe (omodel, z);
	}
	void update ()
	{
		ts.update ();
	}
	const Vec& x()
	{	return ts.x;
	}
	const SymMatrix& X()
	{	return ts.X;
	}
	void dump_state()
	{	// output any additional state variables
	}
};

template <class TestScheme>
Filter<TestScheme>::Filter (const Vec& x_init, const SymMatrix& X_init) :
	ts (x_init.size(), NQ)
{
	ts.init_kalman (x_init, X_init);
}

// Specialise for SIR_kalman
template <>
Filter<SIR_kalman_scheme>::Filter (const Vec& x_init, const SymMatrix& X_init) :
	ts (x_init.size(), NS, ::Random2)
{
	ts.init_kalman (x_init, X_init);
}

// dump_state specialisations
template <>
void Filter<Information_root_scheme>::dump_state()
{	// output any additional state variables
	std::cout << ts.r << ts.R << std::endl;
}

template <>
void Filter<Unscented_scheme>::dump_state()
{	// output any additional state variables
	std::cout << ts.XX << std::endl;
}

/*
 * Compare Two filters
 */

class CCompare
{
public:
	CCompare (const Vec x, const SymMatrix X, unsigned compare_iterations);
	template<class Tf1, class Tf2>
	void compare ();
private:
	// Test configuration
	const Vec x_init; const SymMatrix X_init;
	unsigned nIterations;
	pred_model f;
	uobs_model uh;
	cobs_model ch;

	void dump_state ();
	bool tolerably_equal ();

	// Test state
	Vec ztrue, z;
	Vec t_x;
	Vec f1_x, f2_x;
	SymMatrix f1_X, f2_X;
	Vec f1_xpred, f2_xpred;
	SymMatrix f1_Xpred, f2_Xpred;
};

CCompare::CCompare (const Vec x, const SymMatrix X, unsigned compare_iterations) :
	x_init(x), X_init(X),
	nIterations (compare_iterations),
	uh(),
	ch(uh),
	ztrue(NZ), z(NZ),
	t_x(NX),
	f1_x (NX),
	f2_x (NX),
	f1_X (NX, NX),
	f2_X (NX, NX),
	f1_xpred (NX),
	f2_xpred (NX),
	f1_Xpred (NX, NX),
	f2_Xpred (NX, NX)
{
	// Initialises test variables
	z.clear();

	f1_xpred.clear();
	f2_xpred.clear();
	f1_Xpred.clear();
	f2_Xpred.clear();
}

void CCompare::dump_state ()
{
	// Additional Scheme state
	//f1.dump_state(); f2.dump_state();

	Float zx, zy;
	if (RA_MODEL) {
		zx = t_x[0] + z[0] * cos (z[1]);
		zy = t_x[1] + z[0] * sin (z[1]);
	}
	else {
		zx = t_x[0] + z[0];
		zy = t_x[1] + z[1];
	}

	// Comparison
	{
		using std::cout; using std::endl;
		using boost::format;

		// Comparision and truth line
		//	x(0)diff, x(1)diff, t_x(0), t_x(1), zx, zy)
		cout << format("*%11.4g %11.4g * %10.3f %10.3f  %10.3f %10.3f")
			 	% (f1_x[0]-f2_x[0]) % (f1_x[1]-f2_x[1])
			 	% Float(t_x[0]) % Float(t_x[1]) % zx % zy << endl;	// ISSUE mixed type proxy assignment

		format state(" %11.4g %11.4g * %10.3f %10.3f");
		format covariance(" %12.4e %12.4e %12.4e");

		// Filter f1 performance
		//		x[0]err, x[1]err, x[0], x[1],  Xpred*3, X*3
		cout << state % (f1_x[0]-t_x[0]) % (f1_x[1]-t_x[1])
				% f1_x[0] % f1_x[1] << endl;
		cout << covariance % f1_Xpred(0,0) % f1_Xpred(1,1) % f1_Xpred(1,0);
		cout << covariance % f1_X(0,0) % f1_X(1,1) % f1_X(1,0) << endl;
		// Filter f2 performace
		//		x[0]err, x[1]err, x[0], x[1], x[01]dist,  Xpred*3, X*3
		cout << state % (f2_x[0]-t_x[0]) % (f2_x[1]-t_x[1])
				% f2_x[0] % f2_x[1] << endl;
		cout << covariance % f2_Xpred(0,0) % f2_Xpred(1,1) % f2_Xpred(1,0);
		cout << covariance % f2_X(0,0) % f2_X(1,1) % f2_X(1,0) << endl;
	}
}

bool CCompare::tolerably_equal ()
{
	bool tolerable = true;
#ifdef REMOVED
	// Tolerance bounds
	const Float truth_sigma = 100*5;	// Expect all estimates to lie withing this sigma bound
	const Float epsilon = 100*std::numeric_limits<Float>::epsilon();
	//using boost::test_toolsbox::close_at_tolerance;
	typedef boost::test_toolbox::close_at_tolerance<Float> close_at_tolerance;		// Floating point comparison

	// Filter estimate is close to truth
//	using boost::test_toolbox::check_is_closed
	Float t0 = f1_x[0]- t_x[0];
	Float t1 = f1_x[1]- t_x[1];
	Float T0 = f2_x[0]- t_x[0];
	Float T1 = f2_x[1]- t_x[1];
	tolerable &= close_at_tolerance(truth_sigma * std::sqrt(f1_X(0,0)), boost::test_toolbox::FPC_WEAK)  (f1_x[0], t_x[0]);
	tolerable &= close_at_tolerance(truth_sigma * std::sqrt(f1_X(1,1)), boost::test_toolbox::FPC_WEAK)  (f1_x[1], t_x[1]);
	tolerable &= close_at_tolerance(truth_sigma * std::sqrt(f2_X(0,0)), boost::test_toolbox::FPC_WEAK)  (f2_x[0], t_x[0]);
	tolerable &= close_at_tolerance(truth_sigma * std::sqrt(f2_X(1,1)), boost::test_toolbox::FPC_WEAK)  (f2_x[1], t_x[1]);

	// Filter states equal
	close_at_tolerance state_close(epsilon * 100);
	Float d0 = f1_x[0]- f2_x[0];
	Float d1 = f1_x[1]- f2_x[1];
	tolerable &= state_close (f1_x[0], f2_x[0]);
	tolerable &= state_close (f1_x[1], f2_x[1]);

	// Covraiance toleration
	close_at_tolerance var_close(epsilon * 50);
	close_at_tolerance cov_close(epsilon * 50);
	// Predict covariance equal
	Float x00 = f1_Xpred(0,0)- f2_Xpred(0,0);
	Float x11 = f1_Xpred(1,1)- f2_Xpred(1,1);
	Float x01 = f1_Xpred(1,0)- f2_Xpred(1,0);
	tolerable &= var_close (f1_Xpred(0,0), f2_Xpred(0,0));
	tolerable &= var_close (f1_Xpred(1,1), f2_Xpred(1,1));
	tolerable &= cov_close (f1_Xpred(1,0), f2_Xpred(1,0));
	// Observe covariance equal
	Float X00 = f1_X(0,0)- f2_X(0,0);
	Float X11 = f1_X(1,1)- f2_X(1,1);
	Float X01 = f1_X(1,0)- f2_X(1,0);
	tolerable &= var_close (f1_X(0,0), f2_X(0,0));
	tolerable &= var_close (f1_X(1,1), f2_X(1,1));
	tolerable &= cov_close (f1_X(1,0), f2_X(1,0));
#endif

	return tolerable;
}

template<class Tf1, class Tf2>
void CCompare::compare ()
{
	// Construct truth model
	walk truth (x_init, TRUTH_STATIONARY);
	
	// Construct filter to compare with initial true state
	Tf1 f1(x_init, X_init);
	Tf2 f2(x_init, X_init);

	// Update the filter x,X representation
	f1.update (); f2.update ();
	t_x = truth.x;
	f1_xpred = f1.x(); f2_xpred = f2.x();
	f1_Xpred = f1.X(); f2_Xpred = f2.X();
	f1_x = f1.x(); f2_x = f2.x();
	f1_X = f1.X(); f2_X = f2.X();
	z = ztrue = uh.h(truth.x);
	dump_state ();

	for (unsigned i = 0; i < nIterations; i++ ) {
		truth.predict ();		// Predict truth model
		f1.predict (f);			// Predict filters
		f2.predict (f);
		
		// Update the filter
		f1.update (); f2.update ();

		f1_xpred = f1.x(); f2_xpred = f2.x();
		f1_Xpred = f1.X(); f2_Xpred = f2.X();

		// Observation, true and Randomize based on an uncorrelated noise model
		ztrue = uh.h(truth.x);
		if (NOISE_MODEL) {
			z[0] = Random.normal (ztrue[0], sqrt(uh.Zv[0]) );
			z[1] = Random.normal (ztrue[1], sqrt(uh.Zv[1]) );
		}
		else
			z = ztrue;

		// Observe using model linearised about filter state estimate
		if (Z_CORRELATION == 0.)
		{
			uh.state(f1.x()); f1.observe (uh, z);
			uh.state(f2.x()); f2.observe (uh, z);
		}
		else
		{
			ch.state(f1.x()); f1.observe (ch, z);
			ch.state(f2.x()); f2.observe (ch, z);
		}

		// Update the filter
		f1.update (); f2.update ();
		t_x = truth.x;
		f1_x = f1.x(); f2_x = f2.x();
		f1_X = f1.X(); f2_X = f2.X();
		
		dump_state ();
		tolerably_equal ();
		// DEBUG char c;std::cin>>c;
	}
}


int main()
{
	// Other things I might want to test
	extern void other_tests();
	other_tests();

	// Use a know sequence for comparisons between systems
	Random.seed();
	Random2.seed();
	
	// Setup the test filters
	Vec x_init (NX);
	SymMatrix X_init (NX, NX);

	// Cartessian start position (in meters)
	x_init[0] = INIT_XY[0];
	x_init[1] = INIT_XY[1];
	// Initial state covariance, correlated
	FM::identity (X_init);
	X_init(0,0) = sqr(INIT_X_NOISE);
	X_init(1,1) = sqr(INIT_Y_NOISE);
	X_init(1,0) = X_init(0,1) = INIT_X_NOISE*INIT_Y_NOISE*INIT_XY_NOISE_CORRELATION;
	// Additional state correlation is useful for testing
	X_init(2,2) = sqr(INIT_2_NOISE);
	X_init(2,0) = X_init(0,2) = INIT_X_NOISE*INIT_2_NOISE*INIT_2_NOISE_CORRELATION;
	X_init(2,1) = X_init(1,2) = INIT_Y_NOISE*INIT_2_NOISE*INIT_2_NOISE_CORRELATION;

	// Test iterations
	CCompare test (x_init, X_init, 4);

	// Initialise and do the comparison
	std::cout << "udfilter, ufilter " << "RA_MODEL:" << RA_MODEL << " NOISE_MODEL:" << NOISE_MODEL << " TRUTH_STATIONARY:" << TRUTH_STATIONARY << std::endl;
	Random.seed();
	test.compare<Filter<UD_scheme>, Filter<Unscented_scheme> >();
	std::cout << std::endl;

	std::cout << "cfilter, ifilter " << "RA_MODEL:" << RA_MODEL << " NOISE_MODEL:" << NOISE_MODEL << " TRUTH_STATIONARY:" << TRUTH_STATIONARY << std::endl;
	Random.seed();
	test.compare<Filter<Covariance_scheme>, Filter<Information_linrz_scheme> >();
	std::cout << std::endl;

	std::cout << "irfilter, ilfilter " << "RA_MODEL:" << RA_MODEL << " NOISE_MODEL:" << NOISE_MODEL << " TRUTH_STATIONARY:" << TRUTH_STATIONARY << std::endl;
	Random.seed();
	test.compare<Filter<Information_root_scheme>, Filter<Information_scheme> >();
	std::cout << std::endl;

	std::cout << "sfilter, itfilter " << "RA_MODEL:" << RA_MODEL << " NOISE_MODEL:" << NOISE_MODEL << " TRUTH_STATIONARY:" << TRUTH_STATIONARY << std::endl;
	Random.seed();
	test.compare<Filter<SIR_kalman_scheme>, Filter<Iterated_covariance_scheme> >();
	std::cout << std::endl;

	return 0;
}
