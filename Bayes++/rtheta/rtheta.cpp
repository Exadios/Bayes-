/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

 /*
 * Implement a NON-LINEAR range angle observer
 */

#include "BayesFilter/allFilters.hpp"
#include <angle.hpp>			// Angle arithmatic header from Ms
#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>
#include <boost/format.hpp>
#include <boost/limits.hpp>


using namespace Bayesian_filter;
using namespace Bayesian_filter_matrix;
using namespace angleArith;


const size_t NX = 2+1;			// State dimension (x,y) and some empty dummies so GqG' is singular
const size_t NQ = 2;			// State dimension of noise
const size_t NZ = 2+2;			// Observation dimension (r,a) and some empty dummies
const size_t NS = 1000;			// Number of samples for a sampled representation

const bool RA_MODEL = true;		// Use Range angle NON-linear model (requires normalising angle)
const bool NOISE_MODEL = true;		// Add noise to truth model
const bool TRUTH_STATIONARY = false;// Truth model setup
const Float INIT_XY[2] = {1.,Float(-0.2)};// XY initial position
const Float TARGET[2] = {-11.,0.};	// XY position of target

const Float RANGE_NOISE = NOISE_MODEL ? Float(0.1) : Float(1e-6);
const Float ANGLE_NOISE = NOISE_MODEL ? Float(5. * angle<Float>::Pi / 180.) : Float(1e-6);
const Float Z_CORRELATION = Float(0e-1);	// Correlated observation model
const Float X_NOISE = Float(0.05);		// Prediction model
const Float Y_NOISE = Float(0.09);
const Float XY_NOISE_COUPLING = Float(0.05);
const Float INIT_X_NOISE = Float(0.07);
const Float INIT_Y_NOISE = Float(0.10);
const Float INIT_XY_NOISE_CORRELATION = Float(0.04);


// Square 
template <class scalar>
inline scalar sqr(scalar x)
{
	return x*x;
}

// Random numbers for filters from Boost
class Boost_random : public SIR_random
/*
 * Random number distributions
 */
{
public:
	Boost_random() : gen_normal(rng), gen_uniform(rng)
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
	void reseed()
	{
		rng.seed();
	}
private:
	boost::mt19937 rng;
	boost::normal_distribution<boost::mt19937, Float> gen_normal;
	boost::uniform_01<boost::mt19937, Float> gen_uniform;
} Random, Random2;


/*
 * Linear prediction model
 *  x static and y tending to x
 *  Correlated additive noise
 */
class pred_model : public General_LiInAd_predict_model
{
public:
	pred_model();
private:
	class Predict_random : public General_LiInAd_predict_model::Random
	/*
	 * Uses global Random
	 */
	{
	public:
		void normal(FM::DenseVec& v)
		{
			::Random2.normal(v);
		}
	} localRandom;

	Vec x_pred;
};


pred_model::pred_model () :
	General_LiInAd_predict_model(NX,NQ, localRandom), x_pred(NX)
// Construct constant model
{
	// Build Fx, Identity all except active partition
	FM::identity (Fx);
	Fx(0,0) = 1.;
	Fx(0,1) = 0.;
	Fx(1,0) = Float(0.1);
	Fx(1,1) = Float(0.9);

	// Build q,G, Test addition parts using coupled unity noise
	q = ublas::scalar_vector<Float>(NQ, 1);
	q[0] = sqr(X_NOISE);
	q[1] = sqr(Y_NOISE);
	G = ublas::scalar_matrix<Float>(NX,NQ, 1);
	G(0,0) = 1;
	G(0,1) = XY_NOISE_COUPLING;
	G(1,0) = XY_NOISE_COUPLING;
	G(1,1) = 1;

	// Build inverse Fx, Identity all except active partition
	Information_root_filter::inverse_Fx (inv.Fx, Fx);

	// Build Inverse q
	for (size_t qi=0; qi < q.size(); ++qi)
		inv.q[qi] = Float(1) / q[qi];
}


/*
 * Observation model: gradient linearised versions
 *  Linearisation state must be set before use
 *  Both uncorrelated and correlated representation
 *  Correlated represtation is just the uncorrelated adapted so the
 *  correlation term are added.
 */
class uobs_model : public General_LzUnAd_observe_model
{
public:
	uobs_model();
	const Vec& h(const Vec& x) const;
	void normalise (Vec& z_denorm, const Vec& z_from) const;
	void state (const Vec& x);
private:
	mutable Vec z_pred;
};

class cobs_model : public General_LzCoAd_observe_model
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
	General_LzUnAd_observe_model(NX,NZ),
	z_pred(NZ)
{
	// Observation covariance uncorrelated
	Zv = ublas::scalar_vector<Float>(Zv.size(), Float(1));
	Zv[0] = sqr(RANGE_NOISE);
	Zv[1] = (RA_MODEL ? sqr(ANGLE_NOISE): sqr(RANGE_NOISE));
}

cobs_model::cobs_model(uobs_model& u) :
	General_LzCoAd_observe_model(NX,NZ),
	uobs(u)
{
	FM::identity (Z);
	// Create the correlation in Z
	Z(0,0) = uobs.Zv[0]; Z(1,1) = uobs.Zv[1];
	Z(1,0) = Z(0,1) = sqrt(Z(0,0))*sqrt(Z(1,1))*Z_CORRELATION;
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
	if (RA_MODEL)
		z_denorm[1] = angle<Float>(z_denorm[1]).from (z_from[1]);
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

walk::walk (const Vec start, const bool fixed) : x(NX), x_pred(NX), rootq(NQ), m_base(NX)
{
	m_fixed = fixed;
	m_base = start;
	x = m_base;
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
		x_pred = f(x);
		x = x_pred;
		x += nc;
	}
}


// Special for Information filter using Linrz predict
class Information_linrz_filter : public Information_filter
{
public:
	Information_linrz_filter (size_t x_size, size_t z_initialsize = 0) :
		Information_filter (x_size, z_initialsize)
	{}
	Float predict (Linrz_predict_model& f)
	// Enforce use of Linrz predict
	{
		return Information_filter::predict (f);
	}
};

/*
 * Filter under test. Initialised for state and covariance
 */
template <class TestScheme>
class Filter : public Filter_scheme<TestScheme>
{
public:
	Filter (const Vec& x_init, const SymMatrix& X_init);
};

template <class TestScheme>
Filter<TestScheme>::Filter (const Vec& x_init, const SymMatrix& X_init)
	: Filter_scheme<TestScheme> (x_init.size(), NQ, NZ)
{
	init_kalman (x_init, X_init);
}

// Specialise for SIR_kalman
template <>
class Filter<SIR_kalman_filter> : public SIR_kalman_filter
{
public:
	Filter (const Vec& x_init, const SymMatrix& X_init);
	Float update_resample ()
	// Modifiy Default SIR_filter update
	{
		return SIR_filter::update_resample (Systematic_resampler());
	}
};

Filter<SIR_kalman_filter>::Filter (const Vec& x_init, const SymMatrix& X_init)
	: SIR_kalman_filter (x_init.size(), NS, ::Random2)
{
	init_kalman (x_init, X_init);
}

// Specialise for Iterated_covariance_filter
template <>
class Filter<Iterated_covariance_filter> : public Iterated_covariance_filter
{
public:
	Filter (const Vec& x_init, const SymMatrix& X_init);
	bool observe_iteration_end ();
private:
	unsigned limit, li;
};

Filter<Iterated_covariance_filter>::Filter (const Vec& x_init, const SymMatrix& X_init)
	: Iterated_covariance_filter (x_init.size())
{
	li = limit = 10;
	init_kalman (x_init, X_init);
}

bool Filter<Iterated_covariance_filter>::observe_iteration_end ()
{
	--li;
	if (li == 0)
	{
		li = limit;
		return true;
	}
	else
		return false;
}



/*
 * Compare Two filters
 */

template<class Tf1, class Tf2>
class CCompare
{
public:
	CCompare (const Vec x_init, const SymMatrix X_init, unsigned nIterations);
private:
	void doIt (unsigned nIterations);
	// Attributes
	Tf1 f1;
	Tf2 f2;
	walk truth;
	pred_model f;
	uobs_model uh;
	cobs_model ch;

	// Implementation
	void dumpCompare ();

	Vec ztrue, z;

	DenseVec f1_xpred, f2_xpred;
	DenseSymMatrix f1_Xpred;
	DenseSymMatrix f2_Xpred;
};

template<class Tf1, class Tf2>
CCompare<Tf1,Tf2>::CCompare (const Vec x_init, const SymMatrix X_init, unsigned nIterations) :
	f1(x_init, X_init),
	f2(x_init, X_init),
	truth (x_init, TRUTH_STATIONARY),
	uh(),
	ch(uh),
	ztrue(NZ), z(NZ),
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

	doIt (nIterations);
}

template<class Tf1, class Tf2>
void CCompare<Tf1, Tf2>::dumpCompare ()
{
	Float zx, zy;
	if (RA_MODEL) {
		zx = truth.x[0] + z[0] * cos (z[1]);
		zy = truth.x[1] + z[0] * sin (z[1]);
	}
	else {
		zx = truth.x[0] + z[0];
		zy = truth.x[1] + z[1];
	}

	// Comparison
	{
		using std::cout; using std::endl;
		using boost::format;

		// Compomparision and truth line
		//	x(0)diff, x(1)diff, truth.x(0), truth.x(1), zx, zy)
		cout << format("*%11.4g %11.4g * %10.3f %10.3f  %10.3f %10.3f")
			 	% (f1.x[0]-f2.x[0]) % (f1.x[1]-f2.x[1])
			 	% truth.x[0] % truth.x[1] % zx % zy << endl;

		format state(" %11.4g %11.4g * %10.3f %10.3f");
		format covariance(" %12.4e %12.4e %12.4e");

		// Filter f1 performace
		//		x[0]err, x[1]err, x[0], x[1],  Xpred*3, X*3
		cout << state % (f1.x[0]-truth.x[0]) % (f1.x[1]-truth.x[1])
				% f1.x[0] % f1.x[1] << endl;
		cout << covariance % f1_Xpred(0,0) % f1_Xpred(1,1) % f1_Xpred(1,0);
		cout << covariance % f1.X(0,0) % f1.X(1,1) % f1.X(1,0) << endl;
		// Filter f2 performace
		//		x[0]err, x[1]err, x[0], x[1], x[01]dist,  Xpred*3, X*3
		cout << state % (f2.x[0]-truth.x[0]) % (f2.x[1]-truth.x[1])
				% f2.x[0] % f2.x[1] << endl;
		cout << covariance % f2_Xpred(0,0) % f2_Xpred(1,1) % f2_Xpred(1,0);
		cout << covariance % f2.X(0,0) % f2.X(1,1) % f2.X(1,0) << endl;
	}
}

template<class Tf1, class Tf2>
void CCompare<Tf1, Tf2>::doIt (unsigned nIterations)
{
	// Update the filter x,X representation
	f1.update (); f2.update ();
	dumpCompare();

	for (unsigned i = 0; i < nIterations; i++ ) {
		truth.predict ();		// Predict truth model
		f1.predict (f);			// Predict filters
		f2.predict (f);

		// Update the filter
		f1.update (); f2.update ();

		f1_xpred = f1.x; f2_xpred = f2.x;
		f1_Xpred = f1.X; f2_Xpred = f2.X;

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
			uh.state(f1.x);	f1.observe (uh, z);
			uh.state(f2.x); f2.observe (uh, z);
		}
		else
		{
			ch.state(f1.x);	f1.observe (ch, z);
			ch.state(f2.x);	f2.observe (ch, z);
		}

		// Update the filter
		f1.update (); f2.update ();
		dumpCompare();
		// DEBUG char c;std::cin>>c;
	}
}


int main()
{
	// Other things I might want to test
	extern void other_tests();
	other_tests();

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

	// Initialise and do the comparison
	std::cout << "udfilter, ufilter " << "RA_MODEL:" << RA_MODEL << " NOISE_MODEL:" << NOISE_MODEL << " TRUTH_STATIONARY:" << TRUTH_STATIONARY << std::endl;
	Random.reseed();
	CCompare<Filter<UD_filter>, Filter<Unscented_filter> > test1(x_init, X_init, 4);
	std::cout << std::endl;

	std::cout << "cfilter, ifilter " << "RA_MODEL:" << RA_MODEL << " NOISE_MODEL:" << NOISE_MODEL << " TRUTH_STATIONARY:" << TRUTH_STATIONARY << std::endl;
	Random.reseed();
	CCompare<Filter<Covariance_filter>, Filter<Information_linrz_filter> > test2(x_init, X_init, 4);
	std::cout << std::endl;

	std::cout << "irfilter, ilfilter " << "RA_MODEL:" << RA_MODEL << " NOISE_MODEL:" << NOISE_MODEL << " TRUTH_STATIONARY:" << TRUTH_STATIONARY << std::endl;
	Random.reseed();
	CCompare<Filter<Information_root_filter>, Filter<Information_filter> > test3(x_init, X_init, 4);
	std::cout << std::endl;

	std::cout << "sfilter, itfilter " << "RA_MODEL:" << RA_MODEL << " NOISE_MODEL:" << NOISE_MODEL << " TRUTH_STATIONARY:" << TRUTH_STATIONARY << std::endl;
	Random.reseed();
//	CCompare<Filter<SIR_kalman_filter>, Filter<Iterated_covariance_filter> > test4(x_init, X_init, 4);
//	CCompare<Filter<Covariance_filter>, Filter<Iterated_covariance_filter> > test4(x_init, X_init, 4);
	std::cout << std::endl;

	return 0;
}
