/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * UdU' Factorisation of Covariance Filter.
 *
 * For efficiency UD_filter requires to know the maximum q_size of the predict_model
 * TODO:
 *  Fix efficiency issue in observe for nonlinear z rank 2 or more
 *  Extend Observe function to deal with PSD
 *  Make rcond return overall conditioning not minium of each sequential
 */
#include "UDFlt.hpp"
#include "matSup.hpp"
#include <boost/limits.hpp>

/* Filter namespace */
namespace Bayesian_filter
{


UD_filter::
UD_filter (size_t x_size, size_t q_maxsize, size_t z_initialsize) :
		Linrz_filter(x_size)
		, q_max(q_maxsize)
		, UD(x_size,2*x_size)
		, s(FM::Empty), Sd(FM::Empty)
		, d(x_size+q_max), dv(x_size+q_max), v(x_size+q_max)
		, a(x_size), b(x_size)
		, h1(x_size), w(x_size)
		, zp(FM::Empty)
		, UHx(FM::Empty)
		, zdecol(FM::Empty)
		, Gz(FM::Empty)

/*
 * Initialise filter and set the size of things we know about
 */
{
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

UD_filter&
 UD_filter::operator= (const UD_filter& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Linrz_filter::operator=(a);
	q_max = a.q_max;
	UD = a.UD;
	return *this;
}


void 
 UD_filter::init ()
/*
 * Initialise from a state and state coveriance
 * Computes UD factor from initial covaiance
 * Predcondition:
 *		X is PSD
 * Postcondition:
 *		UdU' is PSD
 */
{
					// Factorise X into left partition of UD
	FM::subcopy(X, UD);
	Float rcond = FM::UdUfactor (UD, UD.size1());
	rclimit.check_PSD(rcond, "Xi not PSD");
}


void 
 UD_filter::update ()
/*
 * Defactor UD back into X
 * Predcondition:
 *		UdU' is PSD (not checked)
 * Postcondition:
 *		X is PSD
 */
{
	size_t i,j,k;
	const size_t n = UD.size1();
						// Build X as lower tri (DU')'
	for (i = 0; i < n; ++i)				// 0..n-1
	{
		FM::SymMatrix::Row Xi(X,i);
		Float d = UD(i,i);
		Xi[i] = d;
		for (j = 0; j < i; ++j)			// 0..i-1
		{
			Xi[j] = UD(j,i)*d;
		}
	}
						// In place U *(DU') working down lower triangle of X */
	for (i = 0; i < n; ++i)				// 0..n-1
	{
		for (j = 0; j <= i; ++j)		// 0..i
		{
			FM::Matrix::Row UDi(UD,i);
			Float s = X(i,j);
			for (k = i+1; k < n; ++k)
				s += UDi[k] * X(k,j);
			X(j,i) = X(i,j) = s;
		}
	}
}


UD_filter::Float
 UD_filter::predict (Linrz_predict_model& f)
/*
 * Prediction using a diagonalised noise q, and its coupling G
 *  q can have order less then x and a matching G so GqG' has order of x
 * Predcondition:
 *		UdU' is PSD (not checked)
 * Postcondition:
 *		UdU' is PSD
 */
{
	x = f.f(x);			// Extended Kalman state predict is f(x) directly

						// Predict UD from model
	Float rcond = predictGq (f.Fx, f.G, f.q);
	rclimit.check_PSD(rcond, "X not PSD in predict");
	return rcond;
}


UD_filter::Float
 UD_filter::predictGq (const FM::Matrix& Fx, const FM::Matrix& G, const FM::Vec& q)
/*
 * MWG-S prediction from Bierman  p.132
 *  q can have order less then x and a matching G so GqG' has order of x
 * Predcondition:
 *		UdU' is PSD (not checked)
 * Postcondition:
 *		UdU' is PSD (see return value)
 *
 * Return:
 *		reciprocal condition number, -1. if negative, 0. if semi-definate (including zero)
 */
{
	size_t i,j,k;
	const size_t n = x.size();
	const size_t Nq = q.size();
	const size_t N = n+Nq;
	Float e;
					// Check preallocated space for q size
	if (Nq > q_max)
		filter_error("Predict model q larger than preallocated space");

	if (n > 0)		// Simplify reverse loop termination
	{
						// Augment d with q, UD with G
		for (i = 0; i < Nq; ++i)		// 0..Nq-1
		{
			d[i+n] = q[i];
		}
		for (j = 0; j < n; ++j)		// 0..n-1
		{
			FM::Matrix::Row UDj(UD,j);
			FM::Matrix::const_Row  Gj(G,j);
			for (i = 0; i < Nq; ++i)		// 0..Nq-1
				UDj[i+n] = Gj[i];
		}

						// U=Fx*U and diagonals retrived
		for (j = n-1; j > 0; --j)		// n-1..1
		{
						// Prepare d(0)..d(j) as temporary
			for (i = 0; i <= j; ++i)	// 0..j
				d[i] = UD(i,j);
			
						// Lower triangle of UD is implicity empty
			for (i = 0; i < n; ++i) 	// 0..n-1
			{
				FM::Matrix::Row UDi(UD,i);
				FM::Matrix::const_Row Fxi(Fx,i);
				UDi[j] = Fxi[j];
				for (k = 0; k < j; ++k)	// 0..j-1
					UDi[j] += Fxi[k] * d[k];
			}
		}
		d[0] = UD(0,0);

						//  Complete U = Fx*U
		for (j = 0; j < n; ++j)			// 0..n-1
		{
			UD(j,0) = Fx(j,0);
		}

						// The MWG-S algorithm on UD transpose
		j = n-1;
		do {							// n-1..0
			FM::Matrix::Row UDj(UD,j);
			e = 0.;
			for (k = 0; k < N; ++k)		// 0..N-1
			{
				v[k] = UDj[k];
				dv[k] = d[k] * v[k];
				e += v[k] * dv[k];
			}
			// Check diagonal element
			if (e > 0.)
			{
				// Positive definate
				UDj[j] = e;

				Float diaginv = 1. / e;
				for (k = 0; k < j; ++k)	// 0..j-1
				{
					FM::Matrix::Row UDk(UD,k);
					e = 0.;
					for (i = 0; i < N; ++i)	// 0..N-1
						e += UDk[i] * dv[i];
					e *= diaginv;
					UDj[k] = e;

					for (i = 0; i < N; ++i)	// 0..N-1
						UDk[i] -= e * v[i];
				}
			}//PD
			else if (e == 0.)
			{
				// Possibly Semidefinate, check not negative
				UDj[j] = e;

				// 1. / e is infinite
				for (k = 0; k < j; ++k)	// 0..j-1
				{
					FM::Matrix::Row UDk(UD,k);
					for (i = 0; i < N; ++i)	// 0..N-1
					{
						e = UDk[i] * dv[i];
						if (e != 0.)
							goto Negative;
					}
					// UD(j,k) uneffected
				}
			}//PD
			else
			{
				// Negative
				goto Negative;
			}
		} while (j-- > 0); //MWG-S loop

						// Transpose and Zero lower triangle
		for (j = 1; j < n; ++j)			// 0..n-1
		{
			FM::Matrix::Row UDj(UD,j);
			for (i = 0; i < j; ++i)
			{
				UD(i,j) = UDj[i];
				UDj[i] = 0.;			// Zeroing unnecessary as lower only used as a scratch
			}
		}

	}

	// Estimate the reciprocal condition number from upper triangular part
	return FM::UdUrcond(UD);

Negative:
	return -1.;
}


void 
 UD_filter::observe_size (size_t z_size)
/*
 * Optimised dyamic observation sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		s.resize(z_size);
		Sd.resize(z_size);
		zp.resize(z_size);
	}
}

Bayes_base::Float 
 UD_filter::observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z)
/* 
 * Standard linrz observe
 *  Uncorrelated observations are applied sequentialy in the order they appear in z
 *  The sequential observation updates state x
 *  Therefore the model of each observation needs to be computed sequentially. Generally this
 *  is inefficient and observe (UD_sequential_observe_model&) should be used instead
 * Predcondition:
 *		UdU' is PD (not checked)
 *		Zv is PSD
 * Postcondition:
 *		UdU' is PD
 * Return: Minimum rcond of all squential observe
 */
{
	size_t o, j;
	const size_t x_size = x.size();
	const size_t z_size = z.size();
	Float s, S;			// Innovation and covariance

								// Dynamic sizing
	observe_size (z_size);
								// Apply observations sequentialy as they are decorrelated
	Float rcondmin = std::numeric_limits<Float>::max();
	for (o = 0; o < z_size; ++o)
	{
								// Observation prediction and model
								// ISSUE inefficient need a observe_model for single z
		zp = h.h(x);			// Observation model
		h.normalise(zp, z);
		
		for (j = 0; j < x_size; ++j)
		{
			h1[j] = h.Hx(o,j);
		}
								// Check Z precondition
		if (h.Zv[o] < 0.) 
			filter_error("Zv not PSD in observe");
								// Update UD and extract gain
		Float rcond = observeUD (w, S, h1, h.Zv[o]);
		rclimit.check_PD(rcond, "X not PD in observe");
		if (rcond < rcondmin) rcondmin = rcond;
								// State update using non-linear innovation
		s = z[o]-zp[o];
		for (j = 0; j < x_size; ++j)
		{
			x[j] += w[j] * s;
		}
								// Copy s and Sd
		UD_filter::s[o] = s;
		UD_filter::Sd[o] = S;
	}
	return rcondmin;
}

Bayes_base::Float
 UD_filter::observe (Linrz_correlated_observe_model& /*h*/, const FM::Vec& /*z*/)
/* No solution for Correlated noise and Linearised model */
{
	filter_error ("observe no Linrz_correlated_observe_model solution");
	return 0.;	// never reached
}

Bayes_base::Float
 UD_filter::observe (Linear_correlated_observe_model& h, const FM::Vec& z)
/* 
 * Special Linear Hx observe for correlated Z
 *  Z must be PD and will be decorrelated
 * Applies observations sequentialy in the order they appear in z
 * Creates temporary Vec and Matrix to decorelate z,Z
 * Predcondition:
 *		UdU' is PD (not checked)
 *		Z is PSD
 * Postcondition:
 *		UdU' is PD
 * Return: Minimum rcond of all squential observe
 */
{
	size_t o, i, j, k;
	const size_t x_size = x.size();
	const size_t z_size = z.size();
	Float s, S;			// Innovation and covariance
						
					// Dynamic sizing
	if (z_size != last_z_size) {
		zdecol.resize(z_size);
		Gz.resize(z_size,z_size);
		UHx.resize(z_size, x_size);
	}
	observe_size (z_size);

					// Factorise process noise as GzG'
	{	Float rcond = FM::UdUfactor (Gz, h.Z);
		rclimit.check_PSD(rcond, "Z not PSD in observe");
	}

								// Observation prediction and model
	zp = h.h(x);
	h.normalise(zp, z);
								// Solve UHx~= Hx
	if (z_size > 0)
	{

		UHx = h.Hx;
		for (j = 0; j < x_size; ++j)
		{
			i = z_size-1;
			do {
				for (k = i+1; k < z_size; ++k)
				{
					if (i != k)
					{
						UHx(i,j) -= Gz(i,k) * UHx(k,j);
					}
					else
					{				// Unit part of U
						UHx(i,j) -= UHx(k,j);
					}
				}
			} while (i-- > 0);
		}
								// Uz~=z, Yzp~=zp  zp is in place
		zdecol.clear();
		i = z_size-1;
		do {
			for (k = i+1; k < z_size; ++k)
			{
				if (i != k)
				{
					zdecol[i] -= Gz(i,k) * z[k];
					zp[i] -= Gz(i,k) * zp[k];
				}
				else
				{				// Unit part of U
					zdecol[i] -= z[k];
					zp[i] -= zp[k];
				}
			}
		} while (i-- > 0);
	}//if (z_size>0)

								// Apply observations sequentialy as they are decorrelated
	Float rcondmin = std::numeric_limits<Float>::max();
	for (o = 0; o < z_size; ++o)
	{
		for (j = 0; j < x_size; ++j)
		{
			h1[j] = UHx(o,j);
		}
								// Update UD and extract gain
		Float rcond = observeUD (w, S, h1, Gz(o,o));
		rclimit.check_PD(rcond, "X not PD in observe linear");
		if (rcond < rcondmin) rcondmin = rcond;
								// State update using linear innovation
		s = zdecol[o]-zp[o];
		for (j = 0; j < x_size; ++j)
		{
			x[j] += w[j] * s;
		}
								// Copy s and Sd
		UD_filter::s[o] = s;
		UD_filter::Sd[o] = S;
	}
	return rcondmin;
}

Bayes_base::Float
 UD_filter::observe (UD_sequential_observe_model& h, const FM::Vec& z)
/*
 * Special observe using observe_model_sequential for fast uncorrelated linrz operation
 * Uncorrelated observations are applied sequentialy in the order they appear in z
 * The sequential observation updates state x. Therefore the model of
 * each observation needs to be computed sequentially
 * Predcondition:
 *		UdU' is PD (not checked)
 *		Z is PSD
 * Postcondition:
 *		UdU' is PD
 * Return: Minimum rcond of all squential observe
 */
{
	size_t o, j;
	const size_t x_size = x.size();
	const size_t z_size = z.size();
	Float s, S;			// Innovation and covariance

								// Dynamic sizing
	observe_size (z_size);
								// Apply observations sequentialy as they are decorrelated
	Float rcondmin = std::numeric_limits<Float>::max();
	for (o = 0; o < z_size; ++o)
	{
								// Observation prediction and model
		zp = h.ho(x, o);
		h.normalise(zp, z);
								// Check Z precondition
		if (h.Zv[o] < 0.) 
			filter_error("Zv not PSD in observe");
								// Update UD and extract gain
		Float rcond = observeUD (w, S, h.Hx_o, h.Zv[o]);
		rclimit.check_PD(rcond, "X not PD in observe");
		if (rcond < rcondmin) rcondmin = rcond;
								// State update using non-linear innovation
		s = z[o]-zp[o];
		for (j = 0; j < x_size; ++j)
		{
			x[j] += w[j] * s;
		}
								// Copy s and Sd
		UD_filter::s[o] = s;
		UD_filter::Sd[o] = S;
	}
	return rcondmin;
}


UD_filter::Float
 UD_filter::observeUD (FM::Vec& gain, Float & alpha, const FM::Vec& h, const Float r)
/*
 * Linear UD factorisation update
 *	Bierman UdU' factorisation update. Bierman p.100
 * Input
 *	h observation coeficients
 *  r observation variance
 * Output
 *	gain  observation Kalman gain
 *  alpha observation innovation variance
 * Variables with physical significance
 *  gamma becomes covariance of innovation
 * Alogirthm resrictions:
 *	Only implemented for PD matricies, assume semi-definate is negative
 * Predcondition:
 *		UdU' is PD (not checked)
 *		r is >=0 (not checked)
 * Postcondition:
 *		UdU' is PD (see return value)
 * Return:
 *		reciprocal condition number, -1. if negative or semi-definate (including zero)
 * TODO return 0 for semi-definate
 */
{
	size_t i,j,k;
	const size_t n = UD.size1();
	Float gamma, alpha_jm1, lamda;
	// a(n) is U'a
	// b(n) is Unweighted Kalman gain

					// Compute b = DU'h, a = U'h
	a = h;
	for (j = n-1; j >= 1; --j)	// n-1..1
	{
		for (k = 0; k < j; ++k)	// 0..j-1
		{
			a[j] += UD(k,j) * a[k];
		}
		b[j] = UD(j,j) * a[j];
	}
	b[0] = UD(0,0) * a[0];

					// Update UD(0,0), d(0) modification
	alpha = r + b[0] * a[0];
	if (alpha <= 0.) goto NotPD;
	gamma = 1. / alpha;
	UD(0,0) *= r * gamma;
					// Update rest of UD and gain b
	for (j = 1; j < n; ++j)		// 1..n-1
	{
					// d modification
		alpha_jm1 = alpha;	// alpha at j-1
		alpha += b[j] * a[j];
		lamda = -a[j] * gamma;
		if (alpha <= 0.) goto NotPD;
		gamma = 1. / alpha;
		UD(j,j) *= alpha_jm1 * gamma;
					// U modification
		for (i = 0; i < j; ++i)		// 0..j-1
		{
			Float UD_jm1 = UD(i,j);
			UD(i,j) = UD_jm1 + lamda * b[i];
			b[i] += b[j] * UD_jm1;
		}
	}
					// Update gain from b
	for (j = 0; j < n; ++j)		// 0..n-1
	{
		gain[j] = b[j] * gamma;
	}
 	// Estimate the reciprocal condition number from upper triangular part
	return FM::UdUrcond(UD);

NotPD:
	return -1.;
}

}//namespace
