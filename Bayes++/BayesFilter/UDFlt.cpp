/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * UdU' Factorisation of Covariance Filter.
 *
 * For efficiency UD_scheme requires to know the maximum q_size of the predict_model
 * ISSUES:
 *  observe functions: returned rcond is the minimum of each sequential update, an overall conditioning would be better
 */
#include "UDFlt.hpp"
#include "matSup.hpp"
#include <boost/limits.hpp>

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


UD_bscheme::
UD_bscheme (std::size_t x_size, std::size_t q_maxsize) :
		Kalman_state(x_size),
		q_max(q_maxsize),
		UD(x_size,x_size+q_max)
/*
 * Initialise filter and set the size of things we know about
 */
{}

UD_scheme::
UD_scheme (std::size_t x_size, std::size_t q_maxsize) :
		Kalman_state(x_size),
		UD_bscheme(x_size,q_maxsize)
{}

UD_bscheme::
Predict_byproduct::Predict_byproduct (std::size_t x_size, std::size_t q_size) :
	d(x_size+q_size), dv(x_size+q_size), v(x_size+q_size)
{}

UD_bscheme::
Observe_innovation_byproduct::Observe_innovation_byproduct (std::size_t x_size, std::size_t z_size) :
		s(z_size),
		Sv(z_size),
		W(x_size, z_size)
{}

UD_bscheme::
Observe_byproduct::Observe_byproduct (std::size_t x_size, std::size_t z_size) :
		a(x_size),
		w(x_size),
		znorm(z_size)
{}

UD_bscheme::
Observe_linear_byproduct::Observe_linear_byproduct (std::size_t x_size, std::size_t z_size) :
		Observe_byproduct(x_size, z_size),
		zpdecol(z_size),
		Gz(z_size, z_size),
		GIHx(z_size, x_size)
{}

UD_bscheme&
 UD_bscheme::operator= (const UD_bscheme& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Kalman_state::operator=(a);
	q_max = a.q_max;
	UD = a.UD;
	return *this;
}


void
 UD_bscheme::init ()
/*
 * Initialise from a state and state coveriance
 * Computes UD factor from initial covaiance
 * Predcond:
 *  X
 * Postcond:
 *  X
 *  UD=X, d is PSD
 */
{
					// Factorise X into left partition of UD
	std::size_t x_size = UD.size1();
	noalias(UD.sub_matrix(0,x_size, 0,x_size)) = X;
	Float rcond = UdUfactor (UD, x_size);
	rclimit.check_PSD(rcond, "Initial X not PSD");
}


void
 UD_bscheme::update ()
/*
 * Defactor UD back into X
 * Precond:
 *  UD
 * Postcond:
 *  X=UD  PSD iff UD is PSD
 */
{
	UdUrecompose (X, UD);
}


UD_scheme::Float
 UD_scheme::predict (Linrz_predict_model& f)
/* Linrz_kalman_filter predict, unused byproduct
 */
{
	Predict_byproduct b(f.G.size1(), f.G.size2());
	return bypredict (f, b);
}

Bayes_base::Float
 UD_bscheme::bypredict (Linrz_predict_model& f, Predict_byproduct& b)
/*
 * Predict using a diagonalised noise q, and its coupling G
 *  q can have order less then x and a matching G so GqG' has order of x
 * Precond:
 *	UD
 * Postcond:
 *  UD is PSD
 */
{
	x = f.f(x);			// Extended Kalman state predict is f(x) directly

						// Predict UD from model
	Float rcond = predictGq (f.Fx, f.G, f.q, b);
	rclimit.check_PSD(rcond, "X not PSD in predict");
	return rcond;
}


Bayes_base::Float
 UD_bscheme::predictGq (const Matrix& Fx, const Matrix& G, const Vec& q, Predict_byproduct& b)
/*
 * MWG-S predict from Bierman  p.132
 *  q can have order less then x and a matching G so GqG' has order of x
 * Precond:
 *  UD
 * Postcond:
 *  UD
 *
 * Return:
 *		reciprocal condition number, -1 if negative, 0 if semi-definite (including zero)
 */
{
	std::size_t i,j,k;
	const std::size_t n = x.size();
	const std::size_t Nq = q.size();
	const std::size_t N = n+Nq;
	Float e;
					// Check preallocated space for q size
	if (Nq > q_max)
		error (Logic_exception("Predict model q larger than preallocated space"));

	if (n > 0)		// Simplify reverse loop termination
	{
						// Augment d with q, UD with G
		for (i = 0; i < Nq; ++i)		// 0..Nq-1
		{
			b.d[i+n] = q[i];
		}
		for (j = 0; j < n; ++j)		// 0..n-1
		{
			Matrix::Row UDj(UD,j);
			Matrix::const_Row  Gj(G,j);
			for (i = 0; i < Nq; ++i)		// 0..Nq-1
				UDj[i+n] = Gj[i];
		}

						// U=Fx*U and diagonals retrived
		for (j = n-1; j > 0; --j)		// n-1..1
		{
						// Prepare d(0)..d(j) as temporary
			for (i = 0; i <= j; ++i)	// 0..j
				b.d[i] = Float(UD(i,j));	// ISSUE mixed type proxy assignment

						// Lower triangle of UD is implicity empty
			for (i = 0; i < n; ++i) 	// 0..n-1
			{
				Matrix::Row UDi(UD,i);
				Matrix::const_Row Fxi(Fx,i);
				UDi[j] = Fxi[j];
				for (k = 0; k < j; ++k)	// 0..j-1
					UDi[j] += Fxi[k] * b.d[k];
			}
		}
		b.d[0] = Float(UD(0,0));	// ISSUE mixed type proxy assignment

						//  Complete U = Fx*U
		for (j = 0; j < n; ++j)			// 0..n-1
		{
			UD(j,0) = Fx(j,0);
		}

						// The MWG-S algorithm on UD transpose
		j = n-1;
		do {							// n-1..0
			Matrix::Row UDj(UD,j);
			e = 0;
			for (k = 0; k < N; ++k)		// 0..N-1
			{
				b.v[k] = Float(UDj[k]);	// ISSUE mixed type proxy assignment
				b.dv[k] = b.d[k] * b.v[k];
				e += b.v[k] * b.dv[k];
			}
			// Check diagonal element
			if (e > 0)
			{
				// Positive definite
				UDj[j] = e;

				Float diaginv = 1 / e;
				for (k = 0; k < j; ++k)	// 0..j-1
				{
					Matrix::Row UDk(UD,k);
					e = 0;
					for (i = 0; i < N; ++i)	// 0..N-1
						e += UDk[i] * b.dv[i];
					e *= diaginv;
					UDj[k] = e;

					for (i = 0; i < N; ++i)	// 0..N-1
						UDk[i] -= e * b.v[i];
				}
			}//PD
			else if (e == 0)
			{
				// Possibly semi-definite, check not negative
				UDj[j] = e;

				// 1 / e is infinite
				for (k = 0; k < j; ++k)	// 0..j-1
				{
					Matrix::Row UDk(UD,k);
					for (i = 0; i < N; ++i)	// 0..N-1
					{
						e = UDk[i] * b.dv[i];
						if (e != 0)
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
			Matrix::Row UDj(UD,j);
			for (i = 0; i < j; ++i)
			{
				UD(i,j) = UDj[i];
				UDj[i] = 0;			// Zeroing unnecessary as lower only used as a scratch
			}
		}

	}

	// Estimate the reciprocal condition number from upper triangular part
	return UdUrcond(UD,n);

Negative:
	return -1;
}


Bayes_base::Float
 UD_scheme::observe (Linrz_uncorrelated_observe_model& h, const Vec& z)
/* Linrz_kalman_filter observe
 */
{
	Observe_innovation_byproduct g(h.Hx.size2(), h.Hx.size1());
	Observe_byproduct b(h.Hx.size2(), h.Hx.size1());
	return byobserve (h, z, g, b);
}

Bayes_base::Float
 UD_bscheme::byobserve (Linrz_uncorrelated_observe_model& h, const Vec& z, Observe_innovation_byproduct& g, Observe_byproduct& b)
/*
 * Linrz observe with byproduct
 *  Uncorrelated observations are applied sequentialy in the order they appear in z
 *  The sequential observation updates state x
 *  Therefore the model of each observation needs to be computed sequentially. Generally this
 *  is inefficient and observe (UD_sequential_observe_model&) should be used instead
 * Precond:
 *  UD
 *  Zv is PSD
 * Postcond:
 *  UD is PSD
 * Return: Minimum rcond of all squential observe
 */
{
	const std::size_t z_size = z.size();
	Float s, S;			// Innovation and covariance

								// Apply observations sequentialy as they are decorrelated
	Float rcondmin = std::numeric_limits<Float>::max();
	for (std::size_t o = 0; o < z_size; ++o)
	{
								// Observation model, extracted for a single z element
		const Vec& zp = h.h(x);
		h.normalise(b.znorm = z, zp);
		noalias(b.a) = row(h.Hx,o);
								// Check Z precondition
		if (h.Zv[o] < 0)
			error (Numeric_exception("Zv not PSD in observe"));
								// Update UD and extract gain
		Float rcond = observeUD (h.Zv[o], b.a, b.w, S);
		rclimit.check_PSD(rcond, "S not PD in observe");	// -1 implies S singular
		if (rcond < rcondmin) rcondmin = rcond;
								// State update using normalised non-linear innovation
		s = b.znorm[o] - zp[o];
		noalias(x) += b.w * s;
								// Innovation and gain byproducts
		g.s[o] = s;
		g.Sv[o] = S;
		noalias(column(g.W, o)) = b.w;
	}
	return rcondmin;
}

Bayes_base::Float
 UD_scheme::observe (Linrz_correlated_observe_model& /*h*/, const Vec& /*z*/)
/* NO solution for Correlated noise and Linearised model */
{
	error (Logic_exception("observe no Linrz_correlated_observe_model solution"));
	return 0;	// never reached
}

Bayes_base::Float
 UD_bscheme::byobserve (Linear_correlated_observe_model& h, const Vec& z, Observe_innovation_byproduct& g, Observe_linear_byproduct& b)
/*
 * Special Linear Hx observe for correlated Z
 *  Z must be PD and will be decorrelated
 * Applies observations sequentialy in the order they appear in z
 * Creates temporary Vec and Matrix to decorelate z,Z
 * Predcondition:
 *  UD
 *  Z is PSD
 * Postcondition:
 *  UD is PSD
 * Return: Minimum rcond of all squential observe
 */
{
	std::size_t i, j, k;
	const std::size_t x_size = x.size();
	const std::size_t z_size = z.size();
	Float s, S;			// Innovation and covariance

					// Factorise process noise as GzG'
	{	Float rcond = UdUfactor (b.Gz, h.Z);
		rclimit.check_PSD(rcond, "Z not PSD in observe");
	}

								// Observation predict and normalised observation
	const Vec& zp = h.h(x);
	h.normalise(b.znorm = z, zp);
	
	if (z_size > 0)
	{							// Solve G* GIHx = Hx for GIHx in-place
		b.GIHx = h.Hx;
		for (j = 0; j < x_size; ++j)
		{
			i = z_size-1;
			do {
				for (k = i+1; k < z_size; ++k)
				{
					b.GIHx(i,j) -= b.Gz(i,k) * b.GIHx(k,j);
				}
			} while (i-- > 0);
		}
					
		b.zpdecol = zp;			// Solve G zp~ = z, G z~ = z  for zp~,z~ in-place
		i = z_size-1;
		do {
			for (k = i+1; k < z_size; ++k)
			{
				b.znorm[i] -= b.Gz(i,k) * b.znorm[k];
				b.zpdecol[i] -= b.Gz(i,k) * b.zpdecol[k];
			}
		} while (i-- > 0);
	}//if (z_size>0)

								// Apply observations sequentialy as they are decorrelated
	Float rcondmin = std::numeric_limits<Float>::max();
	for (std::size_t o = 0; o < z_size; ++o)
	{
		noalias(b.a) = row(b.GIHx,o);
								// Update UD and extract gain
		Float rcond = observeUD (b.Gz(o,o), b.a, b.w, S);
		rclimit.check_PSD(rcond, "S not PD in observe");	// -1 implies S singular
		if (rcond < rcondmin) rcondmin = rcond;
								// State update using linear innovation
		s = b.znorm[o]-b.zpdecol[o];
		noalias(x) += b.w * s;
								// Innovation and gain byproducts
		g.s[o] = s;
		g.Sv[o] = S;
		noalias(column(g.W, o)) = b.w;
	}
	return rcondmin;
}

Bayes_base::Float
 UD_bscheme::byobserve (Sequential_observe_model& h, const Vec& z, Observe_innovation_byproduct& g, Observe_byproduct& b)
/*
 * Special observe using observe_model_sequential for fast uncorrelated linrz operation
 * Uncorrelated observations are applied sequentialy in the order they appear in z
 * The sequential observation updates state x. Therefore the model of
 * each observation needs to be computed sequentially
 * Predcondition:
 *  UD
 *  Z is PSD
 * Postcondition:
 *  UD is PSD
 * Return: Minimum rcond of all squential observe
 */
{
	std::size_t o;
	const std::size_t z_size = z.size();
	Float s, S;			// Innovation and covariance

								// Apply observations sequentialy as they are decorrelated
	Float rcondmin = std::numeric_limits<Float>::max();
	for (o = 0; o < z_size; ++o)
	{
								// Observation predict and model
		Float Zv_o;
		const Vec& zp = h.ho(x, o, Zv_o, b.a);
		h.normalise(b.znorm = z, zp);
								// Check Z precondition
		if (Zv_o < 0)
			error (Numeric_exception("Zv not PSD in observe"));
								// Update UD and extract gain
		Float rcond = observeUD (Zv_o, b.a, b.w, S);
		rclimit.check_PSD(rcond, "S not PD in observe");	// -1 implies S singular
		if (rcond < rcondmin) rcondmin = rcond;
								// State update using non-linear innovation
		s = b.znorm[o]-zp[o];
		noalias(x) += b.w * s;
								// Innovation and gain byproducts
		g.s[o] = s;
		g.Sv[o] = S;
		noalias(column(g.W, o)) = b.w;
	}
	return rcondmin;
}


Bayes_base::Float
 UD_bscheme::observeUD (const Float r, Vec& a, Vec& b, Float& alpha)
/*
 * Linear UD factorisation update
 *  Bierman UdU' factorisation update. Bierman p.100
 * Input
 *  a linrz observation model (h)
 *  r observation variance
 * Output
 *  a  U'h
 *  b  observation Kalman gain
 *  alpha observation innovation variance
 * Predcondition:
 *  UD
 *  r is PSD (not checked)
 * Postcondition:
 *  UD (see return value)
 * Return:
 *  reciprocal condition number of UD, -1 if alpha singular (negative or zero)
 */
{
	std::size_t i,j,k;
	const std::size_t n = UD.size1();
	Float gamma;	// becomes inverse covariance of innovation
	Float alpha_jm1, lamda;

					// Compute b = DU'h, a = U'h
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
	if (alpha <= 0) goto alphaNotPD;
	gamma = 1 / alpha;
	UD(0,0) *= r * gamma;
					// Update rest of UD and gain b
	for (j = 1; j < n; ++j)		// 1..n-1
	{
					// d modification
		alpha_jm1 = alpha;	// alpha at j-1
		alpha += b[j] * a[j];
		lamda = -a[j] * gamma;
		if (alpha <= 0) goto alphaNotPD;
		gamma = 1 / alpha;
		UD(j,j) *= alpha_jm1 * gamma;
					// U modification
		for (i = 0; i < j; ++i)		// 0..j-1
		{
			Float UD_jm1 = UD(i,j);
			UD(i,j) = UD_jm1 + lamda * b[i];
			b[i] += b[j] * UD_jm1;
		}
	}
					// weighted gain
	b *= gamma;
 	// Estimate the reciprocal condition number from upper triangular part
	return UdUrcond(UD,n);

alphaNotPD:
	return -1;
}

}//namespace
