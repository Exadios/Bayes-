/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Unscented Filter.
 */
#include "unsFlt.hpp"
#include "matSup.hpp"
#include "models.hpp"
#include <cmath>


/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


/**
 * Initialise filter and set the size of things we know about
 */
Unscented_scheme::Unscented_scheme (size_t x_size) :
		Kalman_state_filter(x_size), Functional_filter(),
		XX(x_size, 2*x_size+1),
		fXX(x_size, 2*x_size+1)
{
	Unscented_scheme::x_size = x_size;
	Unscented_scheme::XX_size = 2*x_size+1;
}

/** Optimise copy assignment to only copy filter state.
 * @pre matrix size conformance
 */
Unscented_scheme& Unscented_scheme::operator= (const Unscented_scheme& a)
{
	Kalman_state_filter::operator=(a);
	XX = a.XX;
	return *this;
}

/*
 * Generate the unscented point representing a distribution.
 * Fails if scale is negative
 */
void Unscented_scheme::unscented (ColMatrix& XX, const Vec& x, const SymMatrix& X, Float scale)
{
	UTriMatrix Sigma(x_size,x_size);

						// Get a upper Cholesky factoriation
	Float rcond = UCfactor(Sigma, X);
	rclimit.check_PSD(rcond, "X not PSD");
	Sigma *= std::sqrt(scale);

						// Generate XX with the same sample Mean and Covar as before
	column(XX,0) = x;

	for (size_t c = 0; c < x_size; ++c) {
		UTriMatrix::Column SigmaCol = column(Sigma,c);
		noalias(column(XX,c+1)) = x  + SigmaCol;
		noalias(column(XX,x_size+c+1)) = x - SigmaCol;
	}
}

Unscented_scheme::Float Unscented_scheme::predict_Kappa (size_t size) const
// Default Kappa for predict: state augmented with predict noise
{
	// Use the rule to minimise mean squared error of 4 order term
	return Float(3-signed(size));
}

Unscented_scheme::Float Unscented_scheme::observe_Kappa (size_t size) const
// Default Kappa for observation: state on its own
{
	// Use the rule to minimise mean squared error of 4 order term
	return Float(3-signed(size));
}

/** Initialise Kalman state.
 * @pre x,X
 * @post x,X is PSD
 */
void Unscented_scheme::init ()
{
						// Postconditions
	if (!isPSD (X))
		error (Numeric_exception("Initial X not PSD"));
}

/** Initialise from Unscented state.
 * @pre XX, kappa
 * @post x,X is PSD
 */
void Unscented_scheme::init_XX ()
{
	Float x_kappa = Float(x_size) + kappa;
						// Mean of predicted distribution: x
	noalias(x) = column(fXX,0) * kappa;
	for (size_t i = 1; i < XX_size; ++i) {
		noalias(x) += column(fXX,i) / Float(2); // ISSUE uBlas may not be able to promote integer 2
	}
	x /= x_kappa;
						// Covariance of distribution: X
							// Subtract mean from each point in fXX
	for (size_t i = 0; i < XX_size; ++i) {
		noalias(column(fXX,i)) -= x;
	}
							// Center point, premult here by 2 for efficency
    {
		ColMatrix::Column fXX0 = column(fXX,0);
		noalias(X) = FM::outer_prod(fXX0, fXX0);
		X *= 2*kappa;
	}
							// Remaining unscented points
	for (size_t i = 1; i < XX_size; ++i) {
		ColMatrix::Column fXXi = column(fXX,i);
		noalias(X) += FM::outer_prod(fXXi, fXXi);
	}
	X /= 2*x_kappa;
}

/** Update Kalman state.
 * @pre x,X
 * @post x,X
 */
void Unscented_scheme::update ()
{
}

/** Update Unscented state.
 * @pre x,X
 * @post x,X, XX, kappa
 */
void Unscented_scheme::update_XX (Float kappa)
{
	Unscented_scheme::kappa = kappa;
	Float x_kappa = Float(x_size) + kappa;
	unscented (XX, x, X, x_kappa);
}


// ISSUE GCC2.95 cannot link if these are localy defined in member function
// Move them back into member functions for standard compilers
namespace {
	class Adapted_zero_model : public Unscented_predict_model
	{
	public:
		Adapted_zero_model(Functional_predict_model& fm) :
			Unscented_predict_model(),
			fmodel(fm), zeroQ(0,0)
		{}
		const Vec& f(const Vec& x) const
		{
			return fmodel.fx(x);
		}
		const SymMatrix& Q(const Vec& /*x*/) const
		{
			return zeroQ;
		}
	private:
		Functional_predict_model& fmodel;
		SymMatrix zeroQ;
	};

	class Adapted_model : public Unscented_predict_model
	{
	public:
		Adapted_model(Additive_predict_model& am) :
			Unscented_predict_model(),
			amodel(am), QGqG(am.G.size1(),am.G.size1())		// Q gets size from GqG'
		{
			noalias(QGqG) = prod_SPD(am.G, am.q);
		}
		const Vec& f(const Vec& x) const
		{
			return amodel.f(x);
		}
		const SymMatrix& Q(const Vec& /*x*/) const
		{
			return QGqG;
		}
	private:
		Additive_predict_model& amodel;
		mutable SymMatrix QGqG;
	};
}//namespace


/** Adapt model by creating an Unscented predict with zero noise.
 * ISSUE: A simple specialisation is possible, rather then this adapted implemenation
 */
void Unscented_scheme::predict (Functional_predict_model& f)
{
	Adapted_zero_model adaptedmodel(f);
	predict (adaptedmodel);
}


/** Adapt model by creating an Unscented predict with additive noise.
 * Computes noise covariance Q = GqG'
 */
void Unscented_scheme::predict (Additive_predict_model& f)
{
	Adapted_model adaptedmodel(f);
	predict (adaptedmodel);
}


/** Predict forward.
 * @pre x,X
 * @post x,X is PSD
 * 
 * Implementation uses specific model for fast unscented computation
 */
void Unscented_scheme::predict (Unscented_predict_model& f)
{
	const size_t XX_size = XX.size2();

						// Create unscented distribution
	kappa = predict_Kappa(x_size);
	Float x_kappa = Float(x_size) + kappa;
	unscented (XX, x, X, x_kappa);

						// Predict points of XX using supplied predict model
							// State covariance
	for (size_t i = 0; i < XX_size; ++i) {
		noalias(column(fXX,i)) = f.f( column(XX,i) );
	}

	init_XX ();
						// Additive Noise predict, computed about center point
	noalias(X) += f.Q( column(fXX,0) );
}


/** Observation fusion with uncorrelated noise.
 * @pre x,X
 * @post x,X is PSD
 *
 * ISSUE: Simplified implemenation using uncorrelated noise equations
 */
Bayes_base::Float Unscented_scheme::eobserve (Uncorrelated_additive_observe_model& h, const Vec& z,
				State_byproduct& s, Covariance_byproduct& S, Kalman_gain_byproduct& b)
{
	Adapted_Correlated_additive_observe_model hh(h);
	return eobserve (hh, z, s, S, b);
}


/** Observation fusion with correlated noise.
 * @pre x,X
 * @post x,X is PSD
 */
Bayes_base::Float Unscented_scheme::eobserve (Correlated_additive_observe_model& h, const Vec& z,
				State_byproduct& s, Covariance_byproduct& S, Kalman_gain_byproduct& b)
{
	size_t z_size = z.size();
	ColMatrix zXX (z_size, 2*x_size+1);
	Vec zp(z_size);
	SymMatrix Xzz(z_size,z_size);
	Matrix Xxz(x_size,z_size);

						// Create unscented distribution
	kappa = observe_Kappa(x_size);
	Float x_kappa = Float(x_size) + kappa;
	unscented (XX, x, X, x_kappa);

						// Predict points of XX using supplied observation model
	{
		Vec zXXi(z_size), zXX0(z_size);
		zXX0 = h.h( column(XX,0) );
		column(zXX,0) = zXX0;
		for (size_t i = 1; i < XX.size2(); ++i) {
			zXXi = h.h( column(XX,i) );
						// Normalise relative to zXX0
			h.normalise (zXXi, zXX0);
			column(zXX,i) = zXXi;
		}
	}

						// Mean of predicted distribution: zp
	noalias(zp) = column(zXX,0) * kappa;
	for (size_t i = 1; i < zXX.size2(); ++i) {
		noalias(zp) += column(zXX,i) / Float(2); // ISSUE uBlas may not be able to promote integer 2
	}
	zp /= x_kappa;

						// Covariance of observation predict: Xzz
							// Subtract mean from each point in zXX
	for (size_t i = 0; i < XX_size; ++i) {
		column(zXX,i).minus_assign (zp);
	}
							// Center point, premult here by 2 for efficency
	{
		ColMatrix::Column zXX0 = column(zXX,0);
		noalias(Xzz) = FM::outer_prod(zXX0, zXX0);
		Xzz *= 2*kappa;
	}
							// Remaining unscented points
	for (size_t i = 1; i < zXX.size2(); ++i) {
		ColMatrix::Column zXXi = column(zXX,i);
		noalias(Xzz) += FM::outer_prod(zXXi, zXXi);
	}
	Xzz /= 2*x_kappa;

						// Correlation of state with observation: Xxz
							// Center point, premult here by 2 for efficency
	{
		noalias(Xxz) = FM::outer_prod(column(XX,0) - x, column(zXX,0));
		Xxz *= 2*kappa;
	}
							// Remaining unscented points
	for (size_t i = 1; i < zXX.size2(); ++i) {
		noalias(Xxz) += FM::outer_prod(column(XX,i) - x, column(zXX,i));
	}
	Xxz /= 2* (Float(x_size) + kappa);

						// Innovation covariance
	S = Xzz;
	noalias(S) += h.Z;
						// Inverse innovation covariance
	Float rcond = UdUinversePD (b.SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");
						// Kalman gain
	noalias(b.W) = prod(Xxz,b.SI);

						// Normalised innovation
	h.normalise(s = z, zp);
	noalias(s) -= zp;

						// Filter update
	noalias(x) += prod(b.W,s);
	RowMatrix WStemp(b.W.size1(), S.size2());
	noalias(X) -= prod_SPD(b.W,S, WStemp);

	return rcond;
}


}//namespace
