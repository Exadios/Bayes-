/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
 * $Header$
 * $NoKeywords: $
 *
 * Unscented Filter.
 *	A Julier-Uhlmann unscented non-linear Kalman filter as an Abstract class
 * TODO
 *	Update can be done directly in factorised form
 */
#include "matSup.h"
#include <cmath>

#include "unsFlt.h"
#include "models.h"

/* Filter namespace */
namespace Bayesian_filter
{
	using namespace Bayesian_filter_matrix;


Unscented_filter::Unscented_filter (Subscript x_size, Subscript z_initialsize) :
		Extended_filter(x_size),
		XX(x_size, 2*x_size+1),
		s(Empty), S(Empty), SI(Empty),
		fXX(x_size, 2*x_size+1)
/*
 * Initialise filter and set the size of things we know about
 */
{
	Unscented_filter::x_size = x_size;
	Unscented_filter::XX_size = 2*x_size+1;
	last_z_size = 0;	// Leave z_size dependants Empty if z_initialsize==0
	observe_size (z_initialsize);
}

Unscented_filter& Unscented_filter::operator= (const Unscented_filter& a)
/* Optimise copy assignment to only copy filter state
 * Precond: matrix size conformance
 */
{
	Extended_filter::operator=(a);
	XX = a.XX;
	return *this;
}

void Unscented_filter::unscented (ColMatrix& XX, const Vec& x, const SymMatrix& X, Float Scale)
/*
 * Generate the unscented point representing a distribution
 * Fails if scale is negative
 */
{
	using namespace std;
	UTriMatrix Sigma(x_size,x_size);

						// Get a upper Cholesky factoriation
	Float rcond = UCfactor(Sigma, X);
	rclimit.check_PSD(rcond, "X not PSD");
	Sigma *= sqrt(Scale);

						// Generate XX with the same sample Mean and Covar as before
	XX[0] = x;

	Vec SigmaCol(x_size);
	for (Subscript c = 0; c < x_size; ++c) {
		// ISSUE: Copy a column
		// uBLAS	SigmaCol = UTriMatrix::Column(Sigma,c);
		// MTL		column(sigma,c) may be buggy
		Subscript r;
		for (r = 0; r <= c; ++r)
			SigmaCol[r] = Sigma(r,c);
		for (; r < x_size; ++r)
			SigmaCol[r] = 0.;
		XX[c+1] = x + SigmaCol;
		XX[x_size+c+1] = x - SigmaCol;
	}
}

Unscented_filter::Float Unscented_filter::predict_Kappa (unsigned size) const
// Default Kappa for prediction: state augmented with predict noise
{
	// Use the rule to minimise mean squared error of 4 order term
	return Float(3-signed(size));
}

Unscented_filter::Float Unscented_filter::observe_Kappa (unsigned size) const
// Default Kappa for observation: state on its own
{
	// Use the rule to minimise mean squared error of 4 order term
	return Float(3-signed(size));
}

void Unscented_filter::init ()
/*
 * Initialise unscented state
 *		Pre : x,X
 *		Post: x,X
 */
{
					// Precondition that unscented distribution can be created
	unscented (XX, x, X, 0.);		// Kappa can be zero for check
}

void Unscented_filter::update ()
/*
 * Update state varibles
 *		Pre : x,X
 *		Post: x,X
 */
{
	// TODO: Remove check once proved
	assert_isPSD (X);
}

void Unscented_filter::predict (Functional_predict_model& f)
/*
 * Adapt Unscented noise model by creating a noise model with zero noise
 * ISSUE: A simple specialisation is possible, rather then this adapted implemenation
 */
{
	class Adapted_model : public Unscented_predict_model
	{
	public:
		Adapted_model(Functional_predict_model& fm) :
			Unscented_predict_model(0),
			fmodel(fm), zeroQ(0,0)
		{}
		const Vec& f(const Vec& x) const
		{
			return fmodel.fx(x);
		}
		const SymMatrix& Q(const FM::Vec& /*x*/) const
		{
			return zeroQ;
		}
	private:
		Functional_predict_model& fmodel;
		SymMatrix zeroQ;
	};
	Adapted_model adaptedmodel(f);
	predict (adaptedmodel);
}

void Unscented_filter::predict (Addative_predict_model& f)
/*
 * Adapt Unscented noise model by creating a noise model with zero noise
 *  Computes noise covariance Q = GqG'
 */
{
	class Adapted_model : public Unscented_predict_model
	{
	public:
		Adapted_model(Subscript x_size, Addative_predict_model& am) :
			Unscented_predict_model(0),
			amodel(am), QGqG(x_size,x_size)		// Q gets size from GqG'
		{
			QGqG.clear();
			mult_SPD (am.G, am.q, QGqG);
		}
		const Vec& f(const Vec& x) const
		{
			return amodel.f(x);
		}
		const SymMatrix& Q(const FM::Vec& /*x*/) const
		{
			return QGqG;
		}
	private:
		Addative_predict_model& amodel;
		mutable SymMatrix QGqG;
	};

	Adapted_model adaptedmodel(f.G.size1(), f);
	predict (adaptedmodel);
}


void Unscented_filter::predict (Unscented_predict_model& f)
/*
 * Predict forward
 *		Pre : x,X represent the prior distribution
 *		Post: x,X represent the predicted distribution
 * Implementation uses specific model for fast unscented compuation
 */
{
	Subscript i;
	const Subscript XX_size = XX.size2();
				
						// Create unscented distribution
	Float Kappa = predict_Kappa(x_size);
	Float x_Kappa = Float(x_size) + Kappa;
	unscented (XX, x, X, x_Kappa);

						// Predict points of XX using supplied predict model
							// State covariance
	for (i = 0; i < XX_size; ++i) {
		fXX[i].assign (f.f(XX[i]));
	}
						
						// Mean of predicted distribution: x
	x.assign (ublas::operator*(fXX[0], Kappa));	// VC6 ambiguity workaround
	for (i = 1; i < XX_size; ++i) {
		x.plus_assign (fXX[i] / 2.0);
	}
	x /= x_Kappa;
						// Covariance of distribution: X
							// Subtract mean from each point in fXX
	for (i = 0; i < XX_size; ++i) {
		fXX[i].minus_assign (x);
	}
							// Center point, premult here by 2 for efficency
	X.clear();
	X.plus_assign (FM::outer_prod(fXX[0], fXX[0]));
	X *= 2.*Kappa;
							// Remaining unscented points
	for (i = 1; i < XX_size; ++i) {
		X.plus_assign (FM::outer_prod(fXX[i], fXX[i]));
	}
	X /= 2.*x_Kappa;
						// Addative Noise Prediction, computed about center point
	X.plus_assign (f.Q(fXX[0]));

	// TODO: Remove check once proved
	assert_isPSD (X);
}


void Unscented_filter::observe_size (Subscript z_size)
/*
 * Optimised dynamic observation sizing
 */
{
	if (z_size != last_z_size) {
		last_z_size = z_size;

		s.resize(z_size);
		S.resize(z_size,z_size);
		SI.resize(z_size,z_size);
	}
}


Bayes_base::Float Unscented_filter::observe (Uncorrelated_addative_observe_model& h, const Vec& z)
/*
 * Observation fusion
 *		Pre : x,X represent the predicted distribution
 *		Post: x,X represent the fused distribution
 *
 * Uncorrelated noise
 * ISSUE: Simplify implmenation using uncorrelated noise equations
 */
{
	Adapted_Correlated_addative_observe_model hh(h);
	return observe (hh, z);
}


Bayes_base::Float Unscented_filter::observe (Correlated_addative_observe_model& h, const Vec& z)
/*
 * Observation fusion
 *		Pre : x,X represent the predicted distribution
 *		Post: x,X represent the fused distribution
 */
{
	Subscript i;
	Subscript z_size = z.size();
	ColMatrix zXX (z_size, 2*x_size+1);
	Vec zp(z_size);
	SymMatrix Xzz(z_size,z_size);
	Matrix Xxz(x_size,z_size);
	Matrix W(x_size,z_size);

	observe_size (z.size());	// Dynamic sizing

						// Create unscented distribution
	Float Kappa = observe_Kappa(x_size);
	Float x_Kappa = Float(x_size) + Kappa;
	unscented (XX, x, X, x_Kappa);

						// Predict points of XX using supplied observation model
	{
		Vec zXXi(z_size), zXX0(z_size);
		for (i = 0; i < XX.size2(); ++i) {
			zXXi = h.h(XX[i]);
						// Normalise relative to zXX0 which is normalised to z
			if (i > 0) {
				h.normalise (zXXi, zXX0);
			}
			else {
				h.normalise (zXXi, z);
				zXX0 = zXXi;		// Normalise to zXX0
			}
			zXX[i].assign (zXXi);
		}
	}

						// Mean of predicted distribution: zp
	zp.assign (ublas::operator*(zXX[0], Kappa));	// VC6 ambiguity workaround
	for (i = 1; i < zXX.size2(); ++i) {
		zp.plus_assign (zXX[i] / 2.0);
	}
	zp /= x_Kappa;

						// Covariance of observation prediction: Xzz
							// Subtract mean from each point in zX
	for (i = 0; i < XX_size; ++i) {
		zXX[i].minus_assign (zp);
	}
							// Center point, premult here by 2 for efficency
	{
		Vec tzXX0 = zXX[0];		// TODO Make this a reference
		Xzz.clear();
		Xzz.plus_assign (FM::outer_prod(tzXX0, tzXX0));
		Xzz *= 2.*Kappa;
	}
							// Remaining unscented points
	for (i = 1; i < zXX.size2(); ++i) {
		Vec tzXXi = zXX[i];		// TODO Make this a reference
		Xzz.plus_assign (FM::outer_prod(tzXXi, tzXXi));
	}
	Xzz /= 2.*x_Kappa;

						// Correlation of state with observation: Xxz
							// Center point, premult here by 2 for efficency
	Xxz.clear();
	Xxz.plus_assign (FM::outer_prod(XX[0] - x, zXX[0]));
	Xxz *= 2.*Kappa;
							// Remaining unscented points
	for (i = 1; i < zXX.size2(); ++i) {
		Xxz.plus_assign (FM::outer_prod(XX[i] -x, zXX[i]));
	}
	Xxz /= 2.* (Float(x_size) + Kappa);

						// Innovation covariance
	S = Xzz;
	S.plus_assign (h.Z);
						// Inverse innovation covariance
	Float rcond = UdUinversePD (SI, S);
	rclimit.check_PD(rcond, "S not PD in observe");
						// Kalman gain
	W.assign (prod(Xxz,SI));

						// Filter update
	s = z; s.minus_assign (zp);
	x.plus_assign (prod(W,s));
	X.minus_assign (mult_SPD(W, S));

	// TODO: Remove check once proved
	assert_isPSD (X);

	return rcond;
}


}//namespace
