/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
 * $Header$
 * $NoKeywords: $
 */
 
/*
 * Bayesian_filter Implemention:
 *  algorithm that require additional resources to implement
 */
#include "matSup.h"
#include <vector>		// Only for unique_samples

#include "bayesFlt.h"
#include "models.h"

namespace {

template <class scalar>
inline scalar sqr(scalar x)
// Square 
{
	return x*x;
}
};//namespace


/* Filter namespace */
namespace Bayesian_filter
{

Bayes_base::Float
 Extended_filter::observe (Linrz_correlated_observe_model& h, const FM::Vec& z)
/*
 * Extended linrz correlated observe, compute innovation for observe_innovation
 */
{
	FM::Vec zp = h.h(x);	// Observation model
	h.normalise(zp, z);

	zp.assign (z-zp);	// Innovation
	return observe_innovation (h, zp);
}

Bayes_base::Float
 Extended_filter::observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z)
/*
 * Extended linrz uncorrelated observe, compute innovation for observe_innovation
 */
{
	FM::Vec zp = h.h(x);	// Observation model
	h.normalise(zp, z);

	zp.assign (z-zp);	// Innovation
	return observe_innovation (h, zp);
}


Simple_addative_predict_model::Simple_addative_predict_model (Function_model& f_init, const FM::Matrix& G_init, const FM::Vec& q_init) :
	Addative_predict_model (G_init.size1(), q_init.size()),
	ff(f_init)
/* Addative predict model initialised from function and model matricies
   Precondition:
	q and G are conformantly dimensioned (not checked)
 */
{
	G = G_init;
	q = q_init;
}

Simple_linrz_predict_model::Simple_linrz_predict_model (Function_model& f_init, const FM::Matrix& Fx_init, const FM::Matrix& G_init, const FM::Vec& q_init) :
	Linrz_predict_model (Fx_init.size1(), q_init.size()),
	ff(f_init)
/* Linrz predict model initialised from function and model matricies
   Precondition:
	Fx, q and G are conformantly dimensioned (not checked)
 */
{
	Fx = Fx_init;
	G = G_init;
	q = q_init;
}

Simple_linear_predict_model::Simple_linear_predict_model (const FM::Matrix& Fx_init, const FM::Matrix& G_init, const FM::Vec& q_init) :
	Linear_predict_model (Fx_init.size1(), q_init.size())
/* Linear predict model initialised from model matricies
   Precondition:
	Fx, q and G are conformantly dimensioned (not checked)
 */
{
	Fx = Fx_init;
	G = G_init;
	q = q_init;
}

Simple_linrz_correlated_observe_model::Simple_linrz_correlated_observe_model (Function_model& f_init, const FM::Matrix& Hx_init, const FM::SymMatrix& Z_init) :
	Linrz_correlated_observe_model (Hx_init.size2(), Hx_init.size1()),
	ff(f_init)
/* Linrz observe model initialised from function and model matricies
   Precondition:
   Hx, Z are conformantly dimensioned (not checked)
 */
{
	Hx = Hx_init;
	Z = Z_init;
};

Simple_linrz_uncorrelated_observe_model::Simple_linrz_uncorrelated_observe_model (Function_model& f_init, const FM::Matrix& Hx_init, const FM::Vec& Zv_init) :
	Linrz_uncorrelated_observe_model (Hx_init.size2(), Hx_init.size1()),
	ff(f_init)
/* Linrz observe model initialised from function and model matricies
   Precondition:
   Hx, Z are conformantly dimensioned (not checked)
 */
{
	Hx = Hx_init;
	Zv = Zv_init;
};

Simple_linear_correlated_observe_model::Simple_linear_correlated_observe_model (const FM::Matrix& Hx_init, const FM::SymMatrix& Z_init) :
	Linear_correlated_observe_model (Hx_init.size2(), Hx_init.size1())
/* Linear observe model initialised from model matricies
   Precondition:
   Hx, Z are conformantly dimensioned (not checked)
 */
{
	Hx = Hx_init;
	Z = Z_init;
};

Simple_linear_uncorrelated_observe_model::Simple_linear_uncorrelated_observe_model (const FM::Matrix& Hx_init, const FM::Vec& Zv_init) :
	Linear_uncorrelated_observe_model (Hx_init.size2(), Hx_init.size1())
/* Linear observe model initialised from model matricies
   Precondition:
   Hx, Z are conformantly dimensioned (not checked)
 */
{
	Hx = Hx_init;
	Zv = Zv_init;
};



Bayes_base::Float
 General_LzUnAd_observe_model::L(const FM::Vec& x) const
/*
 * Definition of likelihood given an addative Gaussian observation model:
 *  p(z|x) = exp(-0.5*(z-h(x))'*inv(Z)*(z-h(x))) / sqrt(2pi^nz*det(Z));
 *  L(x) the the Likelihood L(x) doesn't depend on / sqrt(2pi^nz) for constant z size
 * Precond: Observation Information: z,Zv_inv,detZterm
 */
{
	using namespace FM;
	if (!zset)
		Bayes_filter_exception ("General_observe_model used without Lz set");
					// Predict observations using supplied observe model
	const Vec& zp = h(x);
					// Normalise Innovation about z
	Vec zInnov(z.size());
	zInnov = zp;;
	normalise (zInnov, z);
	zInnov -= z;

	// Likelihood w of observation z given particlar state xi is true state
	// The state, xi, defines a predicted observation with a gaussian 
	// distribution with variance Zd. Thus, the likelihood can be determined directly from the gaussian
	Float logL = 0.;

	Vec::iterator zInnovi = zInnov.begin();
	for (Vec::const_iterator ZIi = Zv_inv.begin(); ZIi != Zv_inv.end(); ++ZIi) {
		logL += sqr(*zInnovi) * (*ZIi);
		++zInnovi;
	}
	using namespace std;
	return exp(-0.5*(logL + logdetZ));
}

void General_LzUnAd_observe_model::Lz (const FM::Vec& zz)
/* Set the observation zz and Zv about which to evaluate the Likelihood function
 * Postcond: Observation Information: z,Zv_inv,detZterm
 */
{
	z = zz;
	zset = true;
					// Compute inverse of Zv and its reciprocal condition number
	Float rcond = FM::UdUrcond_vec(Zv);
	rclimit.check_PD(rcond, "Z not PD in observe");

	Float detZ = 1.;
	FM::Vec::iterator ZIi = Zv_inv.begin();
	for (FM::Vec::const_iterator zi = Zv.begin(); zi != Zv.end(); ++zi) {
		detZ *= *zi;
		*ZIi = 1./ (*zi);		// Protected from /0 by rcond check
		++ZIi;
	}
	using namespace std;
	logdetZ = log(detZ);		// Protected from ln(0) by rcond check
}

Bayes_base::Float
 General_LzCoAd_observe_model::L(const FM::Vec& x) const
/*
 * Definition of likelihood given an addative Gaussian observation model:
 *  p(z|x) = exp(-0.5*(z-h(x))'*inv(Z)*(z-h(x))) / sqrt(2pi^nz*det(Z));
 *  L(x) the the Likelihood L(x) doesn't depend on / sqrt(2pi^nz) for constant z size
 * Precond: Observation Information: z,Z_inv,detZterm
 */
{
	if (!zset)
		Bayes_filter_exception ("General_observe_model used without Lz set");
						// Predict observations using supplied observe model
	const FM::Vec& zp = h(x);
						// Normalise Innovation about z
	FM::Vec zInnov(z.size());
	zInnov = zp;
	normalise (zInnov, z);
	zInnov -= z;

	Float logL = scaled_vector_square(zInnov, Z_inv);
	using namespace std;
	return exp(-0.5*(logL + logdetZ));
}

void General_LzCoAd_observe_model::Lz (const FM::Vec& zz)
/* Set the observation zz and Z about which to evaluate the Likelihood function
 * Postcond: Observation Information: z,Z_inv,detZterm
 */
{
	z = zz;
	zset = true;
						// Compute inverse of Z and its reciprocal condition number
	Float detZ;
	Float rcond = FM::UdUinversePD (Z_inv, detZ, Z);
	rclimit.check_PD(rcond, "Z not PD in observe");
	logdetZ = log(detZ);		// Protected from ln(0) by rcond check
}


Bayes_base::Float
 General_LzCoAd_observe_model::scaled_vector_square(const FM::Vec& v, const FM::SymMatrix& V)
/*
 * Compute covariance scaled square inner product of a Vector: v'*V*v
 */
{
	Float p = 0.;

	FM::Vec::const_iterator vi = v.begin();
	for (FM::Subscript i = 0; i < V.size1(); ++i)
	{
		p += *vi * FM::inner_prod(V[i], v);
		++vi;
	}
	return p;
}


Adapted_Correlated_addative_observe_model::Adapted_Correlated_addative_observe_model (Uncorrelated_addative_observe_model& adapt) :
		Correlated_addative_observe_model(adapt.Zv.size()),
		unc(adapt)
{
	Z.clear();
	for (FM::Subscript i = 0; i < unc.Zv.size(); ++i)
		Z(i,i) = unc.Zv[i];
}

Adapted_Linrz_correlated_observe_model::Adapted_Linrz_correlated_observe_model (Linrz_uncorrelated_observe_model& adapt) :
		Linrz_correlated_observe_model(adapt.Hx.size2(), adapt.Hx.size1()),
		unc(adapt)
{
	Hx = unc.Hx;
	Z.clear();
	for (FM::Subscript i = 0; i < unc.Zv.size(); ++i)
		Z(i,i) = unc.Zv[i];
}


Sample_filter::~Sample_filter()
/*
 * Default definition required for a pure virtual distructor
 * Should be in BayesFlt but cannot be defined if Matrix has
 * private distructor
 */
{
}

void Sample_filter::predict (Functional_predict_model& f)
/*
 * Predict samples forward
 *		Pre : S represent the prior distribution
 *		Post: S represent the predicted distribution
 */
{
						// Predict particles S using supplied predict model
	for (FM::ColMatrix::iterator si = S.begin(); si != S.end(); ++si) {
		(*si).assign (f.fx(*si));
	}
}


FM::Subscript Sample_filter::unique_samples () const
/*
 * Count the number of unique samples in S
 */
{
	typedef FM::ColMatrix::const_iterator Sref;
	// Provide a ordering on samples
	struct order {
		static bool less(const Sref a, const Sref b)
		{
			FM::ColMatrix::OneD::const_iterator sai = (*a).begin();
			FM::ColMatrix::OneD::const_iterator sbi = (*b).begin();
			do
			{
				if (*sai < *sbi)
					return true;
				else if (*sai > *sbi)
					return false;

				++sai; ++sbi;
			} while (sai != (*a).end());
			// Equal
			return false;
		}
	};

	// TODO: remove the temporary vector. Requires a uBlas column iterator
						// Sorted reference container
	typedef std::vector<Sref> SRContainer;
	SRContainer sortR(S.size2());

						// Reference each element in S
	{	Sref elem = S.begin();
		SRContainer::iterator ssi = sortR.begin();
		for (; ssi < sortR.end(); ++ssi)
		{
			*ssi = elem; ++elem;
		}
	}

	std::sort (sortR.begin(), sortR.end(), order::less);

						// Count element changes, precond: sortS not empty
	FM::Subscript u = 1;
	SRContainer::const_iterator ssi= sortR.begin();
	SRContainer::const_iterator ssp = ssi;
	++ssi;
	while (ssi < sortR.end())
	{
		if (order::less(*ssp, *ssi))
			++u;
		ssp = ssi;
		++ssi;
	}
	return u;
}

}//namespace
