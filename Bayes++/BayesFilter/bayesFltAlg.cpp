/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */
 
/*
 * Bayesian_filter Implemention:
 *  algorithm that require additional resources to implement
 */
#include "bayesFlt.hpp"
#include "models.hpp"
#include "matSup.hpp"
#include <vector>		// Only for unique_samples


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
	FM::Vec s = h.h(x);		// Observation model, s is predicted observation
	h.normalise(s, z);

	s = z - s;
	return observe_innovation (h, s);
}

Bayes_base::Float
 Extended_filter::observe (Linrz_uncorrelated_observe_model& h, const FM::Vec& z)
/*
 * Extended linrz uncorrelated observe, compute innovation for observe_innovation
 */
{
	FM::Vec s = h.h(x);		// Observation model, s is predicted observation
	h.normalise(s, z);

	s = z - s;
	return observe_innovation (h, s);
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

	Vec::iterator zi = zInnov.begin(), zi_end = zInnov.end();
	for (; zi != zi_end; ++zi) {
		*zi *= *zi;
	}
	Float logL = inner_prod(zInnov, Zv_inv);

	using namespace std;
	return exp(Float(-0.5)*(logL + logdetZ));
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
	
	for (FM::Vec::const_iterator zi = Zv.begin(), zi_end = Zv.end(); zi != zi_end; ++zi) {
		detZ *= *zi;
		Zv_inv[zi.index()] = 1 / (*zi);		// Protected from /0 by rcond check
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
	return exp(Float(-0.5)*(logL + logdetZ));
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
	return FM::inner_prod(v, FM::prod(V,v));
}


Adapted_Correlated_addative_observe_model::Adapted_Correlated_addative_observe_model (Uncorrelated_addative_observe_model& adapt) :
		Correlated_addative_observe_model(adapt.Zv.size()),
		unc(adapt)
{
	Z.clear();
	for (size_t i = 0; i < unc.Zv.size(); ++i)
		Z(i,i) = unc.Zv[i];
}

Adapted_Linrz_correlated_observe_model::Adapted_Linrz_correlated_observe_model (Linrz_uncorrelated_observe_model& adapt) :
		Linrz_correlated_observe_model(adapt.Hx.size2(), adapt.Hx.size1()),
		unc(adapt)
{
	Hx = unc.Hx;
	Z.clear();
	for (size_t i = 0; i < unc.Zv.size(); ++i)
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
	const size_t nSamples = S.size2();
	for (size_t i = 0; i != nSamples; ++i) {
		FM::ColMatrix::Column Si(S,i);
		Si.assign (f.fx(Si));
	}
}


namespace {
	// Column proxy so S can be sorted indirectly
	struct ColProxy
	{
		const FM::ColMatrix* cm;
		size_t col;
		ColProxy& operator=(ColProxy& a)
		{
			col = a.col;
			return a;
		}
		// Provide a ordering on columns
		static bool less(const ColProxy& a, const ColProxy& b)
		{
			FM::ColMatrix::const_iterator1 sai = a.cm->find_first1(1, 0,a.col);
			FM::ColMatrix::const_iterator1 sai_end = a.cm->find_last1(1, a.cm->size1(),a.col); 
			FM::ColMatrix::const_iterator1 sbi = b.cm->find_first1(1,0, b.col);
			while (sai != sai_end)
			{
				if (*sai < *sbi)
					return true;
				else if (*sai > *sbi)
					return false;

				++sai; ++sbi;
			} ;
			return false;		// Equal
		}
	};
}//namespace

size_t Sample_filter::unique_samples () const
/*
 * Count number of unique (unequal value) samples in S
 * Implementation requires std::sort on sample column references
 */
{
						// Temporary container to Reference each element in S
	typedef std::vector<ColProxy> SRContainer;
	SRContainer sortR(S.size2());
	size_t col_index = 0;
	for (SRContainer::iterator si = sortR.begin(); si != sortR.end(); ++si) {
		(*si).cm = &S; (*si).col = col_index++;
	}
						// Sort the column proxies
	std::sort (sortR.begin(), sortR.end(), ColProxy::less);

						// Count element changes, precond: sortS not empty
	size_t u = 1;
	SRContainer::const_iterator ssi = sortR.begin();
	SRContainer::const_iterator ssp = ssi;
	++ssi;
	while (ssi < sortR.end())
	{
		if (ColProxy::less(*ssp, *ssi))
			++u;
		ssp = ssi;
		++ssi;
	}
	return u;
}

}//namespace
