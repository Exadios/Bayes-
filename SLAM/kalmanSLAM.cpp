/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens, Australian Centre for Field Robotics
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 */

/*
 * SLAM : Simultaneous Locatization and Mapping
 *  FastSLAM augmented particle algorithm
 */

// Include Bayesian filtering library
#include "BayesFilter/bayesFlt.hpp"
#include "SLAM.hpp"
#include "kalmanSLAM.hpp"


namespace SLAM_filter
{

template <class Base>
inline void clear(FM::ublas::matrix_range<Base> A)
// Clear a matrix_range as clear not supported for uBLAS proxies
{	// Note A cannot be a reference
	typedef typename Base::value_type Base_value_type;
	A.assign(FM::ublas::scalar_matrix<Base_value_type>(A.size1(),A.size2(), Base_value_type()) );
}

Kalman_SLAM::Kalman_SLAM( Bayesian_filter::Linrz_kalman_filter& full_filter, unsigned full_nL ) :
	SLAM(), full(full_filter), nL(full_nL)
{
	nM = 0;
}

void Kalman_SLAM::predict (BF::Linear_predict_model& lpred)
{
	// lpred is only prediction model for Location. Create a augmented sparse model for full states
	BF::Linear_predict_model all(full.x.size(),full.x.size());
	FM::identity(all.Fx);
	all.G.clear();
	all.q.clear();
	all.Fx.sub_matrix(0,nL, 0,nL) = lpred.Fx;
	all.G.sub_matrix(0,nL, 0,lpred.G.size2()) = lpred.G;
	all.q(0,lpred.q.size()) = lpred.q;

	full.predict(all);
}

void Kalman_SLAM::observe (unsigned feature, const Feature_observe& fom, const FM::Vec& z)
{
	// Assume features added sequenatialy
	if (feature >= nM) {
		error (BF::Logic_exception("Observe non existing feature"));
		return;
	}
	// TODO Implement nonlinear form
	// Create a augmented sparse observe model for full states
	BF::Linear_uncorrelated_observe_model fullm(nL+nM, 1);
	fullm.Hx.clear();
	fullm.Hx.sub_matrix(0,nL, 0,nL) = fom.Hx.sub_matrix(0,nL, 0,nL);
	fullm.Hx(0,nL+feature) = fom.Hx(0,nL);
	fullm.Zv = fom.Zv;
	full.observe(fullm, z);
}

void Kalman_SLAM::observe_new (unsigned feature, const Feature_observe_inverse& fom, const FM::Vec& z)
// fom: must have a the special from required for SLAM::obeserve_new
{
	if (nL+feature >= full.x.size()) {
		error (BF::Logic_exception("Observe_new no feature space"));
		return;
	}
	++nM;

	FM::Vec sz(nL+1);
	sz(0,nL) = full.x(0,nL);
	sz[nL] = z[0];
	FM::Vec t = fom.h(sz);

	full.x[nL+feature] = t[0];

	clear( full.X.sub_matrix(0,full.X.size1(), nL+feature,nL+feature+1) );
	clear( full.X.sub_matrix(nL+feature,nL+feature+1, 0,full.X.size1()) );
	full.X(nL+feature,nL+feature) = fom.Zv[0];
}

void Kalman_SLAM::observe_new( unsigned feature, const FM::Vec& t, const FM::Vec& T)
{
	if (nL+feature >= full.x.size()) {
		error (BF::Logic_exception("Observe_new no feature space"));
		return;
	}
	++nM;

	full.x[nL+feature] = t[0];
	full.X(nL+feature,nL+feature) = T[0];
}

void Kalman_SLAM::forget( unsigned feature, bool must_exist)
{
	full.x[feature] = 0.;
	clear( full.X.sub_matrix(0,full.X.size1(), nL+feature,nL+feature+1) );
	clear( full.X.sub_matrix(nL+feature,nL+feature+1, 0,full.X.size1()) );
	full.init();
	--nM;
}

void Kalman_SLAM::decorrelate( Bayesian_filter::Bayes_base::Float d)
// Reduce correlation by scaling cross-correlation terms
{
	size_t i,j;
	const size_t n = full.X.size1();
	for (i = 1; i < n; ++i)
	{
		FM::SymMatrix::Row Xi(full.X,i);
		for (j = 0; j < i; ++j)
		{
			Xi[j] *= d;
		}
		for (j = i+1; j < n; ++j)
		{
			Xi[j] *= d;
		}
	}
}

}//namespace SLAM
