/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 */

/*
 * SLAM : Simultaneous Locatization and Mapping
 *  Kalman filter representing representation of SLAM
 */

// Include Bayesian filtering library
#include "BayesFilter/bayesFlt.hpp"
#include "SLAM.hpp"
#include "kalmanSLAM.hpp"

namespace SLAM_filter
{

template <class Base>
inline void zero(FM::ublas::matrix_range<Base> A)
// Zero a matrix_rangec
{	// Note A cannot be a reference
	typedef typename Base::value_type Base_value_type;
	A .assign (FM::ublas::scalar_matrix<Base_value_type>(A.size1(),A.size2(), Base_value_type()) );
}

Kalman_SLAM::Kalman_SLAM( Bayesian_filter::Linrz_kalman_filter& location_filter, Full_filter& filter_generator ) :
	SLAM(),
	loc(location_filter),
	fgenerator(filter_generator), nL(location_filter.x.size())
{
		// generate full with only location states
	full = filter_generator.generate(nL);
		// place location in full
	full->init_kalman (loc.x, loc.X);
	nM = 0;
}

void Kalman_SLAM::predict( BF::Linrz_predict_model& lpred )
{
		// extract location part of full
	loc.x = full->x.sub_range(0,nL);
	loc.X = full->X.sub_matrix(0,nL,0,nL);
		// predict location, independant of map
	loc.init();
	loc.predict (lpred);
	loc.update();
		// return location to full
	full->x.sub_range(0,nL) = loc.x;
	full->X.sub_matrix(0,nL,0,nL) = loc.X;
	full->init();
}

void Kalman_SLAM::observe( unsigned feature, const Feature_observe& fom, const FM::Vec& z )
{
	// Assume features added sequenatialy
	if (feature >= nM) {
		error (BF::Logic_exception("Observe non existing feature"));
		return;
	}
	// TODO Implement nonlinear form
	// Create a augmented sparse observe model for full states
	BF::Linear_uncorrelated_observe_model fullm(full->x.size(), 1);
	fullm.Hx.clear();
	fullm.Hx.sub_matrix(0,nL, 0,nL) = fom.Hx.sub_matrix(0,nL, 0,nL);
	fullm.Hx(0,nL+feature) = fom.Hx(0,nL);
	fullm.Zv = fom.Zv;
	full->observe(fullm, z);
}

void Kalman_SLAM::observe_new( unsigned feature, const Feature_observe_inverse& fom, const FM::Vec& z )
// fom: must have a the special from required for SLAM::obeserve_new
{
	FM::Vec sz(nL+1);
	sz.sub_range(0,nL) = full->x.sub_range(0,nL);
	sz[nL] = z[0];
	FM::Vec t = fom.h(sz);

		// Make space in scheme for feature, requires the scheme can deal with resized state
	if (feature >= nM)
	{
		nM = feature+1;	
		Full_filter::Filter_type* nf = fgenerator.generate(nL+nM);
		nf->x.sub_range(0,full->x.size()) .assign (full->x);
		nf->X.sub_matrix(0,full->x.size(),0,full->x.size()) .assign (full->X);

		nf->x[nL+feature] = t[0];
		nf->X(nL+feature,nL+feature) = fom.Zv[0];
		nf->init ();
		fgenerator.dispose(full);
		full = nf;
	}
	else
	{
		full->x[nL+feature] = t[0];
		full->X(nL+feature,nL+feature) = fom.Zv[0];
			// ISSUE uBLAS has problems accessing the lower symmetry via a sub_matrix proxy, there two two parts seperately
		zero( full->X.sub_matrix(0,nL+feature, nL+feature,nL+feature+1) );
		zero( full->X.sub_matrix(nL+feature,nL+feature+1, nL+feature,full->X.size1()) );
		full->init ();
	}

	full->X(nL+feature,nL+feature) = fom.Zv[0];
}

void Kalman_SLAM::observe_new( unsigned feature, const FM::Vec& t, const FM::Vec& T )
{
		// Make space in scheme for feature, requires the scheme can deal with resized state
	if (feature >= nM)
	{
		nM = feature+1;
		Full_filter::Filter_type* nf = fgenerator.generate(nL+nM);
		nf->x.sub_range(0,full->x.size()) .assign (full->x);
		nf->X.sub_matrix(0,full->x.size(),0,full->x.size()) .assign (full->X);

		nf->x[nL+feature] = t[0];
		nf->X(nL+feature,nL+feature) = T[0];
		nf->init ();
		fgenerator.dispose(full);
		full = nf;
	}
	else
	{
		full->x[nL+feature] = t[0];
		full->X(nL+feature,nL+feature) = T[0];
		full->init ();
	}
}

void Kalman_SLAM::forget( unsigned feature, bool must_exist )
{
	full->x[nL+feature] = 0.;
			// ISSUE uBLAS has problems accessing the lower symmetry via a sub_matrix proxy, there two two parts seperately
	zero( full->X.sub_matrix(0,nL+feature, nL+feature,nL+feature+1) );
	zero( full->X.sub_matrix(nL+feature,nL+feature+1, nL+feature,full->X.size1()) );
	full->init();
}

void Kalman_SLAM::decorrelate( Bayesian_filter::Bayes_base::Float d )
// Reduce correlation by scaling cross-correlation terms
{
	size_t i,j;
	const size_t n = full->X.size1();
	for (i = 1; i < n; ++i)
	{
		FM::SymMatrix::Row Xi(full->X,i);
		for (j = 0; j < i; ++j)
		{
			Xi[j] *= d;
		}
		for (j = i+1; j < n; ++j)
		{
			Xi[j] *= d;
		}
	}
	full->init();
}

}//namespace SLAM
