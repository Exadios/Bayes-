/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2004 Michael Stevens
 * See accompanying Bayes++.html for terms and conditions of use.
 *
 * $Header$
 */

/*
 * SLAM : Simultaneous Locatization and Mapping
 *  Kalman filter representing representation of SLAM
 */

		// Bayes++ Bayesian filtering schemes
#include "BayesFilter/bayesFlt.hpp"
#include <iostream>
#include "boost/numeric/ublas/io.hpp"
		// Bayes++ SLAM
#include "SLAM.hpp"
#include "kalmanSLAM.hpp"

namespace SLAM_filter
{

template <class Base>
inline void zero(FM::ublas::matrix_range<Base> A)
// Zero a matrix_range
{	// Note A cannot be a reference
	typedef typename Base::value_type Base_value_type;
	FM::noalias(A) = FM::ublas::scalar_matrix<Base_value_type>(A.size1(),A.size2(), Base_value_type());
}

Kalman_SLAM::Kalman_SLAM( Kalman_filter_generator& filter_generator ) :
	SLAM(),
	fgenerator(filter_generator),
	full(NULL)
{
	nL = 0;
	nM = 0;
}

Kalman_SLAM::~Kalman_SLAM()
{
	fgenerator.dispose (full);
}

void Kalman_SLAM::init_kalman (const FM::Vec& x, const FM::SymMatrix& X)
{
	// TODO maintain map states
	nL = x.size();
	nM = 0;
	if (full) fgenerator.dispose (full);
		// generate full filter
	full = fgenerator.generate(nL);
		// initialise location states
	full->x.sub_range(0,nL) = x;
	full->X.sub_matrix(0,nL,0,nL) = X;
	full->init();
}

void Kalman_SLAM::predict( BF::Linear_predict_model& lpred )
{
		// build full predict model by augmenting location predict with map identity predict model
	BF::Linear_predict_model fullpred (full->x.size(), lpred.q.size());
	FM::identity (fullpred.Fx);
	fullpred.Fx.sub_matrix(0,lpred.Fx.size1(),0,lpred.Fx.size2()) = lpred.Fx;
	fullpred.G.clear ();
	fullpred.G.sub_matrix(0,lpred.G.size1(), 0,lpred.G.size2()) = lpred.G;
	fullpred.q = lpred.q;
		// predict 
	full->predict (fullpred);
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
	FM::Vec sz(nL+z.size());
	sz.sub_range(0,nL) = full->x.sub_range(0,nL);
	sz.sub_range(nL,nL+z.size() )= z;
	FM::Vec t = fom.h(sz);
	FM::SymMatrix Xl (full->X.sub_matrix(0,nL, 0,nL) );
	FM::Matrix Ha (fom.Hx.sub_matrix(0,1, 0,nL) );
	FM::Matrix Hb (fom.Hx.sub_matrix(0,1, nL,nL+z.size()) );
	FM::Matrix temp (1,nL);

		// Make space in scheme for feature, requires the scheme can deal with resized state
	if (feature >= nM)
	{
		nM = feature+1;	
		Kalman_filter_generator::Filter_type* nf = fgenerator.generate(nL+nM);
		FM::noalias(nf->x.sub_range(0,full->x.size())) = full->x;
		FM::noalias(nf->X.sub_matrix(0,full->x.size(),0,full->x.size())) = full->X;

		fgenerator.dispose(full);
		full = nf;
	}
std::cout << full->x <<std::endl;
std::cout << full->X <<std::endl;
	full->x[nL+feature] = t[0];
std::cout << full->x <<std::endl;
	full->X(nL+feature,nL+feature) = ( FM::prod_SPD(Ha,Xl,temp) + FM::prod_SPD(Hb,fom.Zv) )(0,0);
std::cout << full->X <<std::endl;
	//	full->X.sub_matrix(nL+feature,nL+feature+1,0,nL) = prod(Ha,Xl);
	{	// Old uBLAS has problems assigning to symmetric proxy as above, go element by element
		const FM::Matrix cross (prod(Ha,Xl) );
		const FM::Float definate_epsilon = 0.9999;
		for (size_t i = 0; i != nL; ++i)
			full->X(nL+feature, i) = definate_epsilon * cross(0,i);
	}
std::cout << full->X <<std::endl;
		
	full->init ();
}

void Kalman_SLAM::observe_new( unsigned feature, const FM::Float& t, const FM::Float& T )
{
		// Make space in scheme for feature, requires the scheme can deal with resized state
	if (feature >= nM)
	{
		nM = feature+1;
		Kalman_filter_generator::Filter_type* nf = fgenerator.generate(nL+nM);
		FM::noalias(nf->x.sub_range(0,full->x.size())) = full->x;
		FM::noalias(nf->X.sub_matrix(0,full->x.size(),0,full->x.size())) = full->X;

		nf->x[nL+feature] = t;
		nf->X(nL+feature,nL+feature) = T;
		nf->init ();
		fgenerator.dispose(full);
		full = nf;
	}
	else
	{
		full->x[nL+feature] = t;
		full->X(nL+feature,nL+feature) = T;
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
